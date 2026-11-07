/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file communicator.cpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2026
 */

#include <algorithm>
#include <cstdint>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#define HMAP_HAS_WIN32_FIBERS 1
#define HMAP_HAS_UCONTEXT 0
#elif defined(__linux__) || defined(__APPLE__)
#include <ucontext.h>
#define HMAP_HAS_WIN32_FIBERS 0
#define HMAP_HAS_UCONTEXT 1
#else
#define HMAP_HAS_WIN32_FIBERS 0
#define HMAP_HAS_UCONTEXT 0
#endif

#include "highmap/array.hpp"
#include "highmap/math/core.hpp"
#include "highmap/virtual_array/communicator.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

// =====================================
// Free Functions
// =====================================

std::string to_string(ReduceOp op)
{
  switch (op)
  {
  case ReduceOp::Min: return "Min";
  case ReduceOp::Max: return "Max";
  case ReduceOp::Sum: return "Sum";
  }
  return "Unknown";
}

// =====================================
// SingleArrayComm Implementation
// =====================================

float SingleArrayComm::all_reduce(float local_val, ReduceOp /*op*/)
{
  return local_val;
}

std::vector<float> SingleArrayComm::all_union(
    const std::vector<float> &local_vals)
{
  return local_vals;
}

void SingleArrayComm::barrier()
{
}

int SingleArrayComm::rank() const
{
  return 0;
}

int SingleArrayComm::size() const
{
  return 1;
}

void SingleArrayComm::sync_overlap_buffers(Array & /*local_tile*/,
                                           const TileRegion & /*region*/,
                                           SyncOperation /*op*/)
{
  // Standalone single array has no neighboring tiles or overlap seams
}

SingleArrayComm &SingleArrayComm::instance()
{
  static SingleArrayComm inst;
  return inst;
}

// =====================================
// SequentialTileComm Implementation
// =====================================

namespace
{

#if HMAP_HAS_UCONTEXT || HMAP_HAS_WIN32_FIBERS
struct TileFiber
{
#if HMAP_HAS_WIN32_FIBERS
  LPVOID handle = nullptr;
#else
  ucontext_t        ctx{};
  std::vector<char> stack;
#endif
  TileRegion region;
  bool       started = false;
  bool       finished = false;
};

struct FiberSchedulerState
{
#if HMAP_HAS_WIN32_FIBERS
  LPVOID main_fiber = nullptr;
#else
  ucontext_t main_ctx{};
#endif
  std::vector<TileFiber> fibers;
  int                    active_tile_idx = 0;
  SequentialTileComm    *p_comm = nullptr;
  const std::function<void(const TileRegion &, ICommunicator &)>
      *p_dispatcher = nullptr;
};

#if HMAP_HAS_WIN32_FIBERS
void WINAPI fiber_trampoline_win32(LPVOID lp_param)
{
  auto *sched = static_cast<FiberSchedulerState *>(lp_param);
  int   idx = sched->active_tile_idx;

  (*sched->p_dispatcher)(sched->fibers[idx].region, *sched->p_comm);

  sched->fibers[idx].finished = true;
  SwitchToFiber(sched->main_fiber);
}
#else
void fiber_trampoline(uint32_t low, uint32_t high)
{
  uintptr_t ptr = (static_cast<uintptr_t>(high) << 32) |
                  static_cast<uintptr_t>(low);
  auto *sched = reinterpret_cast<FiberSchedulerState *>(ptr);
  int   idx = sched->active_tile_idx;

  (*sched->p_dispatcher)(sched->fibers[idx].region, *sched->p_comm);

  sched->fibers[idx].finished = true;
}
#endif
#endif

} // namespace

SequentialTileComm::SequentialTileComm(int total_tiles, int halo)
    : total_tiles(total_tiles), halo(halo)
{
}

SequentialTileComm::SequentialTileComm(const VirtualArray          &ref_va,
                                       const std::function<void()> &on_suspend,
                                       const std::function<void()> &on_resume)
    : total_tiles(ref_va.get_ntiles()),
      halo(ref_va.halo),
      p_ref_va(&ref_va),
      on_suspend(on_suspend),
      on_resume(on_resume)
{
}

void SequentialTileComm::yield_to_scheduler()
{
#if HMAP_HAS_UCONTEXT || HMAP_HAS_WIN32_FIBERS
  if (!this->p_fiber_scheduler) return;

  if (this->on_suspend) this->on_suspend();

  auto *sched = static_cast<FiberSchedulerState *>(this->p_fiber_scheduler);
  int   idx = this->current_tile_idx;

#if HMAP_HAS_WIN32_FIBERS
  SwitchToFiber(sched->main_fiber);
#else
  swapcontext(&sched->fibers[idx].ctx, &sched->main_ctx);
#endif
#endif
}

float SequentialTileComm::all_reduce(float local_val, ReduceOp op)
{
  size_t seq = this->call_sequence++;

  if (seq >= this->accumulators.size())
  {
    this->accumulators.resize(seq + 1);
    this->initialized.resize(seq + 1, false);
    this->reduce_ops.resize(seq + 1, op);
  }

  if (!this->initialized[seq])
  {
    this->accumulators[seq] = local_val;
    this->initialized[seq] = true;
    this->reduce_ops[seq] = op;
  }
  else
  {
    switch (op)
    {
    case ReduceOp::Min:
      this->accumulators[seq] = std::min(this->accumulators[seq], local_val);
      break;
    case ReduceOp::Max:
      this->accumulators[seq] = std::max(this->accumulators[seq], local_val);
      break;
    case ReduceOp::Sum: this->accumulators[seq] += local_val; break;
    }
  }

  if (this->p_fiber_scheduler && this->total_tiles > 1)
  {
    this->yield_to_scheduler();
    if (this->on_resume) this->on_resume();
  }

  return this->accumulators[seq];
}

std::vector<float> SequentialTileComm::all_union(
    const std::vector<float> &local_vals)
{
  size_t seq = this->union_sequence++;

  if (seq >= this->union_accumulators.size())
  {
    this->union_accumulators.resize(seq + 1);
  }

  auto &acc = this->union_accumulators[seq];
  for (float v : local_vals)
  {
    acc.push_back(v);
  }

  std::sort(acc.begin(), acc.end());
  auto last = std::unique(acc.begin(), acc.end());
  acc.erase(last, acc.end());

  if (this->p_fiber_scheduler && this->total_tiles > 1)
  {
    this->yield_to_scheduler();
    if (this->on_resume) this->on_resume();
  }

  return this->union_accumulators[seq];
}

void SequentialTileComm::barrier()
{
  if (this->p_fiber_scheduler && this->total_tiles > 1)
  {
    this->yield_to_scheduler();
    if (this->on_resume) this->on_resume();
  }
}

float SequentialTileComm::get_accumulated_value(size_t seq) const
{
  if (seq < this->accumulators.size()) return this->accumulators[seq];
  return 0.f;
}

void SequentialTileComm::next_tile(int tile_idx)
{
  this->current_tile_idx = tile_idx;
  this->call_sequence = 0;
  this->sync_sequence = 0;
  this->union_sequence = 0;
}

int SequentialTileComm::rank() const
{
  return this->current_tile_idx;
}

int SequentialTileComm::size() const
{
  return this->total_tiles;
}

void SequentialTileComm::sync_overlap_buffers(Array            &local_tile,
                                              const TileRegion &region,
                                              SyncOperation     op)
{
  if (this->halo <= 0) return;

  // --- Option A2: Cooperative fiber suspension + disk-backed VirtualArray
  // (VA_DISK_LRU_MIN, max 2 live tiles in memory)
  if (this->p_fiber_scheduler && this->p_ref_va && this->total_tiles > 1)
  {
    if (!this->sync_va)
    {
      this->sync_va = std::make_unique<VirtualArray>(
          this->p_ref_va->shape,
          this->p_ref_va->bbox,
          this->p_ref_va->tile_shape,
          this->p_ref_va->halo,
          StorageMode::VA_DISK_LRU_MIN);
    }

    // 1. Write current tile state to disk-backed sync_va and release
    const glm::ivec2 saved_shape = local_tile.shape;
    {
      Array &staging = this->sync_va->storage->get_tile(region);
      staging = local_tile;
      std::vector<float>().swap(local_tile.vector);
      this->sync_va->storage->release_tile(region);
      this->sync_va->storage->trim();
    }

    this->pending_sync_op = op;
    this->has_pending_sync = true;

    // 2. Suspend this tile fiber until all tiles reach this sync point and
    // the scheduler executes sync_va->sync_overlap_buffers(op)
    this->yield_to_scheduler();

    // 3. Re-acquire input/output tiles FIRST, then reload local_tile so that
    // if local_tile IS one of the input/output tiles in VA_DISK_LRU_MIN,
    // local_tile gets the synchronized data instead of being overwritten by
    // the stale disk copy
    if (this->on_resume) this->on_resume();

    {
      Array &staging = this->sync_va->storage->get_tile(region);
      local_tile.shape = saved_shape;
      local_tile = staging;
      this->sync_va->storage->release_tile(region);
      this->sync_va->storage->trim();
    }

    return;
  }

  // --- Non-cooperative fallback
  const size_t seq = this->sync_sequence++;

  if (seq >= this->east_edges.size())
  {
    this->east_edges.resize(seq + 1);
  }
  if (seq >= this->south_edges.size())
  {
    this->south_edges.resize(seq + 1);
  }

  auto &seq_east_edges = this->east_edges[seq];
  auto &seq_south_edges = this->south_edges[seq];

  // In sequential traversal (raster order tx then ty), a tile can receive
  // incoming overlap data from its west neighbor (tx - 1) and north neighbor
  // (ty - 1), and send its boundary overlap data to its east neighbor (tx + 1)
  // and south neighbor (ty + 1).

  const int tx = region.key.tx;
  const int ty = region.key.ty;

  // 1. Ingest Western overlap from (tx - 1, ty) if present
  if (tx > 0 && !seq_east_edges.empty())
  {
    size_t prev_tile_key = static_cast<size_t>(ty * 10000 + (tx - 1));
    if (prev_tile_key < seq_east_edges.size() &&
        !seq_east_edges[prev_tile_key].empty())
    {
      const auto &src_col = seq_east_edges[prev_tile_key];
      // Western halo region in local_tile: x in [0, halo)
      for (int p = 0; p < this->halo; ++p)
      {
        for (int q = 0; q < local_tile.shape.y; ++q)
        {
          size_t src_idx = p * local_tile.shape.y + q;
          if (src_idx < src_col.size())
          {
            float v_prev = src_col[src_idx];
            float v_curr = local_tile(p, q);

            switch (op)
            {
            case SyncOperation::Average:
            case SyncOperation::Mean:
              local_tile(p, q) = 0.5f * (v_prev + v_curr);
              break;
            case SyncOperation::CopyFirst: local_tile(p, q) = v_prev; break;
            case SyncOperation::CopySecond: break; // keep local_tile
            case SyncOperation::Min:
              local_tile(p, q) = std::min(v_prev, v_curr);
              break;
            case SyncOperation::Max:
              local_tile(p, q) = std::max(v_prev, v_curr);
              break;
            case SyncOperation::SmoothBlend:
            {
              float r = (this->halo > 1) ? float(p) / float(this->halo - 1)
                                         : 0.5f;
              local_tile(p, q) = lerp(v_prev, v_curr, smoothstep5(r));
            }
            break;
            }
          }
        }
      }
    }
  }

  // 2. Ingest Northern overlap from (tx, ty - 1) if present
  if (ty > 0 && !seq_south_edges.empty())
  {
    size_t prev_tile_key = static_cast<size_t>((ty - 1) * 10000 + tx);
    if (prev_tile_key < seq_south_edges.size() &&
        !seq_south_edges[prev_tile_key].empty())
    {
      const auto &src_row = seq_south_edges[prev_tile_key];
      // Northern halo region in local_tile: y in [0, halo)
      for (int q = 0; q < this->halo; ++q)
      {
        for (int p = 0; p < local_tile.shape.x; ++p)
        {
          size_t src_idx = q * local_tile.shape.x + p;
          if (src_idx < src_row.size())
          {
            float v_prev = src_row[src_idx];
            float v_curr = local_tile(p, q);

            switch (op)
            {
            case SyncOperation::Average:
            case SyncOperation::Mean:
              local_tile(p, q) = 0.5f * (v_prev + v_curr);
              break;
            case SyncOperation::CopyFirst: local_tile(p, q) = v_prev; break;
            case SyncOperation::CopySecond: break;
            case SyncOperation::Min:
              local_tile(p, q) = std::min(v_prev, v_curr);
              break;
            case SyncOperation::Max:
              local_tile(p, q) = std::max(v_prev, v_curr);
              break;
            case SyncOperation::SmoothBlend:
            {
              float r = (this->halo > 1) ? float(q) / float(this->halo - 1)
                                         : 0.5f;
              local_tile(p, q) = lerp(v_prev, v_curr, smoothstep5(r));
            }
            break;
            }
          }
        }
      }
    }
  }

  // 3. Stash Eastern boundary cells to pass to (tx + 1, ty)
  size_t curr_key = static_cast<size_t>(ty * 10000 + tx);
  if (curr_key >= seq_east_edges.size())
  {
    seq_east_edges.resize(curr_key + 1);
  }
  auto &east_buf = seq_east_edges[curr_key];
  east_buf.resize(this->halo * local_tile.shape.y);

  for (int p = 0; p < this->halo; ++p)
  {
    int pbuf = local_tile.shape.x - 2 * this->halo + p;
    for (int q = 0; q < local_tile.shape.y; ++q)
    {
      east_buf[p * local_tile.shape.y + q] = local_tile(pbuf, q);
    }
  }

  // 4. Stash Southern boundary cells to pass to (tx, ty + 1)
  if (curr_key >= seq_south_edges.size())
  {
    seq_south_edges.resize(curr_key + 1);
  }
  auto &south_buf = seq_south_edges[curr_key];
  south_buf.resize(this->halo * local_tile.shape.x);

  for (int q = 0; q < this->halo; ++q)
  {
    int qbuf = local_tile.shape.y - 2 * this->halo + q;
    for (int p = 0; p < local_tile.shape.x; ++p)
    {
      south_buf[q * local_tile.shape.x + p] = local_tile(p, qbuf);
    }
  }
}

void SequentialTileComm::run_cooperative_tiles(
    const std::function<void(const TileRegion &, ICommunicator &)> &dispatcher)
{
  if (!this->p_ref_va) return;

  const VirtualArray &ref_va = *this->p_ref_va;
  const int           nx = ceil_div(ref_va.shape.x, ref_va.tile_shape.x);
  const int           ny = ceil_div(ref_va.shape.y, ref_va.tile_shape.y);
  const int           n_tiles = nx * ny;

  if (n_tiles <= 0) return;

#if HMAP_HAS_UCONTEXT || HMAP_HAS_WIN32_FIBERS
  if (n_tiles == 1)
  {
    this->next_tile(0);
    TileRegion region = ref_va.tile_region_from_tile_coords(0, 0);
    dispatcher(region, *this);
    return;
  }

  FiberSchedulerState sched;
  sched.p_comm = this;
  sched.p_dispatcher = &dispatcher;
  sched.fibers.resize(n_tiles);

#if HMAP_HAS_WIN32_FIBERS
  bool converted_thread_to_fiber = false;
  if (IsThreadAFiber())
  {
    sched.main_fiber = GetCurrentFiber();
  }
  else
  {
    sched.main_fiber = ConvertThreadToFiber(nullptr);
    converted_thread_to_fiber = (sched.main_fiber != nullptr);
  }
#endif

  int idx = 0;
  for (int ty = 0; ty < ny; ++ty)
    for (int tx = 0; tx < nx; ++tx)
    {
      sched.fibers[idx].region = ref_va.tile_region_from_tile_coords(tx, ty);
      ++idx;
    }

  this->p_fiber_scheduler = &sched;

#if !HMAP_HAS_WIN32_FIBERS
  const uintptr_t ptr_val = reinterpret_cast<uintptr_t>(&sched);
  const uint32_t  ptr_low = static_cast<uint32_t>(ptr_val & 0xFFFFFFFFULL);
  const uint32_t  ptr_high = static_cast<uint32_t>((ptr_val >> 32) &
                                                  0xFFFFFFFFULL);
#endif

  constexpr size_t    k_fiber_stack_bytes = 256 * 1024;
  std::vector<size_t> saved_call_sequences(n_tiles, 0);
  std::vector<size_t> saved_union_sequences(n_tiles, 0);

  while (true)
  {
    bool any_active = false;
    this->has_pending_sync = false;

    for (int t = 0; t < n_tiles; ++t)
    {
      auto &fiber = sched.fibers[t];
      if (fiber.finished) continue;

      any_active = true;
      sched.active_tile_idx = t;
      this->current_tile_idx = t;
      this->call_sequence = saved_call_sequences[t];
      this->union_sequence = saved_union_sequences[t];

      if (!fiber.started)
      {
        fiber.started = true;
#if HMAP_HAS_WIN32_FIBERS
        fiber.handle = CreateFiber(k_fiber_stack_bytes,
                                   fiber_trampoline_win32,
                                   &sched);
#else
        fiber.stack.resize(k_fiber_stack_bytes);

        getcontext(&fiber.ctx);
        fiber.ctx.uc_stack.ss_sp = fiber.stack.data();
        fiber.ctx.uc_stack.ss_size = fiber.stack.size();
        fiber.ctx.uc_link = &sched.main_ctx;

        makecontext(&fiber.ctx,
                    reinterpret_cast<void (*)()>(fiber_trampoline),
                    2,
                    ptr_low,
                    ptr_high);
#endif
      }

#if HMAP_HAS_WIN32_FIBERS
      SwitchToFiber(fiber.handle);
#else
      swapcontext(&sched.main_ctx, &fiber.ctx);
#endif

      saved_call_sequences[t] = this->call_sequence;
      saved_union_sequences[t] = this->union_sequence;

      if (fiber.finished)
      {
#if HMAP_HAS_WIN32_FIBERS
        if (fiber.handle)
        {
          DeleteFiber(fiber.handle);
          fiber.handle = nullptr;
        }
#else
        std::vector<char>().swap(fiber.stack);
#endif
      }
    }

    if (!any_active) break;

    if (this->has_pending_sync && this->sync_va)
    {
      // Stitch all 4 directions across the VirtualArray with at most 2 live
      // tiles in memory at any point
      this->sync_va->sync_overlap_buffers(this->pending_sync_op);
      this->sync_va->trim_storage();
      this->has_pending_sync = false;
    }
  }

#if HMAP_HAS_WIN32_FIBERS
  if (converted_thread_to_fiber)
  {
    ConvertFiberToThread();
  }
#endif

  this->p_fiber_scheduler = nullptr;
  this->sync_va.reset();
#else
  int tile_idx = 0;
  for (int ty = 0; ty < ny; ++ty)
    for (int tx = 0; tx < nx; ++tx)
    {
      this->next_tile(tile_idx++);
      TileRegion region = ref_va.tile_region_from_tile_coords(tx, ty);
      dispatcher(region, *this);
    }
#endif
}

// =====================================
// ThreadedCollectiveState Implementation
// =====================================

ThreadedCollectiveState::ThreadedCollectiveState(int num_threads)
    : total_threads(num_threads)
{
}

float ThreadedCollectiveState::all_reduce(size_t seq, float val, ReduceOp op)
{
  std::unique_lock<std::mutex> lock(this->mtx);

  auto &slot = this->reduce_slots[seq];
  if (slot.remaining == 0)
  {
    slot.remaining = this->total_threads;
    slot.value = val;
  }
  else
  {
    switch (op)
    {
    case ReduceOp::Min: slot.value = std::min(slot.value, val); break;
    case ReduceOp::Max: slot.value = std::max(slot.value, val); break;
    case ReduceOp::Sum: slot.value += val; break;
    }
  }

  if (--slot.remaining == 0)
  {
    slot.ready = true;
    this->cv.notify_all();
  }
  else
  {
    this->cv.wait(lock, [&slot]() { return slot.ready; });
  }

  return slot.value;
}

std::vector<float> ThreadedCollectiveState::all_union(
    size_t                    seq,
    const std::vector<float> &vals)
{
  std::unique_lock<std::mutex> lock(this->mtx);

  auto &slot = this->union_slots[seq];
  if (slot.remaining == 0)
  {
    slot.remaining = this->total_threads;
  }

  for (float v : vals)
  {
    slot.values.push_back(v);
  }

  if (--slot.remaining == 0)
  {
    std::sort(slot.values.begin(), slot.values.end());
    auto last = std::unique(slot.values.begin(), slot.values.end());
    slot.values.erase(last, slot.values.end());

    slot.ready = true;
    this->cv.notify_all();
  }
  else
  {
    this->cv.wait(lock, [&slot]() { return slot.ready; });
  }

  return slot.values;
}

void ThreadedCollectiveState::barrier(size_t seq)
{
  std::unique_lock<std::mutex> lock(this->mtx);

  auto &slot = this->barrier_slots[seq];
  if (slot.remaining == 0)
  {
    slot.remaining = this->total_threads;
  }

  if (--slot.remaining == 0)
  {
    slot.ready = true;
    this->cv.notify_all();
  }
  else
  {
    this->cv.wait(lock, [&slot]() { return slot.ready; });
  }
}

int ThreadedCollectiveState::get_total_threads() const
{
  return this->total_threads;
}

void ThreadedCollectiveState::register_tile(int               tile_idx,
                                            Array            *tile,
                                            const TileRegion &region)
{
  std::unique_lock<std::mutex> lock(this->mtx);
  this->registered_tiles[tile_idx] = tile;
  this->registered_regions.insert_or_assign(tile_idx, region);
}

Array *ThreadedCollectiveState::get_registered_tile(int tile_idx) const
{
  auto it = this->registered_tiles.find(tile_idx);
  if (it != this->registered_tiles.end())
  {
    return it->second;
  }
  return nullptr;
}

const TileRegion *ThreadedCollectiveState::get_registered_region(
    int tile_idx) const
{
  auto it = this->registered_regions.find(tile_idx);
  if (it != this->registered_regions.end())
  {
    return &it->second;
  }
  return nullptr;
}

// =====================================
// ThreadedTileComm Implementation
// =====================================

ThreadedTileComm::ThreadedTileComm(
    std::shared_ptr<ThreadedCollectiveState> state,
    int                                      thread_rank,
    int                                      total_size,
    int                                      grid_nx)
    : state(state),
      thread_rank(thread_rank),
      total_size(total_size),
      grid_nx(grid_nx)
{
}

float ThreadedTileComm::all_reduce(float local_val, ReduceOp op)
{
  if (this->state)
  {
    return this->state->all_reduce(this->reduce_sequence++, local_val, op);
  }
  return local_val;
}

std::vector<float> ThreadedTileComm::all_union(
    const std::vector<float> &local_vals)
{
  if (this->state)
  {
    return this->state->all_union(this->union_sequence++, local_vals);
  }
  return local_vals;
}

void ThreadedTileComm::barrier()
{
  if (this->state)
  {
    this->state->barrier(this->barrier_sequence++);
  }
}

int ThreadedTileComm::rank() const
{
  return this->thread_rank;
}

int ThreadedTileComm::size() const
{
  return this->total_size;
}

void ThreadedTileComm::sync_overlap_buffers(Array            &local_tile,
                                            const TileRegion &region,
                                            SyncOperation     op)
{
  if (!this->state) return;

  const int halo = std::max(
      {region.halo.x, region.halo.y, region.halo.z, region.halo.w});
  if (halo <= 0) return;

  // Step 1: Every rank registers its local tile and region
  this->state->register_tile(this->thread_rank, &local_tile, region);

  // Step 2: Synchronize so all tiles are registered and ready
  this->barrier();

  // Step 3: Rank #0 performs the global boundary and corner stitching pass
  if (this->thread_rank == 0)
  {
    const int nx = this->grid_nx > 0 ? this->grid_nx : 1;
    const int ny = ceil_div(this->total_size, nx);

    // --- Pass 1: X-direction stitching
    for (int ty = 0; ty < ny; ++ty)
    {
      for (int tx = 0; tx < nx - 1; ++tx)
      {
        int idx0 = ty * nx + tx;
        int idx1 = ty * nx + (tx + 1);

        Array *t0 = this->state->get_registered_tile(idx0);
        Array *t1 = this->state->get_registered_tile(idx1);
        if (!t0 || !t1) continue;

        Array &tile0 = *t0;
        Array &tile1 = *t1;

        switch (op)
        {
        case SyncOperation::Average:
        case SyncOperation::Mean:
          for (int p = 0; p < halo; ++p)
          {
            int pbuf = tile0.shape.x - 2 * halo + p;
            for (int q = 0; q < tile0.shape.y; ++q)
            {
              float avg = 0.5f * (tile0(pbuf, q) + tile1(p, q));
              tile0(pbuf, q) = avg;
              tile1(p, q) = avg;
            }
          }
          break;

        case SyncOperation::CopyFirst:
          for (int p = 0; p < halo; ++p)
          {
            int pbuf = tile0.shape.x - 2 * halo + p;
            for (int q = 0; q < tile0.shape.y; ++q)
            {
              tile1(p, q) = tile0(pbuf, q);
            }
          }
          break;

        case SyncOperation::CopySecond:
          for (int p = 0; p < halo; ++p)
          {
            int pbuf = tile0.shape.x - 2 * halo + p;
            for (int q = 0; q < tile0.shape.y; ++q)
            {
              tile0(pbuf, q) = tile1(p, q);
            }
          }
          break;

        case SyncOperation::Min:
          for (int p = 0; p < halo; ++p)
          {
            int pbuf = tile0.shape.x - 2 * halo + p;
            for (int q = 0; q < tile0.shape.y; ++q)
            {
              float m = std::min(tile0(pbuf, q), tile1(p, q));
              tile0(pbuf, q) = m;
              tile1(p, q) = m;
            }
          }
          break;

        case SyncOperation::Max:
          for (int p = 0; p < halo; ++p)
          {
            int pbuf = tile0.shape.x - 2 * halo + p;
            for (int q = 0; q < tile0.shape.y; ++q)
            {
              float m = std::max(tile0(pbuf, q), tile1(p, q));
              tile0(pbuf, q) = m;
              tile1(p, q) = m;
            }
          }
          break;

        case SyncOperation::SmoothBlend:
          for (int p = 0; p < halo; ++p)
          {
            float r = (halo > 1) ? float(p) / float(halo - 1) : 0.5f;
            r = smoothstep5(r);
            int pbuf = tile0.shape.x - 2 * halo + p;

            for (int q = 0; q < tile0.shape.y; ++q)
            {
              tile1(p, q) = lerp(tile0(pbuf, q), tile1(p, q), r);
              tile0(pbuf, q) = tile1(p, q);
            }
          }
          break;
        }
      }
    }

    // --- Pass 2: Y-direction stitching (operates on x-updated corners as well)
    for (int ty = 0; ty < ny - 1; ++ty)
    {
      for (int tx = 0; tx < nx; ++tx)
      {
        int idx0 = ty * nx + tx;
        int idx1 = (ty + 1) * nx + tx;

        Array *t0 = this->state->get_registered_tile(idx0);
        Array *t1 = this->state->get_registered_tile(idx1);
        if (!t0 || !t1) continue;

        Array &tile0 = *t0;
        Array &tile1 = *t1;

        switch (op)
        {
        case SyncOperation::Average:
        case SyncOperation::Mean:
          for (int q = 0; q < halo; ++q)
          {
            int qbuf = tile0.shape.y - 2 * halo + q;
            for (int p = 0; p < tile0.shape.x; ++p)
            {
              float avg = 0.5f * (tile0(p, qbuf) + tile1(p, q));
              tile0(p, qbuf) = avg;
              tile1(p, q) = avg;
            }
          }
          break;

        case SyncOperation::CopyFirst:
          for (int q = 0; q < halo; ++q)
          {
            int qbuf = tile0.shape.y - 2 * halo + q;
            for (int p = 0; p < tile0.shape.x; ++p)
            {
              tile1(p, q) = tile0(p, qbuf);
            }
          }
          break;

        case SyncOperation::CopySecond:
          for (int q = 0; q < halo; ++q)
          {
            int qbuf = tile0.shape.y - 2 * halo + q;
            for (int p = 0; p < tile0.shape.x; ++p)
            {
              tile0(p, qbuf) = tile1(p, q);
            }
          }
          break;

        case SyncOperation::Min:
          for (int q = 0; q < halo; ++q)
          {
            int qbuf = tile0.shape.y - 2 * halo + q;
            for (int p = 0; p < tile0.shape.x; ++p)
            {
              float m = std::min(tile0(p, qbuf), tile1(p, q));
              tile0(p, qbuf) = m;
              tile1(p, q) = m;
            }
          }
          break;

        case SyncOperation::Max:
          for (int q = 0; q < halo; ++q)
          {
            int qbuf = tile0.shape.y - 2 * halo + q;
            for (int p = 0; p < tile0.shape.x; ++p)
            {
              float m = std::max(tile0(p, qbuf), tile1(p, q));
              tile0(p, qbuf) = m;
              tile1(p, q) = m;
            }
          }
          break;

        case SyncOperation::SmoothBlend:
          for (int q = 0; q < halo; ++q)
          {
            float r = (halo > 1) ? float(q) / float(halo - 1) : 0.5f;
            r = smoothstep5(r);
            int qbuf = tile0.shape.y - 2 * halo + q;

            for (int p = 0; p < tile0.shape.x; ++p)
            {
              tile1(p, q) = lerp(tile0(p, qbuf), tile1(p, q), r);
              tile0(p, qbuf) = tile1(p, q);
            }
          }
          break;
        }
      }
    }
  }

  // Step 4: Final barrier so all ranks resume after Rank #0 finishes stitching
  this->barrier();
}

} // namespace hmap
