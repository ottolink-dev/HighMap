/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file communicator.cpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2026
 */

#include <algorithm>

#include "highmap/array.hpp"
#include "highmap/math/core.hpp"
#include "highmap/virtual_array/communicator.hpp"
#include "highmap/virtual_array/tile_region.hpp"

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

SequentialTileComm::SequentialTileComm(int total_tiles, int halo)
    : total_tiles(total_tiles), halo(halo)
{
}

float SequentialTileComm::all_reduce(float local_val, ReduceOp op)
{
  size_t seq = this->call_sequence++;

  if (seq >= this->accumulators.size())
  {
    this->accumulators.resize(seq + 1);
    this->initialized.resize(seq + 1, false);
  }

  if (!this->initialized[seq])
  {
    this->accumulators[seq] = local_val;
    this->initialized[seq] = true;
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

  return this->accumulators[seq];
}

std::vector<float> SequentialTileComm::all_union(
    const std::vector<float> &local_vals)
{
  for (float v : local_vals)
  {
    this->accumulated_unique.push_back(v);
  }

  std::sort(this->accumulated_unique.begin(), this->accumulated_unique.end());
  auto last = std::unique(this->accumulated_unique.begin(),
                          this->accumulated_unique.end());
  this->accumulated_unique.erase(last, this->accumulated_unique.end());

  return this->accumulated_unique;
}

void SequentialTileComm::barrier()
{
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

  // In sequential traversal (raster order tx then ty), a tile can receive
  // incoming overlap data from its west neighbor (tx - 1) and north neighbor
  // (ty - 1), and send its boundary overlap data to its east neighbor (tx + 1)
  // and south neighbor (ty + 1).

  const int tx = region.key.tx;
  const int ty = region.key.ty;

  // 1. Ingest Western overlap from (tx - 1, ty) if present
  if (tx > 0 && !this->east_edges.empty())
  {
    size_t prev_tile_key = static_cast<size_t>(ty * 10000 + (tx - 1));
    if (prev_tile_key < this->east_edges.size() &&
        !this->east_edges[prev_tile_key].empty())
    {
      const auto &src_col = this->east_edges[prev_tile_key];
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
  if (ty > 0 && !this->south_edges.empty())
  {
    size_t prev_tile_key = static_cast<size_t>((ty - 1) * 10000 + tx);
    if (prev_tile_key < this->south_edges.size() &&
        !this->south_edges[prev_tile_key].empty())
    {
      const auto &src_row = this->south_edges[prev_tile_key];
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
  if (curr_key >= this->east_edges.size())
  {
    this->east_edges.resize(curr_key + 1);
  }
  auto &east_buf = this->east_edges[curr_key];
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
  if (curr_key >= this->south_edges.size())
  {
    this->south_edges.resize(curr_key + 1);
  }
  auto &south_buf = this->south_edges[curr_key];
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

// =====================================
// ThreadedCollectiveState Implementation
// =====================================

ThreadedCollectiveState::ThreadedCollectiveState(int num_threads)
    : total_threads(num_threads),
      remaining_threads(num_threads),
      barrier_generation(0),
      reduced_val(0.f)
{
}

float ThreadedCollectiveState::all_reduce(float val, ReduceOp op)
{
  std::unique_lock<std::mutex> lock(this->mtx);
  int                          gen = this->barrier_generation;

  if (this->remaining_threads == this->total_threads)
  {
    this->reduced_val = val;
  }
  else
  {
    switch (op)
    {
    case ReduceOp::Min:
      this->reduced_val = std::min(this->reduced_val, val);
      break;
    case ReduceOp::Max:
      this->reduced_val = std::max(this->reduced_val, val);
      break;
    case ReduceOp::Sum: this->reduced_val += val; break;
    }
  }

  if (--this->remaining_threads == 0)
  {
    this->remaining_threads = this->total_threads;
    this->barrier_generation++;
    this->cv.notify_all();
  }
  else
  {
    this->cv.wait(lock,
                  [this, gen]() { return this->barrier_generation != gen; });
  }

  return this->reduced_val;
}

std::vector<float> ThreadedCollectiveState::all_union(
    const std::vector<float> &vals)
{
  std::unique_lock<std::mutex> lock(this->mtx);
  int                          gen = this->barrier_generation;

  for (float v : vals)
  {
    this->union_vals.push_back(v);
  }

  if (--this->remaining_threads == 0)
  {
    std::sort(this->union_vals.begin(), this->union_vals.end());
    auto last = std::unique(this->union_vals.begin(), this->union_vals.end());
    this->union_vals.erase(last, this->union_vals.end());

    this->remaining_threads = this->total_threads;
    this->barrier_generation++;
    this->cv.notify_all();
  }
  else
  {
    this->cv.wait(lock,
                  [this, gen]() { return this->barrier_generation != gen; });
  }

  return this->union_vals;
}

void ThreadedCollectiveState::barrier()
{
  std::unique_lock<std::mutex> lock(this->mtx);
  int                          gen = this->barrier_generation;

  if (--this->remaining_threads == 0)
  {
    this->remaining_threads = this->total_threads;
    this->barrier_generation++;
    this->cv.notify_all();
  }
  else
  {
    this->cv.wait(lock,
                  [this, gen]() { return this->barrier_generation != gen; });
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
    return this->state->all_reduce(local_val, op);
  }
  return local_val;
}

std::vector<float> ThreadedTileComm::all_union(
    const std::vector<float> &local_vals)
{
  if (this->state)
  {
    return this->state->all_union(local_vals);
  }
  return local_vals;
}

void ThreadedTileComm::barrier()
{
  if (this->state)
  {
    this->state->barrier();
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
