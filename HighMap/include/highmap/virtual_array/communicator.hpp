/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file communicator.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2026
 */
#pragma once

#include <condition_variable>
#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include <unordered_map>

namespace hmap
{

// =====================================
// Collective Operations
// =====================================

enum class ReduceOp : uint8_t
{
  Min,
  Max,
  Sum
};

std::string to_string(ReduceOp op);

enum class SyncOperation
{
  Average,
  CopyFirst,
  CopySecond,
  Min,
  Max,
  Mean,
  SmoothBlend
};

// Forward declarations
class Array;
struct TileRegion;

// =====================================
// Communicator Interface
// =====================================

/**
 * @brief Abstract communicator interface providing MPI-style collective
 * directives.
 */
class ICommunicator
{
public:
  virtual ~ICommunicator() = default;

  virtual float              all_reduce(float local_val, ReduceOp op) = 0;
  virtual std::vector<float> all_union(
      const std::vector<float> &local_vals) = 0;
  virtual void barrier() = 0;
  virtual int  rank() const = 0;
  virtual int  size() const = 0;
  virtual void sync_overlap_buffers(
      Array            &local_tile,
      const TileRegion &region,
      SyncOperation     op = SyncOperation::SmoothBlend) = 0;
};

// =====================================
// Single Array Communicator (Standalone)
// =====================================

/**
 * @brief Default standalone communicator for single-array computation. All
 * collective calls are identity or no-ops, resulting in zero overhead.
 */
class SingleArrayComm final : public ICommunicator
{
public:
  float              all_reduce(float local_val, ReduceOp op) override;
  std::vector<float> all_union(const std::vector<float> &local_vals) override;
  void               barrier() override;
  int                rank() const override;
  int                size() const override;
  void               sync_overlap_buffers(
                    Array            &local_tile,
                    const TileRegion &region,
                    SyncOperation     op = SyncOperation::SmoothBlend) override;

  static SingleArrayComm &instance();
};

// =====================================
// Sequential Tile Communicator
// =====================================

/**
 * @brief Communicator for single-threaded sequential tile traversal. Tracks
 * reductions per collective sequence index so multiple distinct all_reduce
 * calls per tile (min, max, sum...) accumulate correctly.
 */
class SequentialTileComm final : public ICommunicator
{
public:
  explicit SequentialTileComm(int total_tiles, int halo = 0);

  float              all_reduce(float local_val, ReduceOp op) override;
  std::vector<float> all_union(const std::vector<float> &local_vals) override;
  void               barrier() override;
  int                rank() const override;
  int                size() const override;
  void               sync_overlap_buffers(
                    Array            &local_tile,
                    const TileRegion &region,
                    SyncOperation     op = SyncOperation::SmoothBlend) override;

  float get_accumulated_value(size_t seq = 0) const;
  void  next_tile(int tile_idx);

private:
  int                total_tiles = 1;
  int                current_tile_idx = 0;
  int                halo = 0;
  size_t             call_sequence = 0;
  std::vector<float> accumulators;
  std::vector<bool>  initialized;
  std::vector<float> accumulated_unique;
  // Boundary edge buffers for sequential neighbor exchange (keyed by tile
  // index) east_edges[ty][tx] holds the eastern overlapping halo column
  // south_edges[ty][tx] holds the southern overlapping halo row
  std::vector<std::vector<float>> east_edges;
  std::vector<std::vector<float>> south_edges;
};

// =====================================
// Threaded Collective State & Communicator
// =====================================

/**
 * @brief Shared synchronization state for multi-threaded tile communication.
 */
class ThreadedCollectiveState
{
public:
  explicit ThreadedCollectiveState(int num_threads);

  float              all_reduce(size_t seq, float val, ReduceOp op);
  std::vector<float> all_union(size_t seq, const std::vector<float> &vals);
  void               barrier(size_t seq);
  int                get_total_threads() const;

  void   register_tile(int tile_idx, Array *tile, const TileRegion &region);
  Array *get_registered_tile(int tile_idx) const;
  const TileRegion *get_registered_region(int tile_idx) const;

private:
  struct ReduceSlot
  {
    float value = 0.f;
    int   remaining = 0;
    bool  ready = false;
  };

  struct UnionSlot
  {
    std::vector<float> values;
    int                remaining = 0;
    bool               ready = false;
  };

  struct BarrierSlot
  {
    int  remaining = 0;
    bool ready = false;
  };

  mutable std::mutex                      mtx;
  std::condition_variable                 cv;
  int                                     total_threads;
  std::unordered_map<size_t, ReduceSlot>  reduce_slots;
  std::unordered_map<size_t, UnionSlot>   union_slots;
  std::unordered_map<size_t, BarrierSlot> barrier_slots;
  std::unordered_map<int, Array *>        registered_tiles;
  std::unordered_map<int, TileRegion>     registered_regions;
};

/**
 * @brief Thread-local communicator proxy connected to a
 * ThreadedCollectiveState.
 */
class ThreadedTileComm final : public ICommunicator
{
public:
  ThreadedTileComm(std::shared_ptr<ThreadedCollectiveState> state,
                   int                                      thread_rank,
                   int                                      total_size,
                   int                                      grid_nx = 1);

  float              all_reduce(float local_val, ReduceOp op) override;
  std::vector<float> all_union(const std::vector<float> &local_vals) override;
  void               barrier() override;
  int                rank() const override;
  int                size() const override;
  void               sync_overlap_buffers(
                    Array            &local_tile,
                    const TileRegion &region,
                    SyncOperation     op = SyncOperation::SmoothBlend) override;

private:
  std::shared_ptr<ThreadedCollectiveState> state;
  int                                      thread_rank;
  int                                      total_size;
  int                                      grid_nx;
  size_t                                   reduce_sequence = 0;
  size_t                                   union_sequence = 0;
  size_t                                   barrier_sequence = 0;
};

} // namespace hmap
