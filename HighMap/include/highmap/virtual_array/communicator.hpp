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
  virtual std::vector<float> all_union(const std::vector<float> &local_vals) = 0;
  virtual void               barrier() = 0;
  virtual int                rank() const = 0;
  virtual int                size() const = 0;
};

// =====================================
// Single Array Communicator (Standalone)
// =====================================

/**
 * @brief Default standalone communicator for single-array computation.
 * All collective calls are identity or no-ops, resulting in zero overhead.
 */
class SingleArrayComm final : public ICommunicator
{
public:
  float              all_reduce(float local_val, ReduceOp op) override;
  std::vector<float> all_union(const std::vector<float> &local_vals) override;
  void               barrier() override;
  int                rank() const override;
  int                size() const override;

  static SingleArrayComm &instance();
};

// =====================================
// Sequential Tile Communicator
// =====================================

/**
 * @brief Communicator for single-threaded sequential tile traversal.
 * Tracks reductions per collective sequence index so multiple distinct
 * all_reduce calls per tile (min, max, sum...) accumulate correctly.
 */
class SequentialTileComm final : public ICommunicator
{
public:
  explicit SequentialTileComm(int total_tiles);

  float              all_reduce(float local_val, ReduceOp op) override;
  std::vector<float> all_union(const std::vector<float> &local_vals) override;
  void               barrier() override;
  int                rank() const override;
  int                size() const override;

  float get_accumulated_value(size_t seq = 0) const;
  void  next_tile(int tile_idx);

private:
  int                total_tiles = 1;
  int                current_tile_idx = 0;
  size_t             call_sequence = 0;
  std::vector<float> accumulators;
  std::vector<bool>  initialized;
  std::vector<float> accumulated_unique;
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

  float              all_reduce(float val, ReduceOp op);
  std::vector<float> all_union(const std::vector<float> &vals);
  void               barrier();
  int                get_total_threads() const;

private:
  std::mutex              mtx;
  std::condition_variable cv;
  int                     total_threads;
  int                     remaining_threads;
  int                     barrier_generation;
  float                   reduced_val;
  std::vector<float>      union_vals;
};

/**
 * @brief Thread-local communicator proxy connected to a ThreadedCollectiveState.
 */
class ThreadedTileComm final : public ICommunicator
{
public:
  ThreadedTileComm(std::shared_ptr<ThreadedCollectiveState> state,
                   int                                     thread_rank,
                   int                                     total_size);

  float              all_reduce(float local_val, ReduceOp op) override;
  std::vector<float> all_union(const std::vector<float> &local_vals) override;
  void               barrier() override;
  int                rank() const override;
  int                size() const override;

private:
  std::shared_ptr<ThreadedCollectiveState> state;
  int                                     thread_rank;
  int                                     total_size;
};

} // namespace hmap
