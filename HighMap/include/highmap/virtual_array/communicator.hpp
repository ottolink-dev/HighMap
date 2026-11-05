/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file communicator.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2026
 */
#pragma once

#include <algorithm>
#include <condition_variable>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <string>

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

inline std::string to_string(ReduceOp op)
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

  virtual float all_reduce(float local_val, ReduceOp op) = 0;
  virtual void  barrier() = 0;
  virtual int   rank() const = 0;
  virtual int   size() const = 0;
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
  float all_reduce(float local_val, ReduceOp /*op*/) override
  {
    return local_val;
  }

  void barrier() override
  {
  }

  int rank() const override
  {
    return 0;
  }

  int size() const override
  {
    return 1;
  }

  static SingleArrayComm &instance()
  {
    static SingleArrayComm inst;
    return inst;
  }
};

// =====================================
// Sequential Tile Communicator
// =====================================

/**
 * @brief Communicator for single-threaded sequential tile traversal.
 * Tiles stream one by one so only 1 tile resides in memory at a time.
 */
class SequentialTileComm final : public ICommunicator
{
public:
  SequentialTileComm(int total_tiles)
      : total_tiles(total_tiles)
  {
  }

  float all_reduce(float local_val, ReduceOp op) override
  {
    if (!reduced)
    {
      switch (op)
      {
      case ReduceOp::Min:
        accumulated_val = std::min(accumulated_val, local_val);
        break;
      case ReduceOp::Max:
        accumulated_val = std::max(accumulated_val, local_val);
        break;
      case ReduceOp::Sum:
        accumulated_val += local_val;
        break;
      }
      reduced = true;
    }

    return accumulated_val;
  }

  void barrier() override
  {
  }

  int rank() const override
  {
    return current_tile_idx;
  }

  int size() const override
  {
    return total_tiles;
  }

  void next_tile(int tile_idx)
  {
    current_tile_idx = tile_idx;
    reduced = false;
  }

  float get_accumulated_value() const
  {
    return accumulated_val;
  }

private:
  int   total_tiles = 1;
  int   current_tile_idx = 0;
  bool  reduced = false;
  float accumulated_val = std::numeric_limits<float>::infinity();
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
  explicit ThreadedCollectiveState(int num_threads)
      : total_threads(num_threads),
        remaining_threads(num_threads),
        barrier_generation(0),
        reduced_val(0.f)
  {
  }

  void barrier()
  {
    std::unique_lock<std::mutex> lock(mtx);
    int                          gen = barrier_generation;

    if (--remaining_threads == 0)
    {
      remaining_threads = total_threads;
      barrier_generation++;
      cv.notify_all();
    }
    else
    {
      cv.wait(lock, [this, gen]() { return barrier_generation != gen; });
    }
  }

  float all_reduce(float val, ReduceOp op)
  {
    std::unique_lock<std::mutex> lock(mtx);
    int                          gen = barrier_generation;

    if (remaining_threads == total_threads)
    {
      // First thread in this round initializes reduced_val
      reduced_val = val;
    }
    else
    {
      switch (op)
      {
      case ReduceOp::Min:
        reduced_val = std::min(reduced_val, val);
        break;
      case ReduceOp::Max:
        reduced_val = std::max(reduced_val, val);
        break;
      case ReduceOp::Sum:
        reduced_val += val;
        break;
      }
    }

    if (--remaining_threads == 0)
    {
      remaining_threads = total_threads;
      barrier_generation++;
      cv.notify_all();
    }
    else
    {
      cv.wait(lock, [this, gen]() { return barrier_generation != gen; });
    }

    return reduced_val;
  }

  int get_total_threads() const
  {
    return total_threads;
  }

private:
  std::mutex              mtx;
  std::condition_variable cv;
  int                     total_threads;
  int                     remaining_threads;
  int                     barrier_generation;
  float                   reduced_val;
};

/**
 * @brief Thread-local communicator proxy connected to a ThreadedCollectiveState.
 */
class ThreadedTileComm final : public ICommunicator
{
public:
  ThreadedTileComm(std::shared_ptr<ThreadedCollectiveState> state,
                   int                                     thread_rank,
                   int                                     total_size)
      : state(state),
        thread_rank(thread_rank),
        total_size(total_size)
  {
  }

  float all_reduce(float local_val, ReduceOp op) override
  {
    if (state)
    {
      return state->all_reduce(local_val, op);
    }
    return local_val;
  }

  void barrier() override
  {
    if (state)
    {
      state->barrier();
    }
  }

  int rank() const override
  {
    return thread_rank;
  }

  int size() const override
  {
    return total_size;
  }

private:
  std::shared_ptr<ThreadedCollectiveState> state;
  int                                     thread_rank;
  int                                     total_size;
};

} // namespace hmap
