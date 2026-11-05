/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file communicator.cpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2026
 */

#include "highmap/virtual_array/communicator.hpp"

#include <algorithm>

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

std::vector<float> SingleArrayComm::all_union(const std::vector<float> &local_vals)
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

SingleArrayComm &SingleArrayComm::instance()
{
  static SingleArrayComm inst;
  return inst;
}

// =====================================
// SequentialTileComm Implementation
// =====================================

SequentialTileComm::SequentialTileComm(int total_tiles)
    : total_tiles(total_tiles)
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
    case ReduceOp::Sum:
      this->accumulators[seq] += local_val;
      break;
    }
  }

  return this->accumulators[seq];
}

std::vector<float> SequentialTileComm::all_union(const std::vector<float> &local_vals)
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
    case ReduceOp::Sum:
      this->reduced_val += val;
      break;
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
    this->cv.wait(lock, [this, gen]() { return this->barrier_generation != gen; });
  }

  return this->reduced_val;
}

std::vector<float> ThreadedCollectiveState::all_union(const std::vector<float> &vals)
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
    this->cv.wait(lock, [this, gen]() { return this->barrier_generation != gen; });
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
    this->cv.wait(lock, [this, gen]() { return this->barrier_generation != gen; });
  }
}

int ThreadedCollectiveState::get_total_threads() const
{
  return this->total_threads;
}

// =====================================
// ThreadedTileComm Implementation
// =====================================

ThreadedTileComm::ThreadedTileComm(std::shared_ptr<ThreadedCollectiveState> state,
                                   int                                     thread_rank,
                                   int                                     total_size)
    : state(state),
      thread_rank(thread_rank),
      total_size(total_size)
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

std::vector<float> ThreadedTileComm::all_union(const std::vector<float> &local_vals)
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

} // namespace hmap
