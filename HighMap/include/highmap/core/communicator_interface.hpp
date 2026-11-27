/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file communicator_interface.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Lightweight abstract SPMD communicator interface for Array and
 * VirtualArray collective operations.
 * @copyright Copyright (c) 2026
 */
#pragma once

#include <cstdint>
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

enum class SyncOperation
{
  Average,
  CopyFirst,
  CopySecond,
  Min,
  Max,
  Mean,
  SmoothBlend,
  Sum
};

// Forward declarations
class Array;
struct UvTile;
using UvTile = UvTile;

// =====================================
// Communicator Interface
// =====================================

/**
 * @brief Abstract communicator interface providing MPI-style collective
 * directives for single-source SPMD Array algorithms.
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
      Array        &local_tile,
      const UvTile &tile,
      SyncOperation op = SyncOperation::SmoothBlend) = 0;
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
                    Array        &local_tile,
                    const UvTile &tile,
                    SyncOperation op = SyncOperation::SmoothBlend) override;

  static SingleArrayComm &instance();
};

} // namespace hmap
