/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file flow_graph.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for FlowGraph structure.
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once
#include <cstdint>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/internal/radix_heap.hpp"

namespace hmap
{

enum FlowGraphFlag : uint8_t
{
  F_CLOSED = 1, // already reached by the flood
  F_PIT = 2,    // reached through the pit queue: flat or closed depression
  F_OUTLET = 4  // fixed base-level cell
};

struct FlowGraph
{
  int                  nx = 0, ny = 0, n = 0;
  Array                zf;    // routing surface: z with the depressions filled
  std::vector<int>     order; // flood order, receivers before donors
  std::vector<int>     rank;  // inverse permutation of order
  std::vector<int>     rcv;   // receiver of each cell, itself for an outlet
  std::vector<uint8_t> flag;
  Array                area;  // drainage area in cells, or accumulated rain
  Array                steep; // upstream sum of the relative slope excess

  std::vector<int> lake;        // lake index of a flooded cell, else -1
  std::vector<int> lake_outlet; // spill cell of each lake
  std::vector<int>
      lake_root; // first flooded cell visited upstream-to-downstream
  std::vector<float>  lake_volume, lake_sediment, lake_target;
  std::vector<float>  lake_lo, lake_hi, lake_level;
  std::vector<double> lake_acc;
  std::vector<int>    pitq; // all flooded cells, in flood order
  int                 pit_cells = 0;

  RadixHeap<int>   heap;         // exact flood
  std::vector<int> bhead, bnext; // bucket queue of the fast flood

  void  resize(int nx_, int ny_);
  float dist(int i) const;
  bool  submerged(int i, const Array &z) const;
  void  flood(const Array   &z_arr,
              const uint8_t *outlet,
              bool           fast,
              bool           seed_borders = true);
  void  receivers(float p, uint32_t seed);
  void  accumulate(const Array &z_arr, const Array *p_rain, float talus_cell);
};

} // namespace hmap
