/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <future>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap::va
{

void smooth_cpulse(VirtualArray       &array,
                   int                 ir,
                   const VirtualArray *p_mask,
                   const ComputeMode  &cm)
{
  hmap::for_each_tile(
      {p_mask},
      {&array},
      [&ir](std::vector<const hmap::Array *> p_arrays_in,
            std::vector<hmap::Array *>       p_arrays_out,
            const hmap::TileRegion &)
      {
        const auto [pa_mask] = unpack<1>(p_arrays_in);
        auto [pa_out] = unpack<1>(p_arrays_out);

        hmap::gpu::smooth_cpulse(*pa_out, ir, pa_mask);
      },
      cm);
}

VirtualArray smooth_cpulse(const VirtualArray &array,
                           int                 ir,
                           const VirtualArray *p_mask,
                           const ComputeMode  &cm)
{
  VirtualArray out;
  out.copy_from(array, cm, /* copy_src_data */ ir == 0);

  if (ir == 0) return out;

  hmap::for_each_tile(
      {&array, p_mask},
      {&out},
      [&ir](std::vector<const hmap::Array *> p_arrays_in,
            std::vector<hmap::Array *>       p_arrays_out,
            const hmap::TileRegion &)
      {
        const auto [pa_in, pa_mask] = unpack<2>(p_arrays_in);
        auto [pa_out] = unpack<1>(p_arrays_out);
        *pa_out = *pa_in;

        hmap::gpu::smooth_cpulse(*pa_out, ir, pa_mask);
      },
      cm);

  return out;
}

} // namespace hmap::va
