/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>

#include "highmap/erosion/hydraulic_erosion.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/math.hpp"
#include "highmap/morphology.hpp"

namespace hmap
{

Array flow_fixing_hydraulic_mise(const Array      &z,
                                 const MiseParams &params,
                                 float             flow_threshold,
                                 float             riverbed_talus,
                                 float             merging_distance,
                                 int               flow_fixing_iterations,
                                 int               prefilter_ir)
{
  // --- Stage 1: Pre-conditioning / flow fixing on working copy

  Array z_work = flow_fixing(z,
                             riverbed_talus,
                             flow_fixing_iterations,
                             prefilter_ir,
                             false,
                             merging_distance);

  Array mask_breach(z.shape);

  for (int j = 0; j < z.shape.y; ++j)
    for (int i = 0; i < z.shape.x; ++i)
      if (z_work(i, j) != z(i, j)) mask_breach(i, j) = 1.f;

  // --- Stage 2: Stream erodibility mask construction

  Array flow = flow_accumulation_d8(z_work);
  flow /= static_cast<float>(z.shape.x * z.shape.y);

  Array core_stream(z.shape);

  for (int j = 0; j < z.shape.y; ++j)
    for (int i = 0; i < z.shape.x; ++i)
      if (flow(i, j) >= flow_threshold || mask_breach(i, j) > 0.f)
        core_stream(i, j) = 1.f;

  Array dist = distance_transform(core_stream);
  Array erodibility(z.shape);
  float inv_dist = 1.f / std::max(merging_distance, 1e-6f);

  for (int j = 0; j < z.shape.y; ++j)
    for (int i = 0; i < z.shape.x; ++i)
    {
      float t = 1.f - std::clamp(dist(i, j) * inv_dist, 0.f, 1.f);
      erodibility(i, j) = t * t * (3.f - 2.f * t);
    }

  // --- Stage 3: Channel-focused hydraulic_mise pass using erodibility mask

  hydraulic_mise(z_work,
                 params,
                 nullptr,
                 &erodibility,
                 nullptr,
                 nullptr,
                 nullptr,
                 nullptr,
                 nullptr,
                 nullptr);

  // --- Stage 4: Final monotonicity enforcement

  return flow_fixing(z_work,
                     riverbed_talus,
                     1,
                     prefilter_ir,
                     false,
                     merging_distance);
}

} // namespace hmap
