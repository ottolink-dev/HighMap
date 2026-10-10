/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/erosion/strata_erosion.hpp"
#include "highmap/filters.hpp"
#include "highmap/hydrology/drainage_basin_cell_based.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/core.hpp"
#include "highmap/morphology.hpp"
#include "highmap/transform.hpp"

namespace hmap::gpu
{

Array watershed_ridge(const Array        &z,
                      float               amplitude,
                      float               width,
                      float               edt_exponent,
                      int                 prefilter_ir,
                      FlowDirectionMethod fd_method,
                      const Array        *p_noise_x,
                      const Array        *p_noise_y,
                      const Array        *p_scaling)
{
  if (!validate_non_empty(z)) return Array();
  if (p_noise_x && !validate_same_shape(z, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(z, *p_noise_y)) return Array();
  if (p_scaling && !validate_same_shape(z, *p_scaling)) return Array();

  const glm::ivec2 shape = z.shape;

  // --- distance transform to main channel

  // retrieve main channel (valley bottom)
  Array zf = z;
  if (prefilter_ir > 0) gpu::smooth_cpulse(zf, prefilter_ir);

  // --- stream structure

  auto db = DrainageBasinCellBased(z);
  db.set_outlets(find_flow_sinks_border(z));

  if (fd_method == FlowDirectionMethod::FDM_D8)
  {
    db.compute_receivers();
    auto [subroots, has_lake] = db.find_subroots();
    if (has_lake) db.remove_lakes(subroots);
  }
  else
    db.compute_receivers_priority_flood();

  db.update_traversals();
  auto mcs = db.get_main_channels();

  // distance transform to generate valley shape
  Array edt(shape);

  for (const auto &vec : mcs)
    for (const auto &p : vec)
      edt(p) = 1.f;

  edt = distance_transform(edt);

  gpu::warp(edt, p_noise_x, p_noise_y);

  // --- generate rift profile

  for_each_cell(edt,
                [&](int i, int j, float &d)
                {
                  float sc = p_scaling ? (*p_scaling)(i, j) : 1.f;
                  float r = std::clamp(d / (sc * width), 0.f, 1.f);
                  r = smoothstep3_upper(r);
                  d = sc * (1.f - std::pow(r, edt_exponent));
                });

  // --- apply and output

  Array out = z - amplitude * edt;
  return out;
}

Array watershed_ridge(const Array        &z,
                      const Array        *p_mask,
                      float               amplitude,
                      float               width,
                      float               edt_exponent,
                      int                 prefilter_ir,
                      FlowDirectionMethod fd_method,
                      const Array        *p_noise_x,
                      const Array        *p_noise_y,
                      const Array        *p_scaling)
{
  if (!validate_non_empty(z)) return Array();
  if (p_mask && !validate_same_shape(z, *p_mask)) return Array();
  if (p_noise_x && !validate_same_shape(z, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(z, *p_noise_y)) return Array();
  if (p_scaling && !validate_same_shape(z, *p_scaling)) return Array();

  if (!p_mask)
    return watershed_ridge(z,
                           amplitude,
                           width,
                           edt_exponent,
                           prefilter_ir,
                           fd_method,
                           p_noise_x,
                           p_noise_y,
                           p_scaling);
  else
  {
    Array z_f = watershed_ridge(z,
                                amplitude,
                                width,
                                edt_exponent,
                                prefilter_ir,
                                fd_method,
                                p_noise_x,
                                p_noise_y,
                                p_scaling);
    return lerp(z, z_f, *(p_mask));
  }
}

} // namespace hmap::gpu
