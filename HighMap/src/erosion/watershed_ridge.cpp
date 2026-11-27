/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/erosion/strata_erosion.hpp"
#include "highmap/filters.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/core.hpp"
#include "highmap/morphology.hpp"
#include "highmap/transform.hpp"

namespace hmap::gpu
{

Array watershed_ridge(const Array &z,
                      float        amplitude,
                      float        width,
                      float        edt_exponent,
                      int          prefilter_ir,
                      bool         fast_approx,
                      const Array *p_noise_x,
                      const Array *p_noise_y,
                      const Array *p_scaling)
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

  FlowGraph fg;
  fg.resize(shape.x, shape.y);

  std::vector<uint8_t> outlets(fg.n, 0);
  for (const auto &p : find_flow_sinks_border(z))
    outlets[p.y * shape.x + p.x] = 1;

  fg.flood(z, outlets.data(), fast_approx, false);
  fg.receivers(1.f, 0);
  fg.accumulate(z, nullptr, 1.f);

  std::vector<int> best_child(fg.n, -1);
  for (int i = 0; i < fg.n; ++i)
  {
    const int r = fg.rcv[i];
    if (r == i) continue;

    if (best_child[r] < 0 || fg.area.vector[i] > fg.area.vector[best_child[r]])
      best_child[r] = i;
  }

  // distance transform to generate valley shape
  Array edt(shape);

  for (int i = 0; i < fg.n; ++i)
  {
    if (fg.rcv[i] != i) continue;

    int cur = i;
    while (cur >= 0)
    {
      edt.vector[cur] = 1.f;
      cur = best_child[cur];
    }
  }

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

Array watershed_ridge(const Array &z,
                      const Array *p_mask,
                      float        amplitude,
                      float        width,
                      float        edt_exponent,
                      int          prefilter_ir,
                      bool         fast_approx,
                      const Array *p_noise_x,
                      const Array *p_noise_y,
                      const Array *p_scaling)
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
                           fast_approx,
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
                                fast_approx,
                                p_noise_x,
                                p_noise_y,
                                p_scaling);
    return lerp(z, z_f, *(p_mask));
  }
}

} // namespace hmap::gpu
