/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/morphology.hpp"
#include "highmap/range.hpp"
#include "highmap/selector.hpp"

namespace hmap::gpu
{

Array select_soil_flow(const Array   &z,
                       int            ir_gradient,
                       float          gradient_weight,
                       float          gradient_scaling_factor,
                       float          flow_weight,
                       float          talus_ref,
                       float          clipping_ratio,
                       float          flow_gamma,
                       float          k_smooth,
                       bool           normalize,
                       ICommunicator &comm)
{
  if (!validate_non_empty(z)) return Array();

  if (gradient_scaling_factor <= 0.f) gradient_scaling_factor = z.shape.x;
  if (talus_ref <= 0.f) talus_ref = 1.f / z.shape.x;

  // gradient (scaling is empirical)
  Array dn = gpu::morphological_gradient(z, ir_gradient) *
             gradient_scaling_factor / 32.f / ir_gradient;
  clamp_max_smooth(dn, 1.f, 0.01f);

  Array sr = select_rivers(z, talus_ref, clipping_ratio);

  if (flow_gamma != 1.f) gamma_correction(sr, flow_gamma);

  Array c = hmap::maximum_smooth(gradient_weight * dn, sr, k_smooth);
  c *= (flow_weight + dn) / (flow_weight + 1.f);

  if (normalize)
  {
    normalize_and_clip(c, 2.5f, true, true, comm);
  }

  return c;
}

} // namespace hmap::gpu
