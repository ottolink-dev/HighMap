/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/array.hpp"
#include "highmap/curvature.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/morphology.hpp"
#include "highmap/range.hpp"
#include "highmap/selector.hpp"

namespace hmap::gpu
{

Array select_soil_weathered(const Array   &z,
                            int            ir_curvature,
                            int            ir_gradient,
                            ClampMode      curvature_clamp_mode,
                            float          curvature_clamping,
                            float          curvature_weight,
                            float          gradient_weight,
                            float          gradient_scaling_factor,
                            bool           normalize,
                            ICommunicator &comm)
{
  if (!validate_non_empty(z)) return Array();

  if (gradient_scaling_factor <= 0.f) gradient_scaling_factor = z.shape.x;

  // gradient (scaling is empirical)
  Array dn = gpu::morphological_gradient(z, ir_gradient) *
             gradient_scaling_factor / 32.f / ir_gradient;

  return select_soil_weathered(z,
                               dn,
                               ir_curvature,
                               curvature_clamp_mode,
                               curvature_clamping,
                               curvature_weight,
                               gradient_weight,
                               gradient_scaling_factor,
                               normalize,
                               comm);
}

Array select_soil_weathered(const Array   &z,
                            const Array   &gradient_norm,
                            int            ir_curvature,
                            ClampMode      curvature_clamp_mode,
                            float          curvature_clamping,
                            float          curvature_weight,
                            float          gradient_weight,
                            float          gradient_scaling_factor,
                            bool           normalize,
                            ICommunicator &comm)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, gradient_norm))
    return Array();

  // curvature
  Array cm = z;

  if (ir_curvature) gpu::smooth_cpulse(cm, ir_curvature);

  cm = -gradient_scaling_factor *
       gpu::curvature_quadric(cm, 0, CurvatureType::CT_MEAN);
  clamp(cm, curvature_clamping, curvature_clamp_mode);

  Array out = curvature_weight * cm + gradient_weight * gradient_norm;

  if (normalize)
  {
    normalize_and_clip(out, 2.5f, true, true, comm);
  }

  return out;
}

} // namespace hmap::gpu
