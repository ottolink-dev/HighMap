/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>

#include "highmap/array.hpp"
#include "highmap/curvature.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/range.hpp"
#include "highmap/selector.hpp"

namespace hmap::gpu
{

Array select_soil_rocks(const Array   &z,
                        int            ir_max,
                        int            ir_min,
                        int            steps,
                        float          smaller_scales_weight,
                        ClampMode      curvature_clamp_mode,
                        float          curvature_clamping,
                        bool           normalize,
                        ICommunicator &comm)
{
  if (!validate_non_empty(z)) return Array();

  Array sr(z.shape); // output

  float actual_min = (ir_min == 0 ? 1.f : float(ir_min));
  float di = (std::log(float(ir_max)) - std::log(actual_min)) /
             float(steps - 1);
  float scale = 1.f;

  for (int k = 0; k < steps; ++k)
  {
    // log progression
    int ir = int(std::exp(std::log(actual_min) + k * di));

    if (ir_min == 0 && k == 0) ir = 0;

    // mean curvature on filtered field
    Array zf = z;
    gpu::smooth_cpulse(zf, ir);
    Array cm = -scale * gpu::curvature_quadric(zf, 0, CurvatureType::CT_MEAN) *
               ir;

    clamp(cm, curvature_clamping, curvature_clamp_mode);

    sr = maximum(sr, cm);

    scale /= smaller_scales_weight;
  }

  if (normalize)
  {
    normalize_and_clip(sr, 2.5f, true, true, comm);
  }

  return sr;
}

} // namespace hmap::gpu
