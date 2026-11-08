/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include "highmap/erosion/aeolian_erosion.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/operator.hpp"
#include "highmap/range.hpp"

namespace hmap
{

namespace
{

constexpr float pi = 3.14159265358979323846f;

} // namespace

void aeolian_erosion(Array        &z,
                     int           steps,
                     std::uint32_t seed,
                     glm::vec2     wind,
                     const Array  *p_bedrock,
                     const Array  *p_vegetation,
                     const Array  *p_abrasion_mask,
                     Array        *p_sediment_map,
                     Array        *p_bedrock_map,
                     float         matter_to_move,
                     int           samples_per_step,
                     float         sand_repose_angle,
                     float         shadow_angle_min,
                     float         shadow_angle_max,
                     float         max_shadow_dist,
                     float         abrasion_rate,
                     float         repose_rate)
{
  if (!validate_non_empty(z)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;
  if (p_vegetation && !validate_same_shape(z, *p_vegetation)) return;
  if (p_abrasion_mask && !validate_same_shape(z, *p_abrasion_mask)) return;

  const int nx = z.shape.x;
  const int ny = z.shape.y;
  const int total_cells = nx * ny;

  int n_samples = (samples_per_step > 0) ? samples_per_step : total_cells;

  // Bedrock layer: initialize from input if provided, otherwise zero
  Array bed(z.shape, 0.f);
  if (p_bedrock)
    bed = *p_bedrock;

  // Sediment layer: initial loose material is z - bedrock
  Array sed = z - bed;
  clamp_min(sed, 0.f);

  // Convert angles to tangents
  float tan_repose = std::tan(sand_repose_angle * pi / 180.f);
  float tan_shadow_min = std::tan(shadow_angle_min * pi / 180.f);
  float tan_shadow_max = std::tan(shadow_angle_max * pi / 180.f);

  // Buffer for double-buffered repose stabilization
  Array sed_temp(z.shape, 0.f);

  // Simulation loop
  for (int step = 0; step < steps; ++step)
  {
    // 1. Saltation, wind transport, and abrasion
    {
      auto run = clwrapper::Run("aeolian_transport");

      run.bind_buffer<float>("bedrock", bed.vector);
      run.bind_buffer<float>("sediment", sed.vector);
      gpu::helper_bind_optional_buffer(run, "vegetation", p_vegetation);
      gpu::helper_bind_optional_buffer(run, "abrasion_mask", p_abrasion_mask);

      run.bind_arguments(nx,
                         ny,
                         n_samples,
                         seed,
                         static_cast<std::uint32_t>(step),
                         wind.x,
                         wind.y,
                         matter_to_move,
                         tan_shadow_min,
                         tan_shadow_max,
                         max_shadow_dist,
                         abrasion_rate,
                         p_vegetation ? 1 : 0,
                         p_abrasion_mask ? 1 : 0);

      run.write_buffer("bedrock");
      run.write_buffer("sediment");
      run.execute(n_samples);
      run.read_buffer("bedrock");
      run.read_buffer("sediment");
    }

    // 2. Angle-of-repose stabilization pass
    {
      auto run_stab = clwrapper::Run("aeolian_stabilize_sediment");

      run_stab.bind_buffer<float>("bedrock", bed.vector);
      run_stab.bind_buffer<float>("sed_in", sed.vector);
      run_stab.bind_buffer<float>("sed_out", sed_temp.vector);

      run_stab.bind_arguments(nx, ny, tan_repose, repose_rate);

      run_stab.write_buffer("bedrock");
      run_stab.write_buffer("sed_in");
      run_stab.execute({nx, ny});
      run_stab.read_buffer("sed_out");

      sed = sed_temp;
    }
  }

  // Update total elevation
  z = bed + sed;

  if (p_sediment_map) *p_sediment_map = sed;
  if (p_bedrock_map) *p_bedrock_map = bed;
}

} // namespace hmap
