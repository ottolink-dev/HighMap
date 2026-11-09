/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/erosion/aeolian_erosion.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/operator.hpp"
#include "highmap/range.hpp"

namespace hmap
{

void aeolian_erosion(Array        &z,
                     int           steps,
                     std::uint32_t seed,
                     glm::vec2     wind,
                     const Array  *p_bedrock,
                     const Array  *p_vegetation,
                     const Array  *p_abrasion_mask,
                     Array        *p_sediment_map,
                     Array        *p_bedrock_map,
                     float         cell_size,
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

  if (cell_size <= 0.f) cell_size = 1.f / static_cast<float>(nx);

  int n_samples = (samples_per_step > 0) ? samples_per_step : total_cells;

  // Bedrock layer: initialize from input if provided, otherwise zero
  Array bed(z.shape, 0.f);
  if (p_bedrock) bed = *p_bedrock;

  // Sediment layer: initial loose material is z - bedrock
  Array sed = z - bed;
  clamp_min(sed, 0.f);

  // Convert angles to tangents
  float tan_repose = std::tan(sand_repose_angle * static_cast<float>(M_PI) /
                              180.f);
  float tan_shadow_min = std::tan(shadow_angle_min * static_cast<float>(M_PI) /
                                  180.f);
  float tan_shadow_max = std::tan(shadow_angle_max * static_cast<float>(M_PI) /
                                  180.f);

  // Buffer for double-buffered repose stabilization
  Array sed_temp(z.shape, 0.f);

  // Initialize kernels with double-buffered sediment resident on GPU
  auto run_transport_a = clwrapper::Run("aeolian_transport");

  run_transport_a.bind_buffer<float>("bedrock", bed.vector);
  run_transport_a.bind_buffer<float>("sediment", sed.vector);
  gpu::helper_bind_optional_buffer(run_transport_a, "vegetation", p_vegetation);
  gpu::helper_bind_optional_buffer(run_transport_a,
                                   "abrasion_mask",
                                   p_abrasion_mask);

  run_transport_a.bind_arguments(nx,
                                 ny,
                                 n_samples,
                                 seed,
                                 0u,
                                 wind.x,
                                 wind.y,
                                 matter_to_move,
                                 cell_size,
                                 tan_shadow_min,
                                 tan_shadow_max,
                                 max_shadow_dist,
                                 abrasion_rate,
                                 p_vegetation ? 1 : 0,
                                 p_abrasion_mask ? 1 : 0);

  run_transport_a.write_buffer("bedrock");
  run_transport_a.write_buffer("sediment");

  // run_stab_ab: repose stabilization reading sed_a and writing sed_b
  auto run_stab_ab = clwrapper::Run("aeolian_stabilize_sediment",
                                    run_transport_a.get_queue());

  run_stab_ab.bind_buffer("bedrock", run_transport_a.get_buffer("bedrock"));
  run_stab_ab.bind_buffer("sed_in", run_transport_a.get_buffer("sediment"));
  run_stab_ab.bind_buffer<float>("sed_out", sed_temp.vector);

  run_stab_ab.bind_arguments(nx, ny, cell_size, tan_repose, repose_rate);

  // run_transport_b: saltation pass operating on sed_b (sed_out of stab_ab)
  auto run_transport_b = clwrapper::Run("aeolian_transport",
                                        run_transport_a.get_queue());

  run_transport_b.bind_buffer("bedrock", run_transport_a.get_buffer("bedrock"));
  run_transport_b.bind_buffer("sediment", run_stab_ab.get_buffer("sed_out"));
  run_transport_b.bind_buffer("vegetation",
                              run_transport_a.get_buffer("vegetation"));
  run_transport_b.bind_buffer("abrasion_mask",
                              run_transport_a.get_buffer("abrasion_mask"));

  run_transport_b.bind_arguments(nx,
                                 ny,
                                 n_samples,
                                 seed,
                                 0u,
                                 wind.x,
                                 wind.y,
                                 matter_to_move,
                                 cell_size,
                                 tan_shadow_min,
                                 tan_shadow_max,
                                 max_shadow_dist,
                                 abrasion_rate,
                                 p_vegetation ? 1 : 0,
                                 p_abrasion_mask ? 1 : 0);

  // run_stab_ba: repose stabilization reading sed_b and writing sed_a
  auto run_stab_ba = clwrapper::Run("aeolian_stabilize_sediment",
                                    run_transport_a.get_queue());

  run_stab_ba.bind_buffer("bedrock", run_transport_a.get_buffer("bedrock"));
  run_stab_ba.bind_buffer("sed_in", run_stab_ab.get_buffer("sed_out"));
  run_stab_ba.bind_buffer("sed_out", run_transport_a.get_buffer("sediment"));

  run_stab_ba.bind_arguments(nx, ny, cell_size, tan_repose, repose_rate);

  // Simulation loop: entire state remains on GPU across all steps
  for (int step = 0; step < steps; ++step)
  {
    if (step % 2 == 0)
    {
      run_transport_a.set_argument(8, static_cast<std::uint32_t>(step));
      run_transport_a.execute(n_samples);

      run_stab_ab.execute({nx, ny});
    }
    else
    {
      run_transport_b.set_argument(8, static_cast<std::uint32_t>(step));
      run_transport_b.execute(n_samples);

      run_stab_ba.execute({nx, ny});
    }
  }

  // Read back results
  run_transport_a.read_buffer("bedrock");

  if (steps % 2 == 1)
  {
    run_stab_ab.read_buffer("sed_out");
    sed = sed_temp;
  }
  else if (steps > 0)
  {
    run_transport_a.read_buffer("sediment");
  }

  // Update total elevation
  z = bed + sed;

  if (p_sediment_map) *p_sediment_map = sed;
  if (p_bedrock_map) *p_bedrock_map = bed;
}

} // namespace hmap
