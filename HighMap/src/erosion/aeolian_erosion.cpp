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
#include "highmap/virtual_array/communicator.hpp"
#include "highmap/virtual_array/tile_region.hpp"

namespace hmap
{

void aeolian_erosion(Array            &z,
                     int               steps,
                     std::uint32_t     seed,
                     glm::vec2         wind,
                     const Array      *p_bedrock,
                     const Array      *p_vegetation,
                     const Array      *p_abrasion_mask,
                     Array            *p_sediment_map,
                     Array            *p_bedrock_map,
                     float             cell_size,
                     float             matter_to_move,
                     int               samples_per_step,
                     float             sand_repose_angle,
                     float             shadow_angle_min,
                     float             shadow_angle_max,
                     float             max_shadow_dist,
                     float             abrasion_rate,
                     float             repose_rate,
                     int               sync_interval,
                     const TileRegion &region,
                     ICommunicator    &comm)
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

  // Resolution-invariant spatial scaling heuristic calibrated at reference grid
  // 256
  constexpr float ref_nx = 256.f;
  const float     scale = static_cast<float>(nx) / ref_nx;
  glm::vec2       eff_wind_base = wind * scale;
  float           eff_max_shadow = max_shadow_dist * scale;

  // Buffer for double-buffered repose stabilization
  Array sed_temp(z.shape, 0.f);

  // Simulation loop: processed in batches between boundary synchronizations
  // so GPU buffers and host memory are released while other tiles execute
  const int batch_size = (comm.size() > 1 && sync_interval > 0) ? sync_interval
                                                                : steps;

  for (int step_start = 0; step_start < steps; step_start += batch_size)
  {
    const int step_end = std::min(steps, step_start + batch_size);
    const int steps_in_batch = step_end - step_start;

    {
      // Initialize kernels with double-buffered sediment resident on GPU
      auto run_transport_a = clwrapper::Run("aeolian_transport");

      run_transport_a.bind_buffer<float>("bedrock", bed.vector);
      run_transport_a.bind_buffer<float>("sediment", sed.vector);
      gpu::helper_bind_optional_buffer(run_transport_a,
                                       "vegetation",
                                       p_vegetation);
      gpu::helper_bind_optional_buffer(run_transport_a,
                                       "abrasion_mask",
                                       p_abrasion_mask);

      run_transport_a.bind_arguments(nx,
                                     ny,
                                     n_samples,
                                     seed,
                                     0u,
                                     eff_wind_base.x,
                                     eff_wind_base.y,
                                     matter_to_move,
                                     cell_size,
                                     tan_shadow_min,
                                     tan_shadow_max,
                                     eff_max_shadow,
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

      run_transport_b.bind_buffer("bedrock",
                                  run_transport_a.get_buffer("bedrock"));
      run_transport_b.bind_buffer("sediment",
                                  run_stab_ab.get_buffer("sed_out"));
      run_transport_b.bind_buffer("vegetation",
                                  run_transport_a.get_buffer("vegetation"));
      run_transport_b.bind_buffer("abrasion_mask",
                                  run_transport_a.get_buffer("abrasion_mask"));

      run_transport_b.bind_arguments(nx,
                                     ny,
                                     n_samples,
                                     seed,
                                     0u,
                                     eff_wind_base.x,
                                     eff_wind_base.y,
                                     matter_to_move,
                                     cell_size,
                                     tan_shadow_min,
                                     tan_shadow_max,
                                     eff_max_shadow,
                                     abrasion_rate,
                                     p_vegetation ? 1 : 0,
                                     p_abrasion_mask ? 1 : 0);

      // run_stab_ba: repose stabilization reading sed_b and writing sed_a
      auto run_stab_ba = clwrapper::Run("aeolian_stabilize_sediment",
                                        run_transport_a.get_queue());

      run_stab_ba.bind_buffer("bedrock", run_transport_a.get_buffer("bedrock"));
      run_stab_ba.bind_buffer("sed_in", run_stab_ab.get_buffer("sed_out"));
      run_stab_ba.bind_buffer("sed_out",
                              run_transport_a.get_buffer("sediment"));

      run_stab_ba.bind_arguments(nx, ny, cell_size, tan_repose, repose_rate);

      for (int step = step_start; step < step_end; ++step)
      {
        int sub_step = step - step_start;

        if (sub_step % 2 == 0)
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

      run_stab_ab.finish();

      if (abrasion_rate > 0.f)
      {
        run_transport_a.read_buffer("bedrock");
      }

      if (steps_in_batch % 2 == 1)
      {
        run_stab_ab.read_buffer("sed_out");
        sed = sed_temp;
      }
      else if (steps_in_batch > 0)
      {
        run_transport_a.read_buffer("sediment");
      }
    }

    if (step_end < steps)
    {
      comm.sync_overlap_buffers(sed, region, SyncOperation::SmoothBlend);
      if (abrasion_rate > 0.f)
      {
        comm.sync_overlap_buffers(bed, region, SyncOperation::SmoothBlend);
      }
    }
  }

  // Final boundary synchronization
  comm.sync_overlap_buffers(sed, region, SyncOperation::SmoothBlend);
  if (abrasion_rate > 0.f)
  {
    comm.sync_overlap_buffers(bed, region, SyncOperation::SmoothBlend);
  }

  // Update total elevation
  z = bed + sed;

  if (p_sediment_map) *p_sediment_map = sed;
  if (p_bedrock_map) *p_bedrock_map = bed;
}

} // namespace hmap
