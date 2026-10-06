/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cstdint>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/blending.hpp"
#include "highmap/boundary.hpp"
#include "highmap/erosion/sand_erosion.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/operator.hpp"
#include "highmap/range.hpp"

namespace hmap::gpu
{

void sand_dune(Array                &z,
               int                   nparticles,
               const SandDuneParams &params,
               const Array          *p_bedrock,
               const Array          *p_mask,
               Array                *p_deposition_map,
               Array                *p_erosion_map,
               int                   iterations)
{
  if (!validate_non_empty(z)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;

  if (p_mask)
  {
    apply_with_mask(z,
                    p_mask,
                    [&](Array &a)
                    {
                      hmap::gpu::sand_dune(a,
                                           nparticles,
                                           params,
                                           p_bedrock,
                                           nullptr,
                                           p_deposition_map,
                                           p_erosion_map,
                                           iterations);
                    });
    return;
  }

  const glm::ivec2 shape = z.shape;
  Array            z_bckp = Array();
  if ((p_erosion_map != nullptr) || (p_deposition_map != nullptr)) z_bckp = z;

  iterations = std::max(1, iterations);
  int particles_per_pass = std::max(1, nparticles / iterations);

  std::vector<float> z_buf(z.vector.size());

  for (int iter = 0; iter < iterations; ++iter)
  {
    int pass_particles = (iter == iterations - 1)
                             ? nparticles - iter * particles_per_pass
                             : particles_per_pass;
    if (pass_particles <= 0) break;

    // --- Particle transport step

    std::uint32_t pass_seed = params.seed +
                              static_cast<std::uint32_t>(iter) * 10007u;

    auto run_dune = clwrapper::Run("sand_dune_step");

    run_dune.bind_buffer<float>("z", z.vector);
    helper_bind_optional_buffer(run_dune, "bedrock", p_bedrock);

    run_dune.bind_arguments(shape.x,
                            shape.y,
                            pass_particles,
                            pass_seed,
                            params.wind_angle,
                            params.wind_speed,
                            params.crest_speedup,
                            params.hop_length,
                            params.elevation_hop_factor,
                            params.jitter_angle,
                            params.p_sand,
                            params.p_bare,
                            params.shadow_talus,
                            params.slab_height,
                            params.max_hops,
                            params.periodic_boundary ? 1 : 0,
                            params.sand_inflow_rate,
                            p_bedrock ? 1 : 0,
                            iter);

    run_dune.write_buffer("z");
    run_dune.execute(pass_particles);
    run_dune.read_buffer("z");

    // --- Avalanche relaxation pass

    if (params.collapse_rate > 0.f && params.talus > 0.f)
    {
      auto run_avalanche_ab = clwrapper::Run("sand_avalanche_step");
      run_avalanche_ab.bind_buffer<float>("z_in", z.vector);
      run_avalanche_ab.bind_buffer<float>("z_out", z_buf);
      helper_bind_optional_buffer(run_avalanche_ab, "bedrock", p_bedrock);

      run_avalanche_ab.bind_arguments(shape.x,
                                      shape.y,
                                      params.talus,
                                      params.collapse_rate,
                                      p_bedrock ? 1 : 0,
                                      params.periodic_boundary ? 1 : 0);

      run_avalanche_ab.write_buffer("z_in");

      auto run_avalanche_ba =
          clwrapper::Run("sand_avalanche_step", run_avalanche_ab.get_queue());
      run_avalanche_ba.bind_buffer("z_in",
                                   run_avalanche_ab.get_buffer("z_out"));
      run_avalanche_ba.bind_buffer("z_out",
                                   run_avalanche_ab.get_buffer("z_in"));
      helper_bind_optional_buffer(run_avalanche_ba, "bedrock", p_bedrock);

      run_avalanche_ba.bind_arguments(shape.x,
                                      shape.y,
                                      params.talus,
                                      params.collapse_rate,
                                      p_bedrock ? 1 : 0,
                                      params.periodic_boundary ? 1 : 0);

      int avalanche_substeps = 2;
      for (int sub = 0; sub < avalanche_substeps; ++sub)
      {
        if (sub % 2 == 0)
          run_avalanche_ab.execute({shape.x, shape.y});
        else
          run_avalanche_ba.execute({shape.x, shape.y});
      }

      if (avalanche_substeps % 2 == 1)
      {
        run_avalanche_ab.read_buffer("z_out");
        z.vector = z_buf;
      }
      else
      {
        run_avalanche_ab.read_buffer("z_in");
      }
    }
  }

  if (!params.periodic_boundary)
  {
    extrapolate_borders(z);
  }

  // --- Splatmaps generation

  if (p_erosion_map)
  {
    *p_erosion_map = z_bckp - z;
    hmap::clamp_min(*p_erosion_map, 0.f);
  }

  if (p_deposition_map)
  {
    *p_deposition_map = z - z_bckp;
    hmap::clamp_min(*p_deposition_map, 0.f);
  }
}

} // namespace hmap::gpu

