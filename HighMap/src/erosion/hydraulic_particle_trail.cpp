/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cstdint>
#include <memory>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/erosion.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/operator.hpp"
#include "highmap/range.hpp"

namespace hmap::gpu
{

void hydraulic_particle_trail(Array        &z,
                              int           nparticles,
                              std::uint32_t seed,
                              const Array  *p_bedrock,
                              const Array  *p_moisture_map,
                              const Array  *p_elevation_shift,
                              Array        *p_erosion_map,
                              Array        *p_deposition_map,
                              Array        *p_trail_map,
                              float         c_capacity,
                              float         c_erosion,
                              float         c_deposition,
                              float         c_inertia,
                              float         c_gravity,
                              float         drag_rate,
                              float         evap_rate,
                              float         talus_slope,
                              float         collapse_rate,
                              float         c_trail_deposit,
                              float         c_trail_attraction,
                              float         trail_evap_rate,
                              int           iterations)
{
  if (!validate_non_empty(z)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;
  if (p_moisture_map && !validate_same_shape(z, *p_moisture_map)) return;
  if (p_elevation_shift && !validate_same_shape(z, *p_elevation_shift)) return;

  const glm::ivec2 shape = z.shape;

  Array z_bckp = Array();
  if ((p_erosion_map != nullptr) || (p_deposition_map != nullptr)) z_bckp = z;

  Array trail(shape, 0.f);
  if (p_trail_map && p_trail_map->shape == shape)
  {
    trail = *p_trail_map;
  }

  iterations = std::max(1, iterations);
  int particles_per_pass = std::max(1, nparticles / iterations);

  for (int iter = 0; iter < iterations; ++iter)
  {
    // last pass gets the remainder so total == nparticles
    int pass_particles = (iter == iterations - 1)
                             ? nparticles - iter * particles_per_pass
                             : particles_per_pass;
    if (pass_particles <= 0) break;

    // vary seed per pass to avoid repeating the same particle paths
    std::uint32_t pass_seed = seed +
                              static_cast<std::uint32_t>(iter) *
                                  static_cast<std::uint32_t>(pass_particles) *
                                  2u;

    auto run = clwrapper::Run("hydraulic_particle_trail");

    run.bind_buffer<float>("z", z.vector);
    helper_bind_optional_buffer(run, "bedrock", p_bedrock);
    helper_bind_optional_buffer(run, "moisture_map", p_moisture_map);
    helper_bind_optional_buffer(run, "elevation_shift", p_elevation_shift);
    run.bind_buffer<float>("trail", trail.vector);

    float cell_talus = talus_slope / static_cast<float>(shape.x);

    run.bind_arguments(shape.x,
                       shape.y,
                       pass_particles,
                       pass_seed,
                       c_capacity,
                       c_erosion,
                       c_deposition,
                       c_inertia,
                       c_gravity,
                       drag_rate,
                       evap_rate,
                       cell_talus,
                       collapse_rate,
                       c_trail_deposit,
                       c_trail_attraction,
                       p_bedrock ? 1 : 0,
                       p_moisture_map ? 1 : 0,
                       p_elevation_shift ? 1 : 0);

    run.write_buffer("z");
    run.write_buffer("trail");
    run.execute(pass_particles);
    run.read_buffer("z");
    run.read_buffer("trail");

    // decay trail between passes
    if (iter < iterations - 1 && trail_evap_rate > 0.f)
    {
      trail *= (1.f - std::clamp(trail_evap_rate, 0.f, 1.f));
    }
  }

  extrapolate_borders(z);

  // splatmaps
  if (p_erosion_map)
  {
    *p_erosion_map = z_bckp - z;
    clamp_min(*p_erosion_map, 0.f);
  }

  if (p_deposition_map)
  {
    *p_deposition_map = z - z_bckp;
    clamp_min(*p_deposition_map, 0.f);
  }

  if (p_trail_map)
  {
    *p_trail_map = trail;
  }
}

void hydraulic_particle_trail(Array        &z,
                              const Array  *p_mask,
                              int           nparticles,
                              std::uint32_t seed,
                              const Array  *p_bedrock,
                              const Array  *p_moisture_map,
                              const Array  *p_elevation_shift,
                              Array        *p_erosion_map,
                              Array        *p_deposition_map,
                              Array        *p_trail_map,
                              float         c_capacity,
                              float         c_erosion,
                              float         c_deposition,
                              float         c_inertia,
                              float         c_gravity,
                              float         drag_rate,
                              float         evap_rate,
                              float         talus_slope,
                              float         collapse_rate,
                              float         c_trail_deposit,
                              float         c_trail_attraction,
                              float         trail_evap_rate,
                              int           iterations)
{
  apply_with_mask(z,
                  p_mask,
                  [&](Array &a)
                  {
                    gpu::hydraulic_particle_trail(a,
                                                  nparticles,
                                                  seed,
                                                  p_bedrock,
                                                  p_moisture_map,
                                                  p_elevation_shift,
                                                  p_erosion_map,
                                                  p_deposition_map,
                                                  p_trail_map,
                                                  c_capacity,
                                                  c_erosion,
                                                  c_deposition,
                                                  c_inertia,
                                                  c_gravity,
                                                  drag_rate,
                                                  evap_rate,
                                                  talus_slope,
                                                  collapse_rate,
                                                  c_trail_deposit,
                                                  c_trail_attraction,
                                                  trail_evap_rate,
                                                  iterations);
                  });
}

} // namespace hmap::gpu
