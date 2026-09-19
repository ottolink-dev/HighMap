/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
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

void hydraulic_particle(Array        &z,
                        int           nparticles,
                        std::uint32_t seed,
                        const Array  *p_bedrock,
                        const Array  *p_moisture_map,
                        const Array  *p_elevation_shift,
                        Array        *p_erosion_map,
                        Array        *p_deposition_map,
                        float         c_capacity,
                        float         c_erosion,
                        float         c_deposition,
                        float         c_inertia,
                        float         c_gravity,
                        float         drag_rate,
                        float         evap_rate,
                        float         talus_slope,
                        float         collapse_rate,
                        int           iterations)
{
  if (!validate_non_empty(z)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;
  if (p_moisture_map && !validate_same_shape(z, *p_moisture_map)) return;
  if (p_elevation_shift && !validate_same_shape(z, *p_elevation_shift)) return;

  const glm::ivec2 shape = z.shape;

  Array z_bckp = Array();
  if ((p_erosion_map != nullptr) || (p_deposition_map != nullptr)) z_bckp = z;

  // B1: multi-pass — split particles across iterations so each pass
  // sees the terrain modifications of the previous one, allowing
  // later particles to follow channels carved by earlier ones
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

    auto run = clwrapper::Run("hydraulic_particle");

    run.bind_buffer<float>("z", z.vector);
    helper_bind_optional_buffer(run, "bedrock", p_bedrock);
    helper_bind_optional_buffer(run, "moisture_map", p_moisture_map);
    helper_bind_optional_buffer(run, "elevation_shift", p_elevation_shift);

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
                       p_bedrock ? 1 : 0,
                       p_moisture_map ? 1 : 0,
                       p_elevation_shift ? 1 : 0);

    run.write_buffer("z");
    run.execute(pass_particles);
    run.read_buffer("z");
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
}

void hydraulic_particle(Array        &z,
                        const Array  *p_mask,
                        int           nparticles,
                        std::uint32_t seed,
                        const Array  *p_bedrock,
                        const Array  *p_moisture_map,
                        const Array  *p_elevation_shift,
                        Array        *p_erosion_map,
                        Array        *p_deposition_map,
                        float         c_capacity,
                        float         c_erosion,
                        float         c_deposition,
                        float         c_inertia,
                        float         c_gravity,
                        float         drag_rate,
                        float         evap_rate,
                        float         talus_slope,
                        float         collapse_rate,
                        int           iterations)
{
  apply_with_mask(z,
                  p_mask,
                  [&](Array &a)
                  {
                    gpu::hydraulic_particle(a,
                                            nparticles,
                                            seed,
                                            p_bedrock,
                                            p_moisture_map,
                                            p_elevation_shift,
                                            p_erosion_map,
                                            p_deposition_map,
                                            c_capacity,
                                            c_erosion,
                                            c_deposition,
                                            c_inertia,
                                            c_gravity,
                                            drag_rate,
                                            evap_rate,
                                            talus_slope,
                                            collapse_rate,
                                            iterations);
                  });
}

void hydraulic_particle_multiscale(Array                  &z,
                                   std::uint32_t           seed,
                                   const std::vector<int> &steps_per_level,
                                   const Array            *p_bedrock,
                                   const Array            *p_moisture_map,
                                   const Array            *p_elevation_shift,
                                   Array                  *p_erosion_map,
                                   Array                  *p_deposition_map,
                                   float                   particles_ratio,
                                   float                   c_capacity,
                                   float                   c_erosion,
                                   float                   c_deposition,
                                   float                   c_inertia,
                                   float                   c_gravity,
                                   float                   drag_rate,
                                   float                   evap_rate,
                                   float                   talus_slope,
                                   float                   collapse_rate,
                                   float                   mix)
{
  if (!validate_non_empty(z)) return;

  int nlevels = static_cast<int>(steps_per_level.size());
  if (nlevels == 0) return;

  Array z_orig = z;
  Array z_bckp = Array();
  if ((p_erosion_map != nullptr) || (p_deposition_map != nullptr)) z_bckp = z;

  // Build halving pyramid ladder (coarsest first, final level == z.shape)
  std::vector<glm::ivec2> ladder(nlevels);
  for (int i = 0; i < nlevels; ++i)
  {
    int shift = nlevels - 1 - i;
    ladder[i] = {std::max(2, z.shape.x >> shift),
                 std::max(2, z.shape.y >> shift)};
  }

  Array current_z = z_orig.resample_to_shape(ladder[0]);

  for (int i = 0; i < nlevels; ++i)
  {
    if (i > 0)
    {
      current_z = current_z.resample_to_shape(ladder[i]);

      // Blend with the original input heightmap resampled at the current level
      // resolution to preserve high-frequency structures
      if (mix < 1.f)
      {
        Array z_input_level = z_orig.resample_to_shape(ladder[i]);
        current_z = hmap::lerp(z_input_level,
                               current_z,
                               std::clamp(mix, 0.f, 1.f));
      }
    }

    int level_particles = static_cast<int>(particles_ratio * ladder[i].x *
                                           ladder[i].y);
    int level_iterations = std::max(1, steps_per_level[i]);

    // Resample optional input maps to level shape if provided
    Array        level_bedrock, level_moisture, level_shift;
    const Array *p_lvl_bedrock = nullptr;
    const Array *p_lvl_moisture = nullptr;
    const Array *p_lvl_shift = nullptr;

    if (p_bedrock)
    {
      level_bedrock = p_bedrock->resample_to_shape(ladder[i]);
      p_lvl_bedrock = &level_bedrock;
    }
    if (p_moisture_map)
    {
      level_moisture = p_moisture_map->resample_to_shape(ladder[i]);
      p_lvl_moisture = &level_moisture;
    }
    if (p_elevation_shift)
    {
      level_shift = p_elevation_shift->resample_to_shape(ladder[i]);
      p_lvl_shift = &level_shift;
    }

    gpu::hydraulic_particle(current_z,
                            level_particles,
                            seed + static_cast<std::uint32_t>(i * 1000),
                            p_lvl_bedrock,
                            p_lvl_moisture,
                            p_lvl_shift,
                            /* p_erosion_map */ nullptr,
                            /* p_deposition_map */ nullptr,
                            c_capacity,
                            c_erosion,
                            c_deposition,
                            c_inertia,
                            c_gravity,
                            drag_rate,
                            evap_rate,
                            talus_slope,
                            collapse_rate,
                            level_iterations);
  }

  z = current_z;

  // Splatmaps at final full resolution
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
}

void hydraulic_particle_multiscale(Array                  &z,
                                   const Array            *p_mask,
                                   std::uint32_t           seed,
                                   const std::vector<int> &steps_per_level,
                                   const Array            *p_bedrock,
                                   const Array            *p_moisture_map,
                                   const Array            *p_elevation_shift,
                                   Array                  *p_erosion_map,
                                   Array                  *p_deposition_map,
                                   float                   particles_ratio,
                                   float                   c_capacity,
                                   float                   c_erosion,
                                   float                   c_deposition,
                                   float                   c_inertia,
                                   float                   c_gravity,
                                   float                   drag_rate,
                                   float                   evap_rate,
                                   float                   talus_slope,
                                   float                   collapse_rate,
                                   float                   mix)
{
  apply_with_mask(z,
                  p_mask,
                  [&](Array &a)
                  {
                    gpu::hydraulic_particle_multiscale(a,
                                                       seed,
                                                       steps_per_level,
                                                       p_bedrock,
                                                       p_moisture_map,
                                                       p_elevation_shift,
                                                       p_erosion_map,
                                                       p_deposition_map,
                                                       particles_ratio,
                                                       c_capacity,
                                                       c_erosion,
                                                       c_deposition,
                                                       c_inertia,
                                                       c_gravity,
                                                       drag_rate,
                                                       evap_rate,
                                                       talus_slope,
                                                       collapse_rate,
                                                       mix);
                  });
}

} // namespace hmap::gpu
