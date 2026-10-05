/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/rocks/rock_seeding.hpp"
#include "highmap/rocks/rock_simulation.hpp"

namespace hmap
{

// ============================================================================
//  Rock Simulation Functions (Alphabetically Sorted)
// ============================================================================

RockField simulate_rock_emission(size_t                       rock_count,
                                 const Array                 &emission_density,
                                 const Array                 &elevation,
                                 const Array                 *p_friction_map,
                                 const RockDistribution      &dist,
                                 const RockSimulationOptions &options)
{
  if (rock_count == 0 || !validate_non_empty(emission_density) ||
      !validate_non_empty(elevation))
  {
    return RockField();
  }

  RockSeedingOptions seed_opts;
  seed_opts.bbox = options.bbox;
  seed_opts.distribution = dist;
  seed_opts.seed = options.seed;

  // Seed initial rocks at emission sources
  RockField initial_rocks = seed_rock_field(rock_count,
                                            emission_density,
                                            {},
                                            seed_opts);

  return simulate_rock_trajectories(initial_rocks,
                                    elevation,
                                    p_friction_map,
                                    options);
}

RockField simulate_rock_trajectories(const RockField &rocks,
                                     const Array     &elevation,
                                     const Array     *p_friction_map,
                                     const RockSimulationOptions &options)
{
  if (rocks.empty() || !validate_non_empty(elevation)) return rocks;
  if (p_friction_map && !validate_same_shape(elevation, *p_friction_map))
  {
    return rocks;
  }

  int        num_rocks = static_cast<int>(rocks.size());
  glm::ivec2 shape = elevation.shape;
  glm::vec4  bbox = options.bbox;

  // --- Prepare GPU Rock Buffers (SoA / Aligned Structs)

  std::vector<float>    pos_rad(num_rocks * 4);
  std::vector<float>    vel_mass(num_rocks * 4);
  std::vector<uint32_t> status_class(num_rocks * 2);

  float max_radius = 1e-4f;

  for (int i = 0; i < num_rocks; ++i)
  {
    const auto &item = rocks[i];
    float       r = std::max(1e-5f, item.radius);
    max_radius = std::max(max_radius, r);

    // Initial position: if z is zero, place on terrain surface
    float z_pos = item.position.z;
    if (std::abs(z_pos) < 1e-5f)
    {
      float u = (item.position.x - bbox.x) / (bbox.y - bbox.x);
      float v = (item.position.y - bbox.z) / (bbox.w - bbox.z);
      float fx = std::clamp(u * (shape.x - 1),
                            0.0f,
                            static_cast<float>(shape.x - 1));
      float fy = std::clamp(v * (shape.y - 1),
                            0.0f,
                            static_cast<float>(shape.y - 1));
      int   ix = std::clamp(static_cast<int>(fx), 0, std::max(0, shape.x - 2));
      int   iy = std::clamp(static_cast<int>(fy), 0, std::max(0, shape.y - 2));
      float tx = fx - ix;
      float ty = fy - iy;
      float h00 = elevation(ix, iy);
      float h10 = elevation(ix + 1, iy);
      float h01 = elevation(ix, iy + 1);
      float h11 = elevation(ix + 1, iy + 1);
      float h = (1.0f - tx) * (1.0f - ty) * h00 + tx * (1.0f - ty) * h10 +
                (1.0f - tx) * ty * h01 + tx * ty * h11;
      z_pos = h + r;
    }

    pos_rad[i * 4 + 0] = item.position.x;
    pos_rad[i * 4 + 1] = item.position.y;
    pos_rad[i * 4 + 2] = z_pos;
    pos_rad[i * 4 + 3] = r;

    // Mass proportional to volume (4/3 * pi * r^3, rho ~ 2500 kg/m^3)
    float mass = (4.0f / 3.0f) * static_cast<float>(M_PI) * r * r * r * 2500.0f;

    vel_mass[i * 4 + 0] = 0.0f;
    vel_mass[i * 4 + 1] = 0.0f;
    vel_mass[i * 4 + 2] = 0.0f;
    vel_mass[i * 4 + 3] = mass;

    status_class[i * 2 + 0] = 0; // ROCK_STATUS_ACTIVE
    status_class[i * 2 + 1] = item.class_id;
  }

  std::vector<float> init_pos_rad = pos_rad;

  // --- Progressive Spawning Step Calculation
  int total_steps = std::max(1, options.max_steps);

  std::vector<int> spawn_step(num_rocks, 0);
  float            spawn_frac = std::clamp(options.spawn_fraction, 0.0f, 1.0f);
  int              max_spawn_step = static_cast<int>(spawn_frac * total_steps);

  if (max_spawn_step > 0 && num_rocks > 1)
  {
    for (int i = 0; i < num_rocks; ++i)
    {
      // Distribute spawn step linearly / randomly across rock indices
      spawn_step[i] = static_cast<int>(
          (static_cast<float>(i) / (num_rocks - 1)) * max_spawn_step);
      if (spawn_step[i] > 0)
      {
        status_class[i * 2 + 0] = 3; // ROCK_STATUS_UNSPAWNED
      }
    }
  }

  // --- OpenCL Simulation Run Dispatch
  // Configure spatial hash grid based on domain bbox and rock radii
  float dom_w = bbox.y - bbox.x;
  float dom_h = bbox.w - bbox.z;
  float cell_size = std::max(2.0f * max_radius, 1e-4f);

  int   grid_cells_x = std::clamp(static_cast<int>(dom_w / cell_size), 4, 512);
  int   grid_cells_y = std::clamp(static_cast<int>(dom_h / cell_size), 4, 512);
  float cell_size_x = dom_w / static_cast<float>(grid_cells_x);
  float cell_size_y = dom_h / static_cast<float>(grid_cells_y);

  int              total_cells = grid_cells_x * grid_cells_y;
  std::vector<int> grid_heads(total_cells, -1);
  std::vector<int> grid_next(num_rocks, -1);
  std::vector<int> grid_step_tag(total_cells, 0);

  auto run_sim = clwrapper::Run("rock_simulate_physics");
  run_sim.bind_buffer<float>("pos_rad", pos_rad);
  run_sim.bind_buffer<float>("vel_mass", vel_mass);
  run_sim.bind_buffer<uint32_t>("status_class", status_class);
  run_sim.bind_buffer<float>("init_pos_rad", init_pos_rad);
  run_sim.bind_buffer<int>("spawn_step", spawn_step);
  run_sim.bind_buffer<float>("z", elevation.vector);
  gpu::helper_bind_optional_buffer(run_sim, "friction_map", p_friction_map);
  run_sim.bind_buffer<int>("grid_heads", grid_heads);
  run_sim.bind_buffer<int>("grid_next", grid_next);
  run_sim.bind_buffer<int>("grid_step_tag", grid_step_tag);

  run_sim.bind_arguments(num_rocks,
                         shape.x,
                         shape.y,
                         bbox,
                         options.time_step,
                         options.gravity,
                         options.soil_friction,
                         options.rolling_resistance,
                         options.inter_rock_restitution,
                         options.min_velocity,
                         p_friction_map ? 1 : 0,
                         options.respawn_out_of_bounds ? 1 : 0,
                         options.seed,
                         total_steps,
                         grid_cells_x,
                         grid_cells_y,
                         cell_size_x,
                         cell_size_y);

  run_sim.write_buffer("pos_rad");
  run_sim.write_buffer("vel_mass");
  run_sim.write_buffer("status_class");
  run_sim.write_buffer("init_pos_rad");
  run_sim.write_buffer("spawn_step");
  run_sim.write_buffer("z");
  if (p_friction_map) run_sim.write_buffer("friction_map");
  run_sim.write_buffer("grid_heads");
  run_sim.write_buffer("grid_next");
  run_sim.write_buffer("grid_step_tag");

  run_sim.execute_async(num_rocks);
  run_sim.finish();

  run_sim.read_buffer("pos_rad");
  run_sim.read_buffer("status_class");

  std::vector<Rock> final_rocks;
  final_rocks.reserve(num_rocks);

  for (int i = 0; i < num_rocks; ++i)
  {
    uint32_t status = status_class[i * 2 + 0];
    if (status == 2) // ROCK_STATUS_OUT_OF_BOUNDS
    {
      if (options.respawn_out_of_bounds)
      {
        final_rocks.emplace_back(rocks[i]);
      }
      continue;
    }

    glm::vec3 pos(pos_rad[i * 4 + 0], pos_rad[i * 4 + 1], pos_rad[i * 4 + 2]);
    float     r = pos_rad[i * 4 + 3];
    uint32_t  cid = status_class[i * 2 + 1];

    final_rocks.emplace_back(pos, cid, r);
  }

  return RockField(std::move(final_rocks));
}

} // namespace hmap
