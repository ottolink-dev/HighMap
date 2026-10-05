#include <iostream>

#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {512, 512};
  shape = {1024, 1024};
  glm::vec2 kw = {4.f, 4.f};
  int       seed = 42;

  // --- 1. Generate Terrain Elevation

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                  shape,
                                  kw,
                                  seed,
                                  8,
                                  0.7f);
  z = hmap::bulkify(z, hmap::PrimitiveType::PRIM_CUBIC_PULSE, 1.f);

  auto ze = z;
  hmap::hydraulic_mise(ze);
  z = hmap::lerp(z, ze, 0.5f);

  hmap::remap(z);

  // --- 2. Build Rock Emission Source Map

  // compute slope magnitude to identify steep cliff rockfall sources
  hmap::Array slope = hmap::gradient_norm(z);

  // emission originates from upper elevations
  hmap::Array emission_density = hmap::Array(shape, 1.f);
  // emission_density = hmap::threshold_smooth(z, 0.f, 0.7f);

  // --- 3. Configure GPU Rock Rolling Simulation Parameters

  // Pareto power-law alpha
  hmap::RockDistribution dist(0u, 0.0005f, 0.002f, 4.f);

  hmap::RockSimulationOptions options;
  options.bbox = {0.f, 1.f, 0.f, 1.f};
  options.time_step = 0.005f;
  options.max_steps = int(8.f * 128);
  options.gravity = 9.81f;
  options.soil_friction = 0.55f;
  options.rolling_resistance = 0.05f;
  options.min_velocity = 0.f;
  options.seed = static_cast<uint32_t>(seed);
  options.respawn_out_of_bounds = true;
  options.spawn_fraction = 0.8f;
  options.inter_rock_restitution = 0.9f;

  // --- 4. Simulate Rock Emission and Downhill Trajectories on GPU

  size_t target_count = 80000;
  std::cout << "Simulating " << target_count << " rocks on GPU...\n";

  hmap::RockField rocks = hmap::simulate_rock_emission(target_count,
                                                       emission_density,
                                                       z,
                                                       nullptr,
                                                       dist,
                                                       options);

  std::cout << "\n=== Simulated Rock Field ===\n";
  std::cout << rocks.to_string() << "\n\n";

  rocks.resolve_collisions();
  rocks.prune_collisions();

  std::cout << "\n=== Simulated Rock Field after collisions resolution ===\n";
  std::cout << rocks.to_string() << "\n\n";

  // --- 5. Convert Rock Field into Heightmap and Accumulate

  std::cout << "Converting RockField into heightmap on GPU...\n";

  hmap::Array rock_map(shape);
  hmap::Array rock_elevation = rocks.to_heightmap(shape,
                                                  hmap::SCATTER_SHAPE_POLYGON,
                                                  2.f,
                                                  std::nullopt,
                                                  &rock_map);

  // accumulate rock elevation contributions directly onto terrain
  hmap::Array z_with_rocks = z + rock_elevation;

  // --- 6. Export Visualization

  rocks.to_png("rock_simulation_density.png", shape, z);
  rocks.to_csv("rock_simulation.csv");
  hmap::export_banner_png("rock_simulation_heightmap.png",
                          {z, rock_elevation, z_with_rocks},
                          hmap::Cmap::TERRAIN);

  rock_elevation.dump("rock_elevation.png");
  rock_map.dump("rock_map.png");
  z_with_rocks.dump("out.png");

  std::cout << "Exported 'rock_simulation_density.png', 'rock_simulation.csv', "
               "and 'rock_simulation_heightmap.png'.\n";
  return 0;
}
