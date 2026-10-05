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
                                  0.f);
  z = hmap::bulkify(z, hmap::PrimitiveType::PRIM_CUBIC_PULSE, 2.f);
  hmap::hydraulic_mise(z);
  hmap::remap(z);

  // --- 2. Build Rock Emission Source Map

  // compute slope magnitude to identify steep cliff rockfall sources
  hmap::Array slope = hmap::gradient_norm(z);

  // emission originates from steep upper ridges and cliffs
  hmap::Array cliff_mask = hmap::select_range(slope,
                                              1.5f / shape.x,
                                              4.0f / shape.x,
                                              0.5f / shape.x);
  hmap::Array high_elevation = hmap::threshold_smooth(z, 0.f, 0.7f);
  hmap::Array emission_density = high_elevation;
  // emission_density = 1.f;

  // --- 3. Configure GPU Rock Rolling Simulation Parameters

  hmap::RockDistribution dist(0u,
                              0.0005f,
                              0.005f,
                              2.0f); // Pareto power-law alpha = 2.0

  hmap::RockSimulationOptions options;
  options.bbox = {0.f, 1.f, 0.f, 1.f};
  options.time_step = 0.005f;
  options.max_steps = int(16.f * 128);
  options.gravity = 9.81f;
  options.soil_friction = 0.3f;
  options.rolling_resistance = 0.05f;
  options.min_velocity = 0.f;
  options.seed = static_cast<uint32_t>(seed);
  options.respawn_out_of_bounds = false;
  options.spawn_fraction = 0.9f;

  // --- 4. Simulate Rock Emission and Downhill Trajectories on GPU

  size_t target_count = 5000;
  std::cout << "Simulating " << target_count << " rocks on GPU...\n";

  hmap::RockField rocks = hmap::simulate_rock_emission(target_count,
                                                       emission_density,
                                                       z,
                                                       nullptr,
                                                       dist,
                                                       options);

  std::cout << "\n=== Simulated Rock Field ===\n";
  std::cout << rocks.to_string() << "\n\n";

  rocks.prune_collisions();

  std::cout << "\n=== Simulated Rock Field after collusion pruning ===\n";
  std::cout << rocks.to_string() << "\n\n";

  // --- 5. Export Visualization

  rocks.to_png("rock_simulation_density.png", shape, z);
  rocks.to_csv("rock_simulation.csv");

  std::cout
      << "Exported 'rock_simulation_density.png' and 'rock_simulation.csv'.\n";
  return 0;
}
