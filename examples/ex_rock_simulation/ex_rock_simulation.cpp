#include <iostream>

#include "highmap.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {512, 512};
  glm::vec2  kw = {3.f, 3.f};
  int        seed = 42;

  // --- 1. Generate Terrain Elevation

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::SIMPLEX2,
                                  shape,
                                  kw,
                                  seed,
                                  4);
  z = hmap::bulkify(z, hmap::PrimitiveType::PRIM_CUBIC_PULSE, 2.f);
  hmap::remap(z);

  // --- 2. Build Rock Emission Source Map

  // compute slope magnitude to identify steep cliff rockfall sources
  hmap::Array slope = hmap::gradient_norm(z);

  // emission originates from steep upper ridges and cliffs
  hmap::Array cliff_mask = hmap::select_range(slope,
                                              1.5f / shape.x,
                                              4.0f / shape.x,
                                              0.5f / shape.x);
  hmap::Array high_elevation = hmap::threshold_smooth(z, 0.4f, 0.7f);
  hmap::Array emission_density = cliff_mask * high_elevation;

  // --- 3. Configure GPU Rock Rolling Simulation Parameters

  hmap::RockDistribution dist(0u,
                              0.001f,
                              0.005f,
                              2.0f); // Pareto power-law alpha = 2.0

  hmap::RockSimulationOptions options;
  options.bbox = {0.f, 1.f, 0.f, 1.f};
  options.time_step = 0.005f;
  options.max_steps = int(2.5f * 128);
  options.sub_steps = 32;
  options.gravity = 9.81f;
  options.rolling_resistance = 0.5f;
  options.min_velocity = 0.01f;
  options.seed = static_cast<uint32_t>(seed);
  options.respawn_out_of_bounds = false;

  // --- 4. Simulate Rock Emission and Downhill Trajectories on GPU

  size_t target_count = 2000;
  std::cout << "Simulating " << target_count << " rocks on GPU...\n";

  hmap::RockField rocks = hmap::simulate_rock_emission(target_count,
                                                       emission_density,
                                                       z,
                                                       nullptr,
                                                       dist,
                                                       options);

  std::cout << "\n=== Simulated Rock Field ===\n";
  std::cout << rocks.to_string() << "\n\n";

  // --- 5. Export Visualization

  rocks.to_png("rock_simulation_density.png", shape, z);
  rocks.to_csv("rock_simulation.csv");

  std::cout
      << "Exported 'rock_simulation_density.png' and 'rock_simulation.csv'.\n";
  return 0;
}
