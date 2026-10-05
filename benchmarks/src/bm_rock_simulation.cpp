#include "highmap.hpp"

#include <benchmark/benchmark.h>

using namespace hmap;

// ------------------------------------------------------------
// Args helper
// ------------------------------------------------------------

static void rock_sim_args(benchmark::internal::Benchmark *b)
{
  std::vector<int> counts = {500, 2000, 10000, 20000};
  for (int n : counts)
    b->Args({n});
}

static void BM_rock_simulation_physics(benchmark::State &state)
{
  int rock_count = state.range(0);

  // Setup terrain elevation heightmap
  glm::ivec2 shape = {512, 512};
  Array      elevation = white(shape, 0.0f, 1.0f, 42);

  Array emission_density = elevation;

  RockDistribution dist(0u, 0.0005f, 0.005f, 4.0f);

  RockSeedingOptions seed_opts;
  seed_opts.bbox = {0.f, 1.f, 0.f, 1.f};
  seed_opts.distribution = dist;
  seed_opts.seed = 42;

  RockField rocks = seed_rock_field(rock_count,
                                    emission_density,
                                    {},
                                    seed_opts);

  RockSimulationOptions sim_opts;
  sim_opts.bbox = {0.f, 1.f, 0.f, 1.f};
  sim_opts.time_step = 0.005f;
  sim_opts.max_steps = 500;
  sim_opts.gravity = 9.81f;
  sim_opts.soil_friction = 0.55f;
  sim_opts.rolling_resistance = 0.05f;
  sim_opts.min_velocity = 0.0f;
  sim_opts.seed = 42;
  sim_opts.respawn_out_of_bounds = true;
  sim_opts.spawn_fraction = 0.5f;
  sim_opts.inter_rock_restitution = 0.9f;

  for (auto _ : state)
  {
    RockField current_rocks = rocks;
    RockField result = simulate_rock_trajectories(current_rocks,
                                                  elevation,
                                                  nullptr,
                                                  sim_opts);
    benchmark::DoNotOptimize(result);
  }
}
BENCHMARK(BM_rock_simulation_physics)
    ->Apply(rock_sim_args)
    ->UseRealTime()
    ->Unit(benchmark::kMillisecond);
