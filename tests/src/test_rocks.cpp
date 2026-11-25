#include <cmath>
#include <fstream>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/world/rocks.hpp"

#include <gtest/gtest.h>

using namespace hmap;

// --- Helper Functions

static constexpr float eps = 1e-5f;

static bool float_eq(float a, float b, float tol = eps)
{
  return std::abs(a - b) < tol;
}

// --- RockField Unit Tests

TEST(RockTest, BasicRockItem)
{
  Rock r(glm::vec3(1.f, 2.f, 3.f), 5u, 0.2f);
  EXPECT_TRUE(float_eq(r.position.x, 1.f));
  EXPECT_TRUE(float_eq(r.position.y, 2.f));
  EXPECT_TRUE(float_eq(r.position.z, 3.f));
  EXPECT_EQ(r.class_id, 5u);
  EXPECT_TRUE(float_eq(r.radius, 0.2f));
}

TEST(RockFieldTest, PowerLawDistribution)
{
  RockField field;
  for (int i = 0; i < 200; ++i)
  {
    field.push_back(Rock(glm::vec3(float(i), 0.f, 0.f), 0u, 0.01f));
  }

  RockDistribution dist(0u, 0.01f, 1.0f, 2.0f);
  field.apply_power_law_distribution(dist, 1234);

  float  min_r = 1e9f;
  float  max_r = -1e9f;
  size_t small_count = 0;

  for (const auto &r : field)
  {
    min_r = std::min(min_r, r.radius);
    max_r = std::max(max_r, r.radius);
    if (r.radius < 0.2f) small_count++;
  }

  EXPECT_GE(min_r, 0.01f);
  EXPECT_LE(max_r, 1.0f);
  // Heavy-tailed power law: the vast majority should be small rocks
  EXPECT_GT(small_count, 150u);
}

TEST(RockFieldTest, SlopeSorting)
{
  // 10x10 array with slope increasing linearly with x
  Array slope({10, 10}, 0.f);
  for (int j = 0; j < 10; ++j)
    for (int i = 0; i < 10; ++i)
      slope(i, j) = static_cast<float>(i) * 0.1f;

  RockField field;
  // Place rocks along x with identical initial radius
  for (int i = 0; i < 10; ++i)
  {
    field.push_back(Rock(glm::vec3(static_cast<float>(i), 0.5f, 0.f),
                         0u,
                         0.1f * float(i + 1)));
  }

  glm::vec4 bbox = {0.f, 9.f, 0.f, 9.f};
  field.apply_slope_sorting(slope, 1.0f, bbox);

  // Rock at low slope (x=0) should have larger radius than rock at high slope
  // (x=9)
  EXPECT_GT(field.front().radius, field.back().radius);
}

TEST(RockFieldTest, InterstitialPacking)
{
  RockField field;
  // Place two large boulders
  field.push_back(Rock(glm::vec3(0.2f, 0.5f, 0.f), 0u, 0.1f));
  field.push_back(Rock(glm::vec3(0.8f, 0.5f, 0.f), 0u, 0.1f));

  size_t    original_count = field.size();
  glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f};

  // 1. Pack without density (uniform)
  field.pack_interstitial(20, {}, 0.01f, 0.02f, 1u, 42, bbox);

  EXPECT_GT(field.size(), original_count);
  EXPECT_EQ(field[0].radius, 0.1f);
  EXPECT_EQ(field[1].radius, 0.1f);

  // 2. Pack with density map
  Array  density({32, 32}, 1.0f);
  size_t count_before_density = field.size();
  field.pack_interstitial(20, density, 0.005f, 0.01f, 2u, 101, bbox);
  EXPECT_GT(field.size(), count_before_density);
}

TEST(RockFieldTest, Seeding)
{
  Array density({32, 32}, 1.0f);
  Array slope({32, 32}, 0.5f);

  RockSeedingOptions opts;
  opts.seed = 123;
  opts.distribution = RockDistribution(1u, 0.005f, 0.05f, 2.0f);

  RockField field = seed_rock_field(50, density, {}, opts);
  EXPECT_GT(field.size(), 0u);
  EXPECT_LE(field.size(), 50u);
}

TEST(RockFieldTest, VirtualArraySeeding)
{
  glm::ivec2 shape = {64, 64};
  glm::ivec2 tile_shape = {32, 32};
  int        halo = 4;
  glm::vec4  bbox = {0.0f, 100.0f, 0.0f, 100.0f};

  VirtualArray density(shape, tile_shape, halo);
  ComputeMode  cm;
  density.fill(1.0f, cm);

  VirtualArray exclusion;

  RockSeedingOptions options;
  options.seed = 42;
  options.bbox = bbox;
  options.distribution = RockDistribution(1u, 0.5f, 2.0f, 2.0f);

  size_t target_count = 60;

  RockField field = hmap::va::seed_rock_field(target_count,
                                              density,
                                              exclusion,
                                              options,
                                              cm);

  EXPECT_GT(field.size(), 0u);
  EXPECT_LE(field.size(), target_count);

  for (const auto &rock : field)
  {
    EXPECT_GE(rock.position.x, 0.0f);
    EXPECT_LE(rock.position.x, 100.0f);
    EXPECT_GE(rock.position.y, 0.0f);
    EXPECT_LE(rock.position.y, 100.0f);
    EXPECT_GE(rock.radius, 0.5f);
    EXPECT_LE(rock.radius, 2.0f);
  }
}

TEST(RockFieldTest, VirtualArrayExclusionMask)
{
  glm::ivec2 shape = {64, 64};
  glm::ivec2 tile_shape = {32, 32};
  int        halo = 4;
  glm::vec4  bbox = {0.0f, 1.0f, 0.0f, 1.0f};

  VirtualArray density(shape, tile_shape, halo);
  VirtualArray exclusion(shape, tile_shape, halo);
  ComputeMode  cm;
  density.fill(1.0f, cm);
  exclusion.fill(0.0f, cm);

  // Exclude top half (y >= 0.5)
  Array excl_array(shape, 0.0f);
  for (int j = 32; j < 64; ++j)
  {
    for (int i = 0; i < 64; ++i)
    {
      excl_array(i, j) = 1.0f;
    }
  }
  exclusion.from_array(excl_array, cm);

  RockSeedingOptions options;
  options.seed = 123;
  options.bbox = bbox;
  options.exclusion_threshold = 0.5f;

  RockField field = hmap::va::seed_rock_field(60,
                                              density,
                                              exclusion,
                                              options,
                                              cm);

  EXPECT_GT(field.size(), 0u);

  for (const auto &rock : field)
  {
    EXPECT_LT(rock.position.y, 0.55f);
  }
}

TEST(RockFieldTest, VirtualArrayEmptyAndZeroDensity)
{
  VirtualArray empty_va;
  ComputeMode  cm;

  RockField f1 = hmap::va::seed_rock_field(50, empty_va);
  EXPECT_EQ(f1.size(), 0u);

  glm::ivec2   shape = {32, 32};
  glm::ivec2   tile_shape = {16, 16};
  VirtualArray zero_density(shape, tile_shape, 2);
  zero_density.fill(0.0f, cm);

  RockField f2 = hmap::va::seed_rock_field(50, zero_density);
  EXPECT_EQ(f2.size(), 0u);
}

// --- Rock Simulation GPU Tests

TEST(RockSimulationTest, DownhillMotionAndFrictionRest)
{
  // Create a slope in the first half that flattens into a valley
  Array elev(glm::ivec2(64, 64));
  for (int j = 0; j < 64; ++j)
  {
    for (int i = 0; i < 64; ++i)
    {
      float x = static_cast<float>(i) / 63.0f;
      elev(i, j) = (x < 0.4f) ? (0.4f - x) * 1.5f : 0.0f;
    }
  }

  // Release a rock at top of the slope
  RockField field;
  field.push_back(Rock(glm::vec3(0.05f, 0.5f, 0.0f), 0u, 0.02f));

  RockSimulationOptions opts;
  opts.time_step = 0.005f;
  opts.max_steps = 1500;
  opts.soil_friction = 0.2f;
  opts.rolling_resistance = 0.05f;

  field.simulate_physics(elev, nullptr, opts);

  ASSERT_EQ(field.size(), 1u);
  // Rock must have moved downhill and come to rest in the valley
  EXPECT_GT(field[0].position.x, 0.2f);
  EXPECT_LT(field[0].position.x, 0.9f);
  // Rock elevation should follow the terrain surface + radius
  EXPECT_NEAR(field[0].position.z, field[0].radius, 0.01f);
}

TEST(RockSimulationTest, DiscardOutOfBounds)
{
  // Constant steep slope that drives the rock across the boundary
  Array elev(glm::ivec2(64, 64));
  for (int j = 0; j < 64; ++j)
  {
    for (int i = 0; i < 64; ++i)
    {
      elev(i, j) = 1.0f - static_cast<float>(i) / 63.0f;
    }
  }

  RockField field;
  field.push_back(Rock(glm::vec3(0.1f, 0.5f, 0.0f), 0u, 0.02f));

  RockSimulationOptions opts;
  opts.time_step = 0.005f;
  opts.max_steps = 2000;
  opts.soil_friction = 0.1f;

  field.simulate_physics(elev, nullptr, opts);

  // Rock that rolled out of bounds must be discarded
  EXPECT_EQ(field.size(), 0u);
}

TEST(RockSimulationTest, RespawnOutOfBounds)
{
  // Constant steep slope that drives the rock across the boundary
  Array elev(glm::ivec2(64, 64));
  for (int j = 0; j < 64; ++j)
  {
    for (int i = 0; i < 64; ++i)
    {
      elev(i, j) = 1.0f - static_cast<float>(i) / 63.0f;
    }
  }

  RockField field;
  field.push_back(Rock(glm::vec3(0.1f, 0.5f, 0.0f), 0u, 0.02f));
  field.push_back(Rock(glm::vec3(0.2f, 0.5f, 0.0f), 0u, 0.02f));

  RockSimulationOptions opts;
  opts.time_step = 0.005f;
  opts.max_steps = 1500;
  opts.soil_friction = 0.1f;
  opts.respawn_out_of_bounds = true;

  field.simulate_physics(elev, nullptr, opts);

  // When respawn is enabled, rock count remains constant
  EXPECT_EQ(field.size(), 2u);
  for (const auto &rock : field)
  {
    EXPECT_GE(rock.position.x, 0.0f);
    EXPECT_LE(rock.position.x, 1.0f);
    EXPECT_GE(rock.position.y, 0.0f);
    EXPECT_LE(rock.position.y, 1.0f);
  }
}

TEST(RockSimulationTest, EmissionSimulation)
{
  // Terrain high in left quarter, sloping down to valley
  Array elev(glm::ivec2(64, 64));
  Array emission(glm::ivec2(64, 64), 0.0f);

  for (int j = 0; j < 64; ++j)
  {
    for (int i = 0; i < 64; ++i)
    {
      float x = static_cast<float>(i) / 63.0f;
      elev(i, j) = (x < 0.3f) ? (0.3f - x) * 2.0f : 0.0f;
      if (i < 10) emission(i, j) = 1.0f; // Cliff top emission
    }
  }

  RockDistribution      dist(1u, 0.01f, 0.03f, 2.0f);
  RockSimulationOptions opts;
  opts.time_step = 0.005f;
  opts.max_steps = 800;
  opts.rolling_resistance = 0.08f;
  opts.seed = 42;

  RockField result = simulate_rock_emission(20,
                                            emission,
                                            elev,
                                            nullptr,
                                            dist,
                                            opts);

  EXPECT_GT(result.size(), 0u);
  // All simulated resting rocks should stay strictly within the domain
  for (const auto &rock : result)
  {
    EXPECT_GE(rock.position.x, 0.0f);
    EXPECT_LE(rock.position.x, 1.0f);
    EXPECT_GE(rock.position.y, 0.0f);
    EXPECT_LE(rock.position.y, 1.0f);
  }
}

TEST(RockSimulationTest, CollisionNonOverlap)
{
  // Flat ground with 2 rocks initialized overlapping
  Array elev(glm::ivec2(32, 32), 0.0f);

  RockField field;
  // Two rocks of radius 0.05 placed at distance 0.02 (overlapping by 0.08)
  field.push_back(Rock(glm::vec3(0.50f, 0.5f, 0.0f), 0u, 0.05f));
  field.push_back(Rock(glm::vec3(0.52f, 0.5f, 0.0f), 0u, 0.05f));

  RockSimulationOptions opts;
  opts.time_step = 0.005f;
  opts.max_steps = 200;
  opts.soil_friction = 0.5f;

  field.simulate_physics(elev, nullptr, opts);

  ASSERT_EQ(field.size(), 2u);
  float dx = field[0].position.x - field[1].position.x;
  float dy = field[0].position.y - field[1].position.y;
  float dist = std::sqrt(dx * dx + dy * dy);
  float min_allowed_dist = field[0].radius + field[1].radius;

  // Rocks must separate and respect each other's radii (distance >= sum of
  // radii - small tolerance)
  EXPECT_GE(dist, min_allowed_dist - 1e-4f);
}

TEST(RockSimulationTest, ProgressiveSpawn)
{
  Array elev(glm::ivec2(64, 64));
  for (int j = 0; j < 64; ++j)
  {
    for (int i = 0; i < 64; ++i)
    {
      float x = static_cast<float>(i) / 63.0f;
      elev(i, j) = (x < 0.4f) ? (0.4f - x) * 1.5f : 0.0f;
    }
  }

  RockField field;
  for (int i = 0; i < 10; ++i)
  {
    field.push_back(
        Rock(glm::vec3(0.05f, 0.1f + 0.08f * float(i), 0.0f), 0u, 0.02f));
  }

  RockSimulationOptions opts;
  opts.time_step = 0.005f;
  opts.max_steps = 1500;
  opts.soil_friction = 0.2f;
  opts.rolling_resistance = 0.05f;
  opts.spawn_fraction =
      0.5f; // rocks spawn over the first 50% of the simulation

  field.simulate_physics(elev, nullptr, opts);

  ASSERT_EQ(field.size(), 10u);
  for (const auto &rock : field)
  {
    // Verify that every rock (including later-spawned rocks) has moved downhill
    EXPECT_GT(rock.position.x, 0.2f);
    EXPECT_LE(rock.position.x, 1.0f);
    EXPECT_GE(rock.position.y, 0.0f);
    EXPECT_LE(rock.position.y, 1.0f);
  }
}

// --- Procedural Placement

TEST(RockProceduralTest, EmptyInput)
{
  RockField field = place_rocks_procedural(100, Array());
  EXPECT_EQ(field.size(), 0u);

  field = place_rocks_procedural(0, Array(glm::ivec2(16, 16), 1.f));
  EXPECT_EQ(field.size(), 0u);
}

TEST(RockProceduralTest, DeterministicAndBounded)
{
  glm::ivec2 shape = {64, 64};
  Array      z(shape);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      z(i, j) = float(i) / float(shape.x - 1);
    }
  }

  RockSeedingOptions options;
  options.seed = 7;
  options.distribution = RockDistribution(0u, 0.005f, 0.03f);

  RockField a = place_rocks_procedural(300,
                                       z,
                                       0.f,
                                       1.f,
                                       1.f,
                                       {},
                                       {},
                                       0.f,
                                       1.f,
                                       1.f,
                                       options);
  RockField b = place_rocks_procedural(300,
                                       z,
                                       0.f,
                                       1.f,
                                       1.f,
                                       {},
                                       {},
                                       0.f,
                                       1.f,
                                       1.f,
                                       options);

  EXPECT_GT(a.size(), 0u);
  EXPECT_LE(a.size(), 300u);
  EXPECT_EQ(a.size(), b.size());
}

TEST(RockProceduralTest, SizeBroadensAtLowElevation)
{
  glm::ivec2 shape = {64, 64};
  Array      z(shape);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      z(i, j) = float(i) / float(shape.x - 1);
    }
  }

  RockSeedingOptions options;
  options.seed = 3;
  options.distribution = RockDistribution(0u, 0.002f, 0.02f, 1.f);

  RockField field = place_rocks_procedural(400,
                                           z,
                                           0.f,
                                           1.f,
                                           1.f,
                                           {},
                                           {},
                                           0.f,
                                           1.f,
                                           1.f,
                                           options);

  // elevation increases with x: low elevation is x < 0.5
  float max_low = 0.f;
  float max_high = 0.f;
  float min_low = 1e9f;
  int   n_low = 0;
  int   n_high = 0;

  for (const auto &rock : field)
  {
    if (rock.position.x < 0.5f)
    {
      max_low = std::max(max_low, rock.radius);
      min_low = std::min(min_low, rock.radius);
      n_low++;
    }
    else
    {
      max_high = std::max(max_high, rock.radius);
      n_high++;
    }
  }

  ASSERT_GT(n_low, 0);
  ASSERT_GT(n_high, 0);

  // larger rocks at low elevation, but small rocks are still present
  EXPECT_GT(max_low, max_high);
  EXPECT_LT(min_low, 0.004f);
}

// --- Rock Density

TEST(RockDensityTest, EmptyAndInvalidInputs)
{
  Array empty_array;
  Array valid_hmap(glm::ivec2(32, 32), 0.5f);
  Array mismatched_array(glm::ivec2(16, 16), 1.0f);

  Array d1 = build_rock_density(empty_array);
  EXPECT_TRUE(d1.vector.empty());

  Array d2 = build_rock_density(valid_hmap,
                                0.05f,
                                0.01f,
                                1.0f,
                                1.0f,
                                1.0f,
                                0.1f,
                                1000.0f,
                                mismatched_array,
                                1.0f);
  EXPECT_TRUE(d2.vector.empty());

  Array d3 = build_rock_density(valid_hmap,
                                0.05f,
                                0.01f,
                                1.0f,
                                1.0f,
                                1.0f,
                                0.1f,
                                1000.0f,
                                {},
                                1.0f,
                                mismatched_array);
  EXPECT_TRUE(d3.vector.empty());
}

TEST(RockDensityTest, BasicComputationAndExclusion)
{
  glm::ivec2 shape = {32, 32};
  Array      hmap(shape, 0.0f);
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      hmap(i, j) = static_cast<float>(j) / static_cast<float>(shape.y);
    }
  }

  Array density = build_rock_density(hmap);
  ASSERT_EQ(density.shape.x, shape.x);
  ASSERT_EQ(density.shape.y, shape.y);

  for (float v : density.vector)
  {
    EXPECT_GE(v, 0.0f);
    EXPECT_LE(v, 1.0f);
  }

  // Exclusion mask test
  Array mask(shape, 0.0f);
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x / 2; ++i)
    {
      mask(i, j) = 1.0f;
    }
  }

  Array density_masked = build_rock_density(hmap,
                                            0.05f,
                                            0.01f,
                                            1.0f,
                                            1.0f,
                                            1.0f,
                                            0.1f,
                                            1000.0f,
                                            {},
                                            1.0f,
                                            mask,
                                            0.5f);

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x / 2; ++i)
    {
      EXPECT_FLOAT_EQ(density_masked(i, j), 0.0f);
    }
  }
}

TEST(RockProceduralTest, AutomaticDensityAndNoOverlaps)
{
  glm::ivec2 shape = {64, 64};
  Array      z(shape, 0.5f);

  RockSeedingOptions options;
  options.seed = 42;
  options.distribution = RockDistribution(0u, 0.01f, 0.04f, 1.5f);

  // When density is not provided, place_rocks_procedural uses a flat density
  // equal to 1
  RockField field = place_rocks_procedural(200,
                                           z,
                                           0.f,
                                           1.f,
                                           1.f,
                                           {},
                                           {},
                                           0.f,
                                           1.f,
                                           1.f,
                                           options);

  EXPECT_GT(field.size(), 0u);

  field.prune_collisions();

  // Check that collision pruning ensured rocks do not significantly overlap
  for (size_t i = 0; i < field.size(); ++i)
  {
    for (size_t j = i + 1; j < field.size(); ++j)
    {
      float dist = glm::distance(glm::vec2(field[i].position),
                                 glm::vec2(field[j].position));
      EXPECT_GE(dist, field[i].radius + field[j].radius - 1e-4f);
    }
  }
}

TEST(RockProceduralTest, VirtualArrayBasicExecution)
{
  glm::ivec2 shape = {64, 64};
  glm::ivec2 tile_shape = {32, 32};
  int        halo = 4;
  glm::vec4  bbox = {0.0f, 100.0f, 0.0f, 100.0f};

  VirtualArray elevation(shape, tile_shape, halo);
  ComputeMode  cm;
  elevation.fill(0.5f, cm);

  RockSeedingOptions options;
  options.seed = 42;
  options.bbox = bbox;
  options.distribution = RockDistribution(0u, 0.5f, 2.0f, 1.5f);

  size_t rock_count = 80;

  RockField field = hmap::va::place_rocks_procedural(rock_count,
                                                     elevation,
                                                     0.f,
                                                     1.f,
                                                     1.f,
                                                     {},
                                                     {},
                                                     0.f,
                                                     1.f,
                                                     1.f,
                                                     options,
                                                     cm);

  EXPECT_GT(field.size(), 0u);
  EXPECT_NEAR(static_cast<double>(field.size()),
              static_cast<double>(rock_count),
              30.0);

  for (const auto &rock : field)
  {
    EXPECT_GE(rock.position.x, 0.0f);
    EXPECT_LE(rock.position.x, 100.0f);
    EXPECT_GE(rock.position.y, 0.0f);
    EXPECT_LE(rock.position.y, 100.0f);
    EXPECT_NEAR(rock.position.z, 0.5f, 1e-4f);
    EXPECT_GE(rock.radius, 0.5f);
    EXPECT_LE(rock.radius, 2.0f);
  }
}

TEST(RockProceduralTest, VirtualArrayExclusionMask)
{
  glm::ivec2 shape = {64, 64};
  glm::ivec2 tile_shape = {32, 32};
  int        halo = 4;
  glm::vec4  bbox = {0.0f, 1.0f, 0.0f, 1.0f};

  VirtualArray elevation(shape, tile_shape, halo);
  VirtualArray exclusion(shape, tile_shape, halo);
  ComputeMode  cm;
  elevation.fill(0.3f, cm);
  exclusion.fill(0.0f, cm);

  // Exclude top half (y >= 0.5)
  Array excl_array(shape, 0.0f);
  for (int j = 32; j < 64; ++j)
  {
    for (int i = 0; i < 64; ++i)
    {
      excl_array(i, j) = 1.0f;
    }
  }
  exclusion.from_array(excl_array, cm);

  RockSeedingOptions options;
  options.seed = 123;
  options.bbox = bbox;
  options.exclusion_threshold = 0.5f;

  RockField field = hmap::va::place_rocks_procedural(60,
                                                     elevation,
                                                     0.f,
                                                     1.f,
                                                     1.f,
                                                     {},
                                                     exclusion,
                                                     0.f,
                                                     1.f,
                                                     1.f,
                                                     options,
                                                     cm);

  EXPECT_GT(field.size(), 0u);

  for (const auto &rock : field)
  {
    EXPECT_LT(rock.position.y, 0.55f);
  }
}

TEST(RockProceduralTest, VirtualArrayEmptyAndZeroDensity)
{
  VirtualArray empty_va;
  ComputeMode  cm;

  RockField f1 = hmap::va::place_rocks_procedural(50, empty_va);
  EXPECT_EQ(f1.size(), 0u);

  glm::ivec2   shape = {32, 32};
  glm::ivec2   tile_shape = {16, 16};
  VirtualArray elevation(shape, tile_shape, 2);
  elevation.fill(0.5f, cm);

  VirtualArray zero_density(shape, tile_shape, 2);
  zero_density.fill(0.0f, cm);

  RockField f2 = hmap::va::place_rocks_procedural(50,
                                                  elevation,
                                                  0.f,
                                                  1.f,
                                                  1.f,
                                                  zero_density,
                                                  {},
                                                  0.f,
                                                  1.f,
                                                  1.f,
                                                  {},
                                                  cm);
  EXPECT_EQ(f2.size(), 0u);
}

TEST(RockDensityTest, TiledBoundingBox)
{
  glm::ivec2 shape = {32, 32};
  Array      hmap(shape, 0.5f);

  // Compute with two distinct tile bounding boxes
  Array d_tile1 = build_rock_density(hmap,
                                     0.05f,
                                     0.01f,
                                     1.0f,
                                     1.0f,
                                     1.0f,
                                     0.1f,
                                     1000.0f,
                                     {},
                                     1.0f,
                                     {},
                                     0.f,
                                     BBox2D({0.f, 0.5f, 0.f, 0.5f}));

  Array d_tile2 = build_rock_density(hmap,
                                     0.05f,
                                     0.01f,
                                     1.0f,
                                     1.0f,
                                     1.0f,
                                     0.1f,
                                     1000.0f,
                                     {},
                                     1.0f,
                                     {},
                                     0.f,
                                     BBox2D({0.5f, 1.f, 0.5f, 1.f}));

  EXPECT_EQ(d_tile1.shape.x, 32);
  EXPECT_EQ(d_tile2.shape.x, 32);
  EXPECT_FALSE(d_tile1.vector.empty());
  EXPECT_FALSE(d_tile2.vector.empty());
}
