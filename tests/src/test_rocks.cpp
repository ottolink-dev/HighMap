#include <cmath>
#include <fstream>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/rocks.hpp"

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
  field.pack_interstitial_rocks(20, {}, 0.01f, 0.02f, 1u, 42, bbox);

  EXPECT_GT(field.size(), original_count);
  EXPECT_EQ(field[0].radius, 0.1f);
  EXPECT_EQ(field[1].radius, 0.1f);

  // 2. Pack with density map
  Array  density({32, 32}, 1.0f);
  size_t count_before_density = field.size();
  field.pack_interstitial_rocks(20, density, 0.005f, 0.01f, 2u, 101, bbox);
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
  opts.max_steps = 1000;
  opts.sub_steps = 10;
  opts.soil_friction = 0.55f;

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
  opts.sub_steps = 10;
  opts.soil_friction = 0.1f;

  field.simulate_physics(elev, nullptr, opts);

  // Rock that rolled out of bounds must be discarded
  EXPECT_EQ(field.size(), 0u);
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
  opts.sub_steps = 10;
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
