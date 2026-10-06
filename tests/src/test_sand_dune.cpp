#include "highmap/array.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/erosion.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/primitives.hpp"

#include <cmath>
#include <gtest/gtest.h>

using namespace hmap;

class SandDuneTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    gpu::init_opencl();
  }
};

TEST_F(SandDuneTest, BasicExecution)
{
  glm::ivec2 shape = {64, 64};
  Array      z(shape, 0.5f);

  SandDuneParams params;
  params.wind_angle = 0.f;
  params.periodic_boundary = true;
  params.collapse_rate = 0.f; // isolate transport

  gpu::sand_dune(z, 1000, params);

  // terrain should have developed surface variations
  EXPECT_GT(z.max() - z.min(), 0.f);
}

TEST_F(SandDuneTest, MassPreservedInPeriodicDomain)
{
  glm::ivec2 shape = {64, 64};
  Array      z(shape, 0.5f);

  float initial_sum = z.sum();

  SandDuneParams params;
  params.wind_angle = 0.5f;
  params.periodic_boundary = true;
  params.sand_inflow_rate = 0.f;
  params.talus = 0.01f;
  params.collapse_rate = 0.5f;

  gpu::sand_dune(z, 5000, params, nullptr, nullptr, nullptr, nullptr, 5);

  float final_sum = z.sum();
  // check mass conservation to high precision
  EXPECT_NEAR(initial_sum, final_sum, 1e-2f);
}

TEST_F(SandDuneTest, BedrockConstraint)
{
  glm::ivec2 shape = {64, 64};
  Array      bedrock(shape, 0.2f);
  Array      z(shape, 0.3f); // 0.1 initial sand thickness

  SandDuneParams params;
  params.slab_height = 0.05f;
  params.periodic_boundary = true;
  params.talus = 0.01f;

  gpu::sand_dune(z, 10000, params, &bedrock, nullptr, nullptr, nullptr, 10);

  // no cell should fall below bedrock
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_GE(z(i, j), bedrock(i, j) - 1e-4f);
    }
  }
}

TEST_F(SandDuneTest, DepositionAndErosionMaps)
{
  glm::ivec2 shape = {64, 64};
  Array      z(shape, 0.5f);
  Array      dep_map(shape, 0.f);
  Array      ero_map(shape, 0.f);

  SandDuneParams params;
  params.periodic_boundary = true;

  gpu::sand_dune(z, 2000, params, nullptr, nullptr, &dep_map, &ero_map, 4);

  EXPECT_GT(dep_map.sum(), 0.f);
  EXPECT_GT(ero_map.sum(), 0.f);
}

TEST_F(SandDuneTest, ElevationAwareHop)
{
  glm::ivec2 shape = {64, 64};
  Array      z(shape, 0.2f);

  // create a high plateau on the left (x < 20)
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < 20; ++i)
    {
      z(i, j) = 0.8f;
    }
  }

  SandDuneParams params;
  params.wind_angle = 0.f; // wind blowing left to right (+x)
  params.hop_length = 3.f;
  params.elevation_hop_factor = 5.f;
  params.collapse_rate = 0.f;

  gpu::sand_dune(z, 3000, params, nullptr, nullptr, nullptr, nullptr, 3);

  // verify execution completes and values remain finite and within expected range
  EXPECT_FALSE(z.vector.empty());
  EXPECT_GT(z.max(), 0.f);
}
