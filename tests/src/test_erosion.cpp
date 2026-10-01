#include "highmap/dbg/assert.hpp"
#include "highmap/erosion.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/morphology.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/primitives.hpp"
#include "highmap/range.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(Erosion, IdentityWhenRadiusZero)
{
  Array input = Array({{1, 2}, {3, 4}});
  Array cpu = erosion(input, 0);
  Array gpu = gpu::erosion(input, 0);

  EXPECT_TRUE(assert_almost_equal(cpu, input));
  EXPECT_TRUE(assert_almost_equal(gpu, input));
}

TEST(Erosion, SingleMinimumPropagation)
{
  Array input = Array({{5, 5, 5}, {5, 1, 5}, {5, 5, 5}});

  // square kernel for GPU
  Array cpu = erosion(input, 1);
  Array expected_cpu = Array({{1, 1, 1}, {1, 1, 1}, {1, 1, 1}});

  // disk kernel for GPU
  Array gpu = gpu::erosion(input, 1);
  Array expected_gpu = Array({{5, 1, 5}, {1, 1, 1}, {5, 1, 5}});

  EXPECT_TRUE(assert_almost_equal(cpu, expected_cpu));
  EXPECT_TRUE(assert_almost_equal(gpu, expected_gpu));
}

TEST(Erosion, EdgeHandling)
{
  Array input = Array({{3, 3, 3}, {3, 1, 3}, {3, 3, 3}});

  Array cpu = erosion(input, 1);
  Array gpu = gpu::erosion(input, 1);

  EXPECT_EQ(cpu(0, 0), 1);
  EXPECT_EQ(cpu(2, 2), 1);

  EXPECT_EQ(gpu(0, 0), 3);
  EXPECT_EQ(gpu(2, 2), 3);
}

TEST(Erosion, LargerRadius)
{
  Array input = Array({{9, 8, 7, 6}, {8, 5, 4, 7}, {7, 4, 3, 8}, {6, 7, 8, 9}});
  Array expected = Array(
      {{3, 3, 3, 3}, {3, 3, 3, 3}, {3, 3, 3, 3}, {3, 3, 3, 3}});

  Array cpu = erosion(input, 3);
  Array gpu = gpu::erosion(input, 3);

  EXPECT_TRUE(assert_almost_equal(cpu, expected));
  EXPECT_TRUE(assert_almost_equal(gpu, expected));
}

TEST(Erosion, MonotonicDecrease)
{
  Array input = Array({{5, 6}, {7, 8}});

  Array cpu = erosion(input, 1);
  Array gpu = gpu::erosion(input, 1);

  for (int i = 0; i < input.shape.x; ++i)
    for (int j = 0; j < input.shape.y; ++j)
    {
      EXPECT_LE(cpu(i, j), input(i, j));
      EXPECT_LE(gpu(i, j), input(i, j));
    }
}

TEST(Erosion, FlatRegionUnchanged)
{
  Array input = Array({{2, 2, 2}, {2, 2, 2}, {2, 2, 2}});

  Array cpu = erosion(input, 1);
  Array gpu = gpu::erosion(input, 1);

  EXPECT_TRUE(assert_almost_equal(cpu, input));
  EXPECT_TRUE(assert_almost_equal(gpu, input));
}

TEST(Erosion, NonSquareArray)
{
  Array input = Array({{5, 4, 3, 2}});
  Array expected = Array({{4, 3, 2, 2}});

  Array cpu = erosion(input, 1);
  Array gpu = gpu::erosion(input, 1);

  EXPECT_TRUE(assert_almost_equal(cpu, expected));
  EXPECT_TRUE(assert_almost_equal(gpu, expected));
}

TEST(ThermalConserve, MassPreserved)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {64, 64};
  glm::vec2  kw = {4.f, 4.f};
  Array      z0 = noise_fbm(NoiseType::PERLIN, shape, kw, 42);

  Array z = z0;
  gpu::thermal_conserve(z, 0.1f / shape.x, 100);

  EXPECT_NEAR(z.sum(), z0.sum(), 1e-2f);
}

TEST(ThermalGPU, VariantsRunAndModify)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {64, 64};
  glm::vec2  kw = {4.f, 4.f};
  Array      z0 = noise_fbm(NoiseType::PERLIN, shape, kw, 42);
  float      talus = 0.5f / shape.x;
  Array      talus_map(shape, talus);

  // thermal
  {
    Array z = z0;
    Array dep(shape);
    gpu::thermal(z, talus_map, 5, nullptr, &dep);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_auto_bedrock
  {
    Array z = z0;
    Array dep(shape);
    gpu::thermal_auto_bedrock(z, talus_map, 5, &dep);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_flatten
  {
    Array z = z0;
    gpu::thermal_flatten(z, talus_map, 5);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_inflate
  {
    Array z = z0;
    gpu::thermal_inflate(z, talus_map, 5);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_olsen
  {
    Array z = z0;
    gpu::thermal_olsen(z, talus_map, 5);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_rib
  {
    Array z = z0;
    gpu::thermal_rib(z, 5);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_ridge
  {
    Array z = z0;
    Array dep(shape);
    gpu::thermal_ridge(z, talus_map, 5, &dep);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_schott
  {
    Array z = z0;
    Array dep(shape);
    gpu::thermal_schott(z, talus_map, 5, 0.2f, &dep);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }

  // thermal_scree
  {
    Array z = z0;
    Array zmax(shape, 1.f);
    Array dep(shape);
    gpu::thermal_scree(z, talus_map, zmax, 5, &dep);
    EXPECT_FALSE(assert_almost_equal(z, z0));
  }
}

TEST(ConvErosion, EmptyArray)
{
  Array empty;
  gpu::conv_erosion(empty, 42, 5, 10, 1, 2);
  EXPECT_TRUE(empty.vector.empty());
}

TEST(ConvErosion, BasicExecution)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {64, 64};
  glm::vec2  kw = {4.f, 4.f};
  Array      z0 = noise_fbm(NoiseType::PERLIN, shape, kw, 42);
  Array      z = z0;

  gpu::conv_erosion(z, 42, 5, 200, 1, 4, 1.f, 0.05f);

  EXPECT_EQ(z.shape, shape);
  EXPECT_FALSE(assert_almost_equal(z, z0));
}

TEST(DepressionFilling, BoundaryConditions)
{
  // A 5x5 bowl with a pit in the center and high borders
  Array z({{5.f, 5.f, 5.f, 5.f, 5.f},
           {5.f, 1.f, 1.f, 1.f, 5.f},
           {5.f, 1.f, 0.f, 1.f, 5.f},
           {5.f, 1.f, 1.f, 1.f, 5.f},
           {5.f, 5.f, 5.f, 5.f, 5.f}});

  // Create an outlet on the right boundary (i = 4, j = 2) and lower the lip at
  // (3, 2)
  z(3, 2) = 0.2f;
  z(4, 2) = 0.1f;

  // Case 1: outflow allowed on all boundaries
  {
    Array z_test = z;
    depression_filling(z_test, 100, 1e-4f, true, true, true, true);
    // Pit (2, 2) fills up to the outlet level (approx 0.2f)
    EXPECT_NEAR(z_test(2, 2), 0.2f, 0.01f);
  }

  // Case 2: outflow closed on the right boundary (where the outlet is)
  {
    Array z_test = z;
    depression_filling(z_test, 100, 1e-4f, true, false, true, true);
    // Since right boundary is closed and all other boundaries are 5.0f, the pit
    // fills up to 5.0f
    EXPECT_NEAR(z_test(2, 2), 5.0f, 0.01f);
  }
}

TEST(DepressionFillingPriorityFlood, BoundaryConditions)
{
  // A 5x5 bowl with a pit in the center and high borders
  Array z({{5.f, 5.f, 5.f, 5.f, 5.f},
           {5.f, 1.f, 1.f, 1.f, 5.f},
           {5.f, 1.f, 0.f, 1.f, 5.f},
           {5.f, 1.f, 1.f, 1.f, 5.f},
           {5.f, 5.f, 5.f, 5.f, 5.f}});

  // Create an outlet on the right boundary (i = 4, j = 2) and lower the lip at
  // (3, 2)
  z(3, 2) = 0.2f;
  z(4, 2) = 0.1f;

  // Case 1: outflow allowed on all boundaries
  {
    Array z_test = z;
    depression_filling_priority_flood(z_test, false, true, true, true, true);
    // Pit (2, 2) fills up to the outlet level (approx 0.2f)
    EXPECT_NEAR(z_test(2, 2), 0.2f, 0.01f);
  }

  // Case 2: outflow closed on the right boundary (where the outlet is)
  {
    Array z_test = z;
    depression_filling_priority_flood(z_test, false, true, false, true, true);
    // Since right boundary is closed and all other boundaries are 5.0f, the pit
    // fills up to 5.0f
    EXPECT_NEAR(z_test(2, 2), 5.0f, 0.01f);
  }
}

TEST(DepressionFillingPriorityFlood_VirtualArray, SingleTileMatchesArray)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{64, 64};
  const int        halo = 0;

  glm::vec2 kw = {4.f, 4.f};
  Array     input = noise_fbm(NoiseType::PERLIN, shape, kw, 42);

  VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_SEQUENTIAL};
  va.from_array(input, cm);

  VirtualArray va_fill_map(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
  VirtualArray va_out = va::depression_filling_priority_flood(va,
                                                              false,
                                                              &va_fill_map,
                                                              cm);

  Array expected = input;
  depression_filling_priority_flood(expected, false);

  Array result = va_out.to_array(cm);
  EXPECT_TRUE(assert_almost_equal(result, expected, 1e-5f));

  Array result_fill = va_fill_map.to_array(cm);
  Array expected_fill = expected - input;
  EXPECT_TRUE(assert_almost_equal(result_fill, expected_fill, 1e-5f));
}

TEST(DepressionFillingPriorityFlood_VirtualArray, MultiTileWithFillMap)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  glm::vec2 kw = {4.f, 4.f};
  Array     input = noise_fbm(NoiseType::PERLIN, shape, kw, 42);

  VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_DISTRIBUTED};
  va.from_array(input, cm);

  VirtualArray va_fill_map(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
  VirtualArray va_out = va::depression_filling_priority_flood(va,
                                                              false,
                                                              &va_fill_map,
                                                              cm);

  Array result = va_out.to_array(cm);
  Array result_fill = va_fill_map.to_array(cm);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_GE(result(i, j), input(i, j) - 1e-5f);
      EXPECT_GE(result_fill(i, j), -1e-5f);
      EXPECT_NEAR(result_fill(i, j), result(i, j) - input(i, j), 1e-5f);
    }
}

TEST(FloodingLakeSystem_VirtualArray, SingleTileMatchesArray)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{64, 64};
  const int        halo = 0;

  glm::vec2 kw = {4.f, 4.f};
  Array     input = noise_fbm(NoiseType::PERLIN, shape, kw, 42);

  VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_SEQUENTIAL};
  va.from_array(input, cm);

  VirtualArray va_water = va::flooding_lake_system(va, 0.f, cm);

  Array expected = flooding_lake_system(input, 0.f);
  Array result = va_water.to_array(cm);

  EXPECT_TRUE(assert_almost_equal(result, expected, 1e-5f));
}

TEST(FloodingLakeSystem_VirtualArray, MultiTileWithSurfaceThreshold)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  glm::vec2 kw = {4.f, 4.f};
  Array     input = noise_fbm(NoiseType::PERLIN, shape, kw, 42);

  VirtualArray va(shape, bbox, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_DISTRIBUTED};
  va.from_array(input, cm);

  VirtualArray va_water = va::flooding_lake_system(va, 10.f, cm);
  Array        result = va_water.to_array(cm);

  EXPECT_GE(result.min(), 0.f);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      EXPECT_GE(result(i, j), 0.f);
}

TEST(HydraulicMusgrave, BasicExecution)
{
  glm::ivec2 shape = {64, 64};
  glm::vec2  kw = {4.f, 4.f};
  Array      z0 = noise_fbm(NoiseType::PERLIN, shape, kw, 42);
  Array      z = z0;

  hydraulic_musgrave(z, 20);

  EXPECT_EQ(z.shape, shape);
  EXPECT_FALSE(assert_almost_equal(z, z0));

  for (int j = 0; j < z.shape.y; j++)
    for (int i = 0; i < z.shape.x; i++)
    {
      EXPECT_FALSE(std::isnan(z(i, j)));
      EXPECT_FALSE(std::isinf(z(i, j)));
    }
}

TEST(HydraulicMusgrave, FlatRegionUnchanged)
{
  Array z = constant(glm::ivec2(16, 16), 5.f);
  Array z0 = z;

  hydraulic_musgrave(z, 10);

  EXPECT_TRUE(assert_almost_equal(z, z0));
}

TEST(HydraulicMusgraveGPU, BasicExecution)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {64, 64};
  glm::vec2  kw = {4.f, 4.f};
  Array      z0 = noise_fbm(NoiseType::PERLIN, shape, kw, 42);
  Array      z = z0;

  gpu::hydraulic_musgrave(z, 20);

  EXPECT_EQ(z.shape, shape);
  EXPECT_FALSE(assert_almost_equal(z, z0));

  for (int j = 0; j < z.shape.y; j++)
    for (int i = 0; i < z.shape.x; i++)
    {
      EXPECT_FALSE(std::isnan(z(i, j)));
      EXPECT_FALSE(std::isinf(z(i, j)));
    }
}

TEST(HydraulicMusgraveGPU, FlatRegionUnchanged)
{
  hmap::gpu::init_opencl();

  Array z = constant(glm::ivec2(16, 16), 5.f);
  Array z0 = z;

  gpu::hydraulic_musgrave(z, 10);

  EXPECT_TRUE(assert_almost_equal(z, z0));
}

TEST(HydraulicFastScape, BasicExecution)
{
  glm::ivec2 shape = {64, 64};
  glm::vec2  kw = {4.f, 4.f};
  Array      z0 = noise_fbm(NoiseType::PERLIN, shape, kw, 42);
  remap(z0, 0.f, 1.f);
  Array z = z0;

  Array erosion_map(shape, 0.f);
  Array flow_map(shape, 0.f);

  hydraulic_fastscape(z,
                      15,
                      1e-2f,
                      1.f,
                      0.5f,
                      1.f,
                      1e-3f,
                      0.f,
                      true,
                      1.f,
                      1e-3f,
                      nullptr,
                      nullptr,
                      &erosion_map,
                      &flow_map);

  EXPECT_EQ(z.shape, shape);
  EXPECT_FALSE(assert_almost_equal(z, z0));

  for (int j = 0; j < z.shape.y; j++)
    for (int i = 0; i < z.shape.x; i++)
    {
      EXPECT_FALSE(std::isnan(z(i, j)));
      EXPECT_FALSE(std::isinf(z(i, j)));
      EXPECT_GE(flow_map(i, j), 1.f);
      EXPECT_GE(erosion_map(i, j), 0.f);
    }
}

TEST(HydraulicFastScape, FlatRegionUnchanged)
{
  Array z = constant(glm::ivec2(16, 16), 5.f);
  Array z0 = z;

  hydraulic_fastscape(z, 10, 1e-2f, 1.f, 0.5f, 1.f, 0.f, 0.f);

  EXPECT_TRUE(assert_almost_equal(z, z0));
}

TEST(HydraulicFastScape, NonlinearExponent)
{
  glm::ivec2 shape = {32, 32};
  glm::vec2  kw = {4.f, 4.f};
  Array      z0 = noise_fbm(NoiseType::PERLIN, shape, kw, 123);
  remap(z0, 0.f, 1.f);
  Array z = z0;

  hydraulic_fastscape(z, 5, 1e-2f, 1.f, 0.5f, 2.f, 1e-3f, 0.f);

  EXPECT_EQ(z.shape, shape);
  for (int j = 0; j < z.shape.y; j++)
    for (int i = 0; i < z.shape.x; i++)
    {
      EXPECT_FALSE(std::isnan(z(i, j)));
      EXPECT_FALSE(std::isinf(z(i, j)));
    }
}

TEST(HydraulicFastScape, VirtualArrayExecution)
{
  glm::ivec2 shape = {64, 64};
  glm::vec2  kw = {4.f, 4.f};
  Array      z0 = noise_fbm(NoiseType::PERLIN, shape, kw, 42);
  remap(z0, 0.f, 1.f);

  ComputeMode  cm = {ForEachMode::VA_SINGLE_ARRAY};
  VirtualArray va(shape, {32, 32}, 4);
  va.from_array(z0, cm);

  VirtualArray va_eroded = va::hydraulic_fastscape(cm, va, 10);
  Array        z_res = va_eroded.to_array(cm);

  EXPECT_EQ(z_res.shape, shape);
  EXPECT_FALSE(assert_almost_equal(z_res, z0));
}

TEST(HydraulicFastScape, MultiscaleExecution)
{
  glm::ivec2 shape = {64, 64};
  glm::vec2  kw = {4.f, 4.f};
  Array      z0 = noise_fbm(NoiseType::PERLIN, shape, kw, 42);
  remap(z0, 0.f, 1.f);
  Array z = z0;

  Array erosion_map(shape, 0.f);
  Array flow_map(shape, 0.f);

  hydraulic_fastscape_multiscale(z,
                                 {8, 4, 2},
                                 1e-2f,
                                 1.f,
                                 0.5f,
                                 1.f,
                                 1e-3f,
                                 0.f,
                                 true,
                                 1.f,
                                 1e-3f,
                                 nullptr,
                                 nullptr,
                                 &erosion_map,
                                 &flow_map,
                                 0.8f);

  EXPECT_EQ(z.shape, shape);
  EXPECT_FALSE(assert_almost_equal(z, z0));

  for (int j = 0; j < z.shape.y; j++)
    for (int i = 0; i < z.shape.x; i++)
    {
      EXPECT_FALSE(std::isnan(z(i, j)));
      EXPECT_FALSE(std::isinf(z(i, j)));
      EXPECT_GE(flow_map(i, j), 1.f);
      EXPECT_GE(erosion_map(i, j), 0.f);
    }
}
