#include "highmap/filters.hpp"
#include "highmap/internal/assert.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(SmoothCPulse, ConstantFieldPreserved)
{
  Array input = Array({{2, 2, 2}, {2, 2, 2}, {2, 2, 2}});

  Array cpu = input;
  smooth_cpulse(cpu, 1);

  EXPECT_TRUE(assert_almost_equal(cpu, input));
}

TEST(SmoothCPulse, IdentityWhenRadiusZero)
{
  Array input = Array({{1, 2, 3}, {4, 5, 6}});

  Array cpu = input;
  smooth_cpulse(cpu, 0);

  EXPECT_TRUE(assert_almost_equal(cpu, input));
}

TEST(SmoothCPulse, SmoothingReducesContrast)
{
  Array input = Array({{0, 0, 0}, {0, 10, 0}, {0, 0, 0}});

  Array cpu = input;
  smooth_cpulse(cpu, 1);

  // center should decrease due to diffusion
  EXPECT_LT(cpu(1, 1), 10.f);
}

TEST(SmoothCPulse, SymmetryPreserved)
{
  Array input = Array({{1, 2, 3, 2, 1},
                       {2, 3, 4, 3, 2},
                       {3, 4, 5, 4, 3},
                       {2, 3, 4, 3, 2},
                       {1, 2, 3, 2, 1}});

  Array cpu = input;
  smooth_cpulse(cpu, 1);

  // symmetry along both axes should remain approximately valid
  for (int i = 0; i < input.shape.x; ++i)
    for (int j = 0; j < input.shape.y; ++j)
    {
      EXPECT_NEAR(cpu(i, j), cpu(input.shape.x - 1 - i, j), 1e-4);

      EXPECT_NEAR(cpu(i, j), cpu(i, input.shape.y - 1 - j), 1e-4);
    }
}

TEST(SmoothCPulse, NoNewExtremaCreated)
{
  Array input = Array({{0, 5, 0}, {5, 10, 5}, {0, 5, 0}});

  float min_before = input.min();
  float max_before = input.max();

  Array cpu = input;
  smooth_cpulse(cpu, 1);

  float min_after = cpu.min();
  float max_after = cpu.max();

  EXPECT_GE(min_after, min_before - 1e-5f);
  EXPECT_LE(max_after, max_before + 1e-5f);
}

TEST(SmoothCPulse_CPU_GPU, RandomBinaryFields)
{
  std::mt19937                          rng(42);
  std::uniform_real_distribution<float> dist(0.f, 1.f);

  const int nx = 64;
  const int ny = 64;

  for (int test = 0; test < 20; ++test)
  {
    Array input(glm::ivec2(nx, ny));

    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        input(i, j) = dist(rng);

    int ir = int(nx * dist(rng));

    Array cpu = input;
    Array gpu = input;

    smooth_cpulse(cpu, ir);
    gpu::smooth_cpulse(gpu, ir);

    // same arrays are expected
    EXPECT_TRUE(assert_almost_equal(cpu, gpu, 1e-6f));
  }
}

TEST(SmoothCPulse_VirtualArray, SingleTileMatchesGPU)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{64, 64};
  const int        halo = 0;

  Array input(shape);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      input(i, j) = float(i + j * shape.x);

  VirtualArray va(shape, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm_seq{.mode = ForEachMode::VA_SEQUENTIAL};
  va.from_array(input, cm_seq);

  const int ir = 5;
  va::smooth_cpulse(va, ir, nullptr, cm_seq);

  Array expected = input;
  gpu::smooth_cpulse(expected, ir);

  Array result = va.to_array(cm_seq);
  EXPECT_TRUE(assert_almost_equal(result, expected, 1e-5f));
}

TEST(SmoothCPulse_VirtualArray, MaskedFilter)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 8;

  Array input(shape);
  Array mask(shape);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      input(i, j) = float(i * i + j);
      mask(i, j) = (i > 32) ? 1.f : 0.f;
    }

  VirtualArray va_input(shape, tile_shape, halo, StorageMode::VA_RAM);
  VirtualArray va_mask(shape, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm_dist{.mode = ForEachMode::VA_DISTRIBUTED};

  va_input.from_array(input, cm_dist);
  va_mask.from_array(mask, cm_dist);

  const int ir = 4;
  va::smooth_cpulse(va_input, ir, &va_mask, cm_dist);

  Array expected = input;
  gpu::smooth_cpulse(expected, ir, &mask);

  // on unmasked region (left half), values should remain identical to input
  Array result = va_input.to_array(cm_dist);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < 30; ++i)
    {
      EXPECT_NEAR(result(i, j), input(i, j), 1e-5f);
    }
}

TEST(SmoothCPulse_VirtualArray, ReturnsFilteredVirtualArray)
{
  const glm::ivec2 shape{64, 64};
  const glm::vec4  bbox{0.f, 1.f, 0.f, 1.f};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 4;

  Array input(shape);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      input(i, j) = float(i + j * shape.x);

  VirtualArray va_in(shape, tile_shape, halo, StorageMode::VA_RAM);
  ComputeMode  cm{.mode = ForEachMode::VA_DISTRIBUTED};
  va_in.from_array(input, cm);

  const int           ir = 3;
  const VirtualArray &va_in_const = va_in;
  VirtualArray        va_out = va::smooth_cpulse(va_in_const, ir, nullptr, cm);

  Array expected = input;
  gpu::smooth_cpulse(expected, ir);

  Array result = va_out.to_array(cm);
  EXPECT_TRUE(assert_almost_equal(result, expected, 1e-5f));
}
