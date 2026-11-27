#include "highmap/filters.hpp"
#include "highmap/internal/assert.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <cmath>

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
  const Array input = Array({{1, 2, 3}, {4, 5, 6}});
  const Array mask = Array({{0.f, 0.5f, 1.f}, {1.f, 0.5f, 0.f}});
  for (const Array *p_mask : {static_cast<const Array *>(nullptr), &mask})
  {
    Array cpu = input;
    smooth_cpulse(cpu, 0, p_mask);
    EXPECT_EQ(cpu.vector, input.vector);
    for (float value : cpu.vector) EXPECT_TRUE(std::isfinite(value));

    // The identity path does not require an OpenCL device or context.
    Array gpu = input;
    gpu::smooth_cpulse(gpu, 0, p_mask);
    EXPECT_EQ(gpu.vector, input.vector);
    for (float value : gpu.vector) EXPECT_TRUE(std::isfinite(value));
  }
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
