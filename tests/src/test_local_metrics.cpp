#include "highmap/dbg/assert.hpp"
#include "highmap/local_metrics.hpp"

#include <gtest/gtest.h>

using namespace hmap;

// ------------------------------------------------------------
// LOCAL MAX
// ------------------------------------------------------------

TEST(LocalMetrics, LocalMax_BasicPeakPropagation)
{
  Array input = Array({{1, 2, 1}, {1, 5, 1}, {1, 1, 1}});

  Array cpu = local_max(input, 1);
  Array gpu = gpu::local_max(input, 1);
  Array gpu_sq = gpu::local_max_square(input, 1);

  EXPECT_EQ(cpu(0, 0), 5);
  EXPECT_EQ(cpu(2, 2), 5);

  // disk kernel
  EXPECT_EQ(gpu(0, 0), 2);
  EXPECT_EQ(gpu(0, 1), 5);
  EXPECT_EQ(gpu(1, 0), 5);

  // square separable GPU kernel matches CPU square kernel
  EXPECT_TRUE(assert_almost_equal(cpu, gpu_sq));
}

TEST(LocalMetrics, LocalMax_Convergence)
{
  Array input = Array({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});

  Array a_cpu = local_max(local_max(input, 1), 1);
  Array b_cpu = local_max(a_cpu, 1);

  Array a_gpu_sq = gpu::local_max_square(gpu::local_max_square(input, 1), 1);
  Array b_gpu_sq = gpu::local_max_square(a_gpu_sq, 1);

  // more iterations because of the disk kernel
  Array a_gpu = input;
  for (int it = 0; it < 4; ++it)
    a_gpu = gpu::local_max(a_gpu, 1);
  Array b_gpu = gpu::local_max(a_gpu, 1);

  EXPECT_TRUE(assert_almost_equal(a_cpu, b_cpu));
  EXPECT_TRUE(assert_almost_equal(a_cpu, a_gpu_sq));
  EXPECT_TRUE(assert_almost_equal(b_cpu, b_gpu_sq));
  EXPECT_TRUE(assert_almost_equal(a_gpu, b_gpu));
}

TEST(LocalMetrics, LocalMax_Monotonicity)
{
  Array input = Array({{1, 3, 2}, {4, 1, 0}, {2, 5, 1}});

  Array cpu = local_max(input, 1);
  Array gpu = gpu::local_max(input, 1);
  Array gpu_sq = gpu::local_max_square(input, 1);

  for (int i = 0; i < input.shape.x; ++i)
    for (int j = 0; j < input.shape.y; ++j)
    {
      EXPECT_GE(cpu(i, j), input(i, j));
      EXPECT_GE(gpu(i, j), input(i, j));
      EXPECT_GE(gpu_sq(i, j), input(i, j));
      EXPECT_FLOAT_EQ(cpu(i, j), gpu_sq(i, j));
    }
}

// ------------------------------------------------------------
// LOCAL MIN
// ------------------------------------------------------------

TEST(LocalMetrics, LocalMin_DualityWithMax)
{
  Array input = Array({{3, 2, 1}, {4, 0, 5}, {6, 7, 8}});

  Array cpu = local_min(input, 1);
  Array gpu = gpu::local_min(input, 1);
  Array gpu_sq = gpu::local_min_square(input, 1);

  Array cpu_ref = -local_max(-input, 1);
  Array gpu_ref = -gpu::local_max(-input, 1);
  Array gpu_sq_ref = -gpu::local_max_square(-input, 1);

  EXPECT_TRUE(assert_almost_equal(cpu, cpu_ref));
  EXPECT_TRUE(assert_almost_equal(cpu, gpu_sq));
  EXPECT_TRUE(assert_almost_equal(gpu, gpu_ref));
  EXPECT_TRUE(assert_almost_equal(gpu_sq, gpu_sq_ref));
}

TEST(LocalMetrics, LocalMin_Monotonicity)
{
  Array input = Array({{5, 6, 7}, {2, 1, 3}, {4, 8, 9}});

  Array cpu = local_min(input, 1);
  Array gpu = gpu::local_min(input, 1);
  Array gpu_sq = gpu::local_min_square(input, 1);

  for (int i = 0; i < input.shape.x; ++i)
    for (int j = 0; j < input.shape.y; ++j)
    {
      EXPECT_LE(cpu(i, j), input(i, j));
      EXPECT_LE(gpu(i, j), input(i, j));
      EXPECT_LE(gpu_sq(i, j), input(i, j));
      EXPECT_FLOAT_EQ(cpu(i, j), gpu_sq(i, j));
    }
}

// ------------------------------------------------------------
// LOCAL MEAN
// ------------------------------------------------------------

TEST(LocalMetrics, LocalMean_SmoothingEffect)
{
  Array input = Array({{0, 0, 0}, {0, 10, 0}, {0, 0, 0}});

  Array cpu = local_mean(input, 1);
  Array gpu = gpu::local_mean(input, 1);

  EXPECT_LT(cpu(1, 1), 10);
  EXPECT_LT(gpu(1, 1), 10);
}

TEST(LocalMetrics, LocalMean_PreservationOfConstantField)
{
  Array input = Array({{2, 2, 2}, {2, 2, 2}, {2, 2, 2}});

  Array cpu = local_mean(input, 1);
  Array gpu = gpu::local_mean(input, 1);

  EXPECT_TRUE(assert_almost_equal(cpu, input));
  EXPECT_TRUE(assert_almost_equal(gpu, input));
}

TEST(LocalMetrics, LocalMean_Linearity)
{
  Array a = Array({{1, 1, 1}, {1, 1, 1}});
  Array b = Array({{2, 2, 2}, {2, 2, 2}});

  Array mean_a_cpu = local_mean(a, 1);
  Array mean_b_cpu = local_mean(b, 1);

  Array mean_a_gpu = gpu::local_mean(a, 1);
  Array mean_b_gpu = gpu::local_mean(b, 1);

  for (int i = 0; i < a.shape.x; ++i)
    for (int j = 0; j < a.shape.y; ++j)
    {
      EXPECT_NEAR(mean_b_cpu(i, j), 2.f * mean_a_cpu(i, j), 1e-5);
      EXPECT_NEAR(mean_b_gpu(i, j), 2.f * mean_a_gpu(i, j), 1e-5);
    }
}
