#include "highmap/array.hpp"
#include "highmap/convolve.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/filters.hpp"
#include "highmap/kernels.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/primitives.hpp"

#include <gtest/gtest.h>

using namespace hmap;

// ------------------------------------------------------------
// Scalability Tests: Full 2D Convolution and Window Operators
// ------------------------------------------------------------

TEST(ConvolutionScaling, CPU_Convolve2D_DirectScaling)
{
  // Test direct 2D convolution across combinations of resolution and kernel
  // size
  const std::vector<int> sizes = {32, 64};
  const std::vector<int> radii = {2, 4, 8};

  for (int n : sizes)
  {
    Array input = white(glm::vec2(n, n), 0.f, 1.f, 42);

    for (int r : radii)
    {
      int   ksize = 2 * r + 1;
      Array kernel = cubic_pulse(glm::ivec2(ksize, ksize));

      Array out = convolve2d(input, kernel);

      EXPECT_EQ(out.shape.x, n);
      EXPECT_EQ(out.shape.y, n);
      EXPECT_FALSE(out.vector.empty());
    }
  }
}

TEST(ConvolutionScaling, CPU_Convolve2D_SvdScaling)
{
  // SVD convolution should scale gracefully for larger resolutions and kernel
  // sizes
  const std::vector<int> sizes = {64, 128};
  const std::vector<int> radii = {4, 16, 28};

  for (int n : sizes)
  {
    Array input = white(glm::vec2(n, n), 0.f, 1.f, 42);

    for (int r : radii)
    {
      int   ksize = 2 * r + 1;
      Array kernel = cubic_pulse(glm::ivec2(ksize, ksize));

      Array out = convolve2d_svd(input, kernel, 3);

      EXPECT_EQ(out.shape.x, n);
      EXPECT_EQ(out.shape.y, n);
      EXPECT_FALSE(out.vector.empty());
    }
  }
}

TEST(ConvolutionScaling, GPU_BilateralFilterScaling)
{
  // OpenCL 2D brute-force window: test scalability with larger kernel sizes
  // relative to grid
  const std::vector<int> sizes = {64, 128};
  const std::vector<int> radii = {2, 8, 16};

  for (int n : sizes)
  {
    Array input = white(glm::vec2(n, n), 0.f, 1.f, 42);

    for (int r : radii)
    {
      Array out = gpu::bilateral_filter(input, r, 1.f);

      EXPECT_EQ(out.shape.x, n);
      EXPECT_EQ(out.shape.y, n);
      EXPECT_FALSE(out.vector.empty());
    }
  }
}

TEST(ConvolutionScaling, GPU_ExpandScaling)
{
  // OpenCL 2D morphology/dilation with 2D kernel
  const std::vector<int> sizes = {64, 128};
  const std::vector<int> radii = {2, 8, 16};

  for (int n : sizes)
  {
    Array input = white(glm::vec2(n, n), 0.f, 1.f, 42);

    for (int r : radii)
    {
      Array out = input;
      gpu::expand(out, r, 1);

      EXPECT_EQ(out.shape.x, n);
      EXPECT_EQ(out.shape.y, n);
      EXPECT_GE(out.max(), input.max());
    }
  }
}

TEST(ConvolutionScaling, GPU_LocalMetricsScaling)
{
  // OpenCL circular window operators: local_max, ruggedness,
  // topographic_position_index
  const std::vector<int> sizes = {64, 128};
  const std::vector<int> radii = {2, 8, 16};

  for (int n : sizes)
  {
    Array input = white(glm::vec2(n, n), 0.f, 1.f, 42);

    for (int r : radii)
    {
      Array out_max = gpu::local_max(input, r);
      Array out_rug = gpu::ruggedness(input, r);
      Array out_tpi = gpu::topographic_position_index(input, r);

      EXPECT_EQ(out_max.shape.x, n);
      EXPECT_EQ(out_rug.shape.x, n);
      EXPECT_EQ(out_tpi.shape.x, n);
    }
  }
}

TEST(ConvolutionScaling, GPU_SparseMaxConvolutionScaling)
{
  // OpenCL scatter/atomic max-convolution
  const std::vector<int> sizes = {64, 128};
  const std::vector<int> radii = {4, 16};

  for (int n : sizes)
  {
    Array input = white_sparse(glm::ivec2(n, n), 0.f, 1.f, 0.05f, 42);

    for (int r : radii)
    {
      int   ksize = 2 * r + 1;
      Array kernel = cubic_pulse(glm::ivec2(ksize, ksize));

      Array out = gpu::sparse_max_convolution(input, kernel, false);

      EXPECT_EQ(out.shape.x, n);
      EXPECT_EQ(out.shape.y, n);
      EXPECT_FALSE(out.vector.empty());
    }
  }
}
