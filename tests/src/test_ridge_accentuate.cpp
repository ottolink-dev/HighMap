#include <random>

#include "highmap/dbg/assert.hpp"
#include "highmap/filters.hpp"
#include "highmap/primitives.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(RidgeAccentuate, ConstantFieldPreserved)
{
  Array input = Array({{2, 2, 2}, {2, 2, 2}, {2, 2, 2}});

  Array cpu = input;
  ridge_accentuate(cpu, 0.1f, 1);

  EXPECT_TRUE(assert_almost_equal(cpu, input));
}

TEST(RidgeAccentuate, IdentityWhenStrengthZero)
{
  Array input = Array({{1, 2, 3}, {4, 5, 6}, {7, 8, 9}});

  Array cpu = input;
  ridge_accentuate(cpu, 0.f, 1);

  EXPECT_TRUE(assert_almost_equal(cpu, input));
}

TEST(RidgeAccentuate, RidgeSharpeningOnGaussianPeak)
{
  // A symmetric 1D ridge along X (constant in Y)
  // z(x, y) = 1 - (x - 0.5)^2 * 4
  const int nx = 33;
  const int ny = 33;
  Array     input(glm::ivec2(nx, ny));

  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      float x = (float)i / (float)(nx - 1);
      float dx = x - 0.5f;
      input(i, j) = std::max(0.f, 1.f - 4.f * dx * dx);
    }

  Array res = input;
  ridge_accentuate(res, 0.1f, 1);

  // Peak at center should remain at max height (1.0)
  EXPECT_NEAR(res(nx / 2, ny / 2), 1.f, 1e-3f);

  // Ridge accentuation sharpens the peak by steepening/dropping the slopes
  int left_slope = nx / 2 - 4;
  EXPECT_LT(res(left_slope, ny / 2), input(left_slope, ny / 2));
}

TEST(RidgeAccentuate, MeshSizeIndependence)
{
  // Test that applying the filter at different grid resolutions (e.g. 65x65 vs
  // 257x257) produces identical elevation values at corresponding normalized
  // coordinates.
  auto make_pyramid = [](int n)
  {
    Array z(glm::ivec2(n, n));
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i)
      {
        float x = (float)i / (float)(n - 1);
        float y = (float)j / (float)(n - 1);
        z(i, j) = 1.f - std::abs(x - 0.5f) - std::abs(y - 0.5f);
      }
    return z;
  };

  const int n_low = 65;
  const int n_high = 257;

  Array z_low = make_pyramid(n_low);
  Array z_high = make_pyramid(n_high);

  float strength = 0.05f;

  ridge_accentuate(z_low, strength, 0);
  ridge_accentuate(z_high, strength, 0);

  // Sample at matching normalized points along the diagonal
  for (float u : {0.2f, 0.35f, 0.5f, 0.65f, 0.8f})
  {
    int i_low = (int)std::round(u * (n_low - 1));
    int i_high = (int)std::round(u * (n_high - 1));

    float val_low = z_low(i_low, i_low);
    float val_high = z_high(i_high, i_high);

    EXPECT_NEAR(val_low, val_high, 0.02f);
  }
}

TEST(RidgeAccentuate, MaskSupport)
{
  const int nx = 33;
  const int ny = 33;
  Array     input(glm::ivec2(nx, ny));
  Array     mask(glm::ivec2(nx, ny), 0.f);

  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      float x = (float)i / (float)(nx - 1);
      input(i, j) = 1.f - std::abs(x - 0.5f) * 2.f;
      if (j >= ny / 2) mask(i, j) = 1.f;
    }

  Array res = input;
  ridge_accentuate(res, &mask, 0.1f, 1);

  // Unmasked half (j < ny / 2) should remain unchanged
  for (int j = 0; j < ny / 4; ++j)
    for (int i = 0; i < nx; ++i)
      EXPECT_NEAR(res(i, j), input(i, j), 1e-5f);

  // Masked half should be modified on slopes
  int left_slope = nx / 2 - 4;
  EXPECT_LT(res(left_slope, 3 * ny / 4), input(left_slope, 3 * ny / 4));
}

TEST(RidgeAccentuate, ReverseFormsValleys)
{
  const int nx = 33;
  const int ny = 33;
  Array     input(glm::ivec2(nx, ny));

  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      float x = (float)i / (float)(nx - 1);
      input(i, j) = 1.f - std::abs(x - 0.5f) * 2.f;
    }

  Array res = input;
  ridge_accentuate(res, 0.1f, 1, true); // reverse = true

  // Reverse widens slopes, so slope points increase in elevation
  int left_slope = nx / 2 - 4;
  EXPECT_GT(res(left_slope, ny / 2), input(left_slope, ny / 2));
}

TEST(RidgeAccentuate_CPU_GPU, RandomFieldsEquivalence)
{
  std::mt19937                          rng(42);
  std::uniform_real_distribution<float> dist(0.f, 1.f);

  const int nx = 64;
  const int ny = 64;

  for (int test = 0; test < 5; ++test)
  {
    Array input(glm::ivec2(nx, ny));

    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        input(i, j) = dist(rng);

    Array cpu = input;
    Array gpu = input;

    float strength = 0.05f + 0.05f * dist(rng);
    int   ir = 1 + (int)(3 * dist(rng));

    ridge_accentuate(cpu, strength, ir);
    gpu::ridge_accentuate(gpu, strength, ir);

    EXPECT_TRUE(assert_almost_equal(cpu, gpu, 1e-2f));
  }
}
