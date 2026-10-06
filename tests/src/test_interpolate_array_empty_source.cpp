#include "highmap.hpp"

#include <gtest/gtest.h>

// Interpolating from a zero-sized source used to be undefined behaviour: the
// index arithmetic reduces to std::clamp(v, 0, shape.x - 1) == std::clamp(v, 0,
// -1), whose precondition (lo <= hi) is violated and which in practice yields
// -1, so source(-1, -1) reads before the buffer and segfaults.
//
// Reached from Hesiod by opening any project whose Brush node deserialized to
// an empty array (otto-link/Hesiod#658): Array::resample_to_shape_bilinear on a
// 0x0 source. An empty source has no meaningful interpolation, so the contract
// is to leave the target at its zero-filled initial state.

namespace
{

void expect_all_zero(const hmap::Array &z)
{
  for (int j = 0; j < z.shape.y; ++j)
    for (int i = 0; i < z.shape.x; ++i)
      EXPECT_EQ(z(i, j), 0.f);
}

} // namespace

TEST(InterpolateArrayEmptySource, BilinearLeavesTargetZeroed)
{
  hmap::Array source; // default-constructed: 0x0
  hmap::Array target(glm::ivec2(16, 16));

  ASSERT_EQ(source.shape.x, 0);
  ASSERT_EQ(source.shape.y, 0);

  hmap::interpolate_array_bilinear(source, target);

  EXPECT_EQ(target.shape.x, 16);
  EXPECT_EQ(target.shape.y, 16);
  expect_all_zero(target);
}

TEST(InterpolateArrayEmptySource, NearestLeavesTargetZeroed)
{
  hmap::Array source;
  hmap::Array target(glm::ivec2(16, 16));

  hmap::interpolate_array_nearest(source, target);

  expect_all_zero(target);
}

TEST(InterpolateArrayEmptySource, BicubicLeavesTargetZeroed)
{
  hmap::Array source;
  hmap::Array target(glm::ivec2(16, 16));

  hmap::interpolate_array_bicubic(source, target);

  expect_all_zero(target);
}

TEST(InterpolateArrayEmptySource, SingleDegenerateDimensionIsHandled)
{
  // only one dimension collapsed — still degenerate for the index arithmetic
  hmap::Array source(glm::ivec2(0, 8));
  hmap::Array target(glm::ivec2(16, 16));

  hmap::interpolate_array_bilinear(source, target);

  expect_all_zero(target);
}

TEST(InterpolateArrayEmptySource, ResampleToShapeIsSafe)
{
  // the path Hesiod's Brush node actually took
  hmap::Array source;

  hmap::Array out_bilinear = source.resample_to_shape_bilinear(
      glm::ivec2(32, 32));
  EXPECT_EQ(out_bilinear.shape.x, 32);
  EXPECT_EQ(out_bilinear.shape.y, 32);
  expect_all_zero(out_bilinear);

  hmap::Array out_nearest = source.resample_to_shape_nearest(
      glm::ivec2(32, 32));
  EXPECT_EQ(out_nearest.shape.x, 32);
  EXPECT_EQ(out_nearest.shape.y, 32);
  expect_all_zero(out_nearest);
}

TEST(InterpolateArrayEmptySource, NonEmptySourceStillInterpolates)
{
  // guard must not short-circuit the normal path
  hmap::Array source(glm::ivec2(2, 2));
  source(0, 0) = 0.f;
  source(1, 0) = 1.f;
  source(0, 1) = 1.f;
  source(1, 1) = 2.f;

  hmap::Array target(glm::ivec2(8, 8));
  hmap::interpolate_array_bilinear(source, target);

  EXPECT_GT(target.max(), 0.f);
}

TEST(InterpolateArrayEmptySource, ResampleBicubicWarpEmptySourceIsSafe)
{
  hmap::Array source;
  hmap::Array out = hmap::resample_bicubic_warp(source,
                                                glm::ivec2(32, 32),
                                                0.35f,
                                                42u);
  EXPECT_EQ(out.shape.x, 0);
  EXPECT_EQ(out.shape.y, 0);
}

TEST(InterpolateBspline, BsplineWeightsProperties)
{
  float w[4];
  for (float t = 0.f; t < 1.f; t += 0.1f)
  {
    hmap::bspline_weights(t, w);
    float sum = w[0] + w[1] + w[2] + w[3];
    EXPECT_NEAR(sum, 1.f, 1e-5f);
    for (int i = 0; i < 4; ++i)
    {
      EXPECT_GE(w[i], 0.f);
      EXPECT_LE(w[i], 1.f);
    }
  }
}

TEST(InterpolateArray, ResampleBicubicWarp)
{
  hmap::Array src(glm::ivec2(8, 8), 1.f);
  hmap::Array out = hmap::resample_bicubic_warp(src, glm::ivec2(16, 16), 0.f);
  EXPECT_EQ(out.shape.x, 16);
  EXPECT_EQ(out.shape.y, 16);
  EXPECT_NEAR(out.mean(), 1.f, 1e-5f);

  hmap::Array out_warp = hmap::resample_bicubic_warp(src,
                                                     glm::ivec2(16, 16),
                                                     0.5f,
                                                     123u);
  EXPECT_EQ(out_warp.shape.x, 16);
  EXPECT_EQ(out_warp.shape.y, 16);
  EXPECT_NEAR(out_warp.mean(), 1.f, 1e-5f);
}
