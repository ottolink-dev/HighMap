#include "highmap/dbg/assert.hpp"
#include "highmap/gradient.hpp"
#include "highmap/selector.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(SelectAngle, WrapAcrossZeroNorth)
{
  // create planes with known slopes (aspects)
  // aspect is downhill direction = atan2(-dy, -dx)
  auto make_plane = [](float alpha_deg)
  {
    Array z(glm::ivec2(5, 5));
    float alpha_rad = alpha_deg * float(M_PI) / 180.f;
    for (int j = 0; j < z.shape.y; ++j)
    {
      for (int i = 0; i < z.shape.x; ++i)
      {
        z(i, j) = -(std::cos(alpha_rad) * float(i) +
                    std::sin(alpha_rad) * float(j));
      }
    }
    return z;
  };

  // terrain with aspect at 5° (near 0°)
  Array z_5deg = make_plane(5.f);

  // select angle at 350° +/- 20° (should cover [330°, 370°] which includes 5°)
  Array sel = select_angle(z_5deg, 350.f, 20.f);

  // distance from 350° to 5° is 15°, which is inside tolerance 20°
  // smoothstep falloff at d/sigma = 15/20 = 0.75
  float expected_val = 1.f -
                       (3.f * 0.75f * 0.75f - 2.f * 0.75f * 0.75f * 0.75f);
  EXPECT_NEAR(sel(2, 2), expected_val, 1e-4);

  // select angle at 10° +/- 20° on aspect 355°
  Array z_355deg = make_plane(355.f);
  Array sel2 = select_angle(z_355deg, 10.f, 20.f);
  EXPECT_NEAR(sel2(2, 2), expected_val, 1e-4);
}

TEST(SelectAngle, BoundaryEquivalence0And360)
{
  Array z(glm::ivec2(8, 8));
  for (int j = 0; j < z.shape.y; ++j)
    for (int i = 0; i < z.shape.x; ++i)
      z(i, j) = std::sin(float(i)) * std::cos(float(j));

  Array sel_0 = select_angle(z, 0.f, 30.f);
  Array sel_360 = select_angle(z, 360.f, 30.f);
  Array sel_neg360 = select_angle(z, -360.f, 30.f);
  Array sel_720 = select_angle(z, 720.f, 30.f);

  EXPECT_TRUE(assert_almost_equal(sel_0, sel_360, 1e-5f));
  EXPECT_TRUE(assert_almost_equal(sel_0, sel_neg360, 1e-5f));
  EXPECT_TRUE(assert_almost_equal(sel_0, sel_720, 1e-5f));
}

TEST(SelectAngle, ExactMatchAndOutOfTolerance)
{
  auto make_plane = [](float alpha_deg)
  {
    Array z(glm::ivec2(5, 5));
    float alpha_rad = alpha_deg * float(M_PI) / 180.f;
    for (int j = 0; j < z.shape.y; ++j)
      for (int i = 0; i < z.shape.x; ++i)
        z(i, j) = -(std::cos(alpha_rad) * float(i) +
                    std::sin(alpha_rad) * float(j));
    return z;
  };

  Array z_90deg = make_plane(90.f);

  // exact match -> should be 1.0
  Array sel_exact = select_angle(z_90deg, 90.f, 30.f);
  EXPECT_NEAR(sel_exact(2, 2), 1.f, 1e-4);

  // completely outside tolerance -> should be 0.0
  Array sel_outside = select_angle(z_90deg, 150.f, 30.f);
  EXPECT_NEAR(sel_outside(2, 2), 0.f, 1e-4);

  // zero or negative sigma -> should be 0.0
  Array sel_zero_sigma = select_angle(z_90deg, 90.f, 0.f);
  EXPECT_NEAR(sel_zero_sigma(2, 2), 0.f, 1e-4);
}

TEST(SelectAngle, FlatTerrainExclusion)
{
  // perfectly flat array
  Array z_flat(glm::ivec2(6, 6), 5.f);

  // should return all 0s even when target angle is 0°
  Array sel = select_angle(z_flat, 0.f, 30.f);
  for (int j = 0; j < sel.shape.y; ++j)
    for (int i = 0; i < sel.shape.x; ++i)
      EXPECT_EQ(sel(i, j), 0.f);
}

TEST(SelectAngle, ValidationEmpty)
{
  Array empty;
  Array sel = select_angle(empty, 0.f, 10.f);
  EXPECT_EQ(sel.size(), 0);
}

// --- SelectRange Unit Tests

TEST(SelectRange, InsideRange)
{
  Array a(glm::ivec2(5, 5), 0.5f);
  Array sel = select_range(a, 0.3f, 0.7f, 0.1f);

  for (int j = 0; j < sel.shape.y; ++j)
    for (int i = 0; i < sel.shape.x; ++i)
      EXPECT_FLOAT_EQ(sel(i, j), 1.0f);
}

TEST(SelectRange, OutsideRange)
{
  Array a(glm::ivec2(5, 5), 0.1f);
  // Range [0.4, 0.8] with width 0.1 -> transition lower is [0.3, 0.4] -> 0.1 is
  // strictly outside
  Array sel = select_range(a, 0.4f, 0.8f, 0.1f);

  for (int j = 0; j < sel.shape.y; ++j)
    for (int i = 0; i < sel.shape.x; ++i)
      EXPECT_FLOAT_EQ(sel(i, j), 0.0f);

  // Upper outside test
  Array a_high(glm::ivec2(5, 5), 0.95f);
  Array sel_high = select_range(a_high, 0.4f, 0.8f, 0.1f);
  for (int j = 0; j < sel_high.shape.y; ++j)
    for (int i = 0; i < sel_high.shape.x; ++i)
      EXPECT_FLOAT_EQ(sel_high(i, j), 0.0f);
}

TEST(SelectRange, SmoothTransition)
{
  // vmin = 0.4, vmax = 0.6, width = 0.2
  // Lower transition: [0.2, 0.4] -> at 0.3 (midpoint of transition):
  // left = threshold_smooth(0.3, 0.2, 0.4) = smoothstep3(0.5) = 0.5
  // right = 1.0 - threshold_smooth(0.3, 0.6, 0.8) = 1.0 - 0 = 1.0
  // result = 0.5
  Array a(glm::ivec2(3, 3), 0.3f);
  Array sel = select_range(a, 0.4f, 0.6f, 0.2f);
  EXPECT_NEAR(sel(1, 1), 0.5f, 1e-5f);

  // Upper transition midpoint: at 0.7:
  // left = 1.0
  // right = 1.0 - threshold_smooth(0.7, 0.6, 0.8) = 1.0 - smoothstep3(0.5) =
  // 0.5
  Array a2(glm::ivec2(3, 3), 0.7f);
  Array sel2 = select_range(a2, 0.4f, 0.6f, 0.2f);
  EXPECT_NEAR(sel2(1, 1), 0.5f, 1e-5f);
}

TEST(SelectRange, ZeroWidthStep)
{
  Array a(glm::ivec2(3, 3), 0.5f);
  Array sel_inside = select_range(a, 0.4f, 0.6f, 0.0f);
  EXPECT_FLOAT_EQ(sel_inside(1, 1), 1.0f);

  Array sel_outside = select_range(a, 0.6f, 0.8f, 0.0f);
  EXPECT_FLOAT_EQ(sel_outside(1, 1), 0.0f);
}

TEST(SelectRange, InvertedVminVmax)
{
  Array a(glm::ivec2(3, 3), 0.5f);
  Array sel1 = select_range(a, 0.3f, 0.7f, 0.1f);
  Array sel2 = select_range(a, 0.7f, 0.3f, 0.1f);
  EXPECT_FLOAT_EQ(sel1(1, 1), sel2(1, 1));
}

TEST(SelectRange, ValidationEmpty)
{
  Array empty;
  Array sel = select_range(empty, 0.2f, 0.8f, 0.1f);
  EXPECT_EQ(sel.size(), 0);
}
