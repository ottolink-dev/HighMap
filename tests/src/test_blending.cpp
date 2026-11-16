#include "highmap/blending.hpp"
#include "highmap/internal/assert.hpp"
#include "highmap/primitives.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(Blending, BlendExclusion)
{
  // Formula: 0.5 - 2 * (a - 0.5) * (b - 0.5) = a + b - 2 * a * b
  Array a = Array({{0.0f, 0.5f}, {1.0f, 0.2f}});
  Array b = Array({{0.0f, 0.5f}, {0.0f, 0.8f}});

  Array res = blend_exclusion(a, b);

  // (0,0): 0 + 0 - 0 = 0
  EXPECT_NEAR(res(0, 0), 0.0f, 1e-5f);
  // (0,1): 0.5 + 0.5 - 2*0.25 = 0.5
  EXPECT_NEAR(res(0, 1), 0.5f, 1e-5f);
  // (1,0): 1 + 0 - 0 = 1.0
  EXPECT_NEAR(res(1, 0), 1.0f, 1e-5f);
  // (1,1): 0.2 + 0.8 - 2 * 0.16 = 1.0 - 0.32 = 0.68
  EXPECT_NEAR(res(1, 1), 0.68f, 1e-5f);

  // Symmetry check: exclusion(a, b) == exclusion(b, a)
  Array res_rev = blend_exclusion(b, a);
  EXPECT_TRUE(assert_almost_equal(res, res_rev, 1e-5f));
}

TEST(Blending, BlendNegate)
{
  // Formula: a < b ? a : 2*b - a
  Array a = Array({{0.2f, 0.8f}, {0.5f, 0.9f}});
  Array b = Array({{0.5f, 0.5f}, {0.5f, 0.1f}});

  Array res = blend_negate(a, b);

  // (0,0): a < b (0.2 < 0.5) => 0.2
  EXPECT_NEAR(res(0, 0), 0.2f, 1e-5f);
  // (0,1): a >= b (0.8 >= 0.5) => 2*0.5 - 0.8 = 0.2
  EXPECT_NEAR(res(0, 1), 0.2f, 1e-5f);
  // (1,0): a == b (0.5 == 0.5) => 2*0.5 - 0.5 = 0.5
  EXPECT_NEAR(res(1, 0), 0.5f, 1e-5f);
  // (1,1): a >= b (0.9 >= 0.1) => 2*0.1 - 0.9 = -0.7
  EXPECT_NEAR(res(1, 1), -0.7f, 1e-5f);
}

TEST(Blending, BlendOverlay)
{
  // Formula: a < 0.5 ? 2*a*b : 1 - 2*(1-a)*(1-b)
  Array a = Array({{0.25f, 0.75f}, {0.0f, 1.0f}});
  Array b = Array({{0.40f, 0.60f}, {0.5f, 0.5f}});

  Array res = blend_overlay(a, b);

  // (0,0): a < 0.5 => 2 * 0.25 * 0.40 = 0.2
  EXPECT_NEAR(res(0, 0), 0.2f, 1e-5f);
  // (0,1): a >= 0.5 => 1 - 2 * (1 - 0.75) * (1 - 0.60) = 1 - 2 * 0.25 * 0.40 =
  // 0.8
  EXPECT_NEAR(res(0, 1), 0.8f, 1e-5f);
  // (1,0): a == 0.0 < 0.5 => 0
  EXPECT_NEAR(res(1, 0), 0.0f, 1e-5f);
  // (1,1): a == 1.0 >= 0.5 => 1 - 2 * 0 * 0.5 = 1.0
  EXPECT_NEAR(res(1, 1), 1.0f, 1e-5f);
}

TEST(Blending, BlendSoft)
{
  // Formula: (1 - a) * a * b + a * (1 - (1 - a) * (1 - b))
  Array a = Array({{0.0f, 1.0f}, {0.5f, 0.2f}});
  Array b = Array({{0.5f, 0.5f}, {0.5f, 0.8f}});

  Array res = blend_soft(a, b);

  // (0,0): a = 0 => 0
  EXPECT_NEAR(res(0, 0), 0.0f, 1e-5f);
  // (0,1): a = 1 => 1
  EXPECT_NEAR(res(0, 1), 1.0f, 1e-5f);
  // (1,0): a = 0.5, b = 0.5 => 0.5 * 0.5 * 0.5 + 0.5 * (1 - 0.5 * 0.5) = 0.125
  // + 0.5 * 0.75 = 0.5
  EXPECT_NEAR(res(1, 0), 0.5f, 1e-5f);
  // (1,1): a = 0.2, b = 0.8 => 0.8 * 0.2 * 0.8 + 0.2 * (1 - 0.8 * 0.2) = 0.128
  // + 0.2 * 0.84 = 0.128 + 0.168 = 0.296
  EXPECT_NEAR(res(1, 1), 0.296f, 1e-5f);
}

TEST(Blending, Mixer)
{
  glm::ivec2 shape = {16, 16};
  Array      a1(shape, 1.0f);
  Array      a2(shape, 5.0f);
  Array      a3(shape, 10.0f);

  std::vector<const Array *> arrays = {&a1, &a2, &a3};

  // t = 0 -> only array 0 (weight 1)
  Array t0(shape, 0.0f);
  Array res0 = mixer(t0, arrays);
  EXPECT_NEAR(res0(0, 0), 1.0f, 1e-5f);
  EXPECT_NEAR(res0(8, 8), 1.0f, 1e-5f);

  // t = 0.5 -> center array 1 (r0 = 1/2 = 0.5, ta = 1, ts = 1)
  Array t1(shape, 0.5f);
  Array res1 = mixer(t1, arrays);
  EXPECT_NEAR(res1(0, 0), 5.0f, 1e-5f);

  // t = 1.0 -> only array 2 (r0 = 2/2 = 1.0, ta = 1, ts = 1)
  Array t2(shape, 1.0f);
  Array res2 = mixer(t2, arrays);
  EXPECT_NEAR(res2(0, 0), 10.0f, 1e-5f);

  // Test with gain factor != 1.0
  Array res_gain = mixer(t0, arrays, 1.5f);
  EXPECT_GT(res_gain(0, 0), 0.0f);
}

TEST(Blending, Transfer)
{
  glm::ivec2 shape = {32, 32};
  Array      source(shape, 10.0f);
  Array      target(shape, 2.0f);

  // When source is constant, high-pass filter w = 0
  Array res = transfer(source, target, 4, 1.0f, false);
  EXPECT_TRUE(assert_almost_equal(res, target, 1e-5f));

  // With prefiltering on constant target, should also equal target
  Array res_pre = transfer(source, target, 4, 1.0f, true);
  EXPECT_TRUE(assert_almost_equal(res_pre, target, 1e-5f));
}

TEST(Blending, BlendGradients)
{
  glm::ivec2 shape = {32, 32};
  Array      a(shape, 0.0f);
  Array      b(shape, 1.0f);

  // Should run without error and produce an array bounded within [0, 1]
  Array res = blend_gradients(a, b, 2);
  EXPECT_EQ(res.shape.x, shape.x);
  EXPECT_EQ(res.shape.y, shape.y);
  EXPECT_GE(res.min(), -1e-4f);
  EXPECT_LE(res.max(), 1.0f + 1e-4f);
}

TEST(Blending, BlendPowerLaw)
{
  glm::ivec2 shape = {4, 4};
  Array      f1(shape, 2.0f);
  Array      f2(shape, 8.0f);

  // alpha = 0 -> arithmetic mean: (2 + 8) / 2 = 5.0
  Array res_alpha0 = blend_power_law(f1, f2, 0.0f);
  EXPECT_NEAR(res_alpha0(0, 0), 5.0f, 1e-5f);
  EXPECT_NEAR(res_alpha0(2, 2), 5.0f, 1e-5f);

  // alpha = 1 -> (2^2 + 8^2) / (2^1 + 8^1) = (4 + 64) / (2 + 8) = 68 / 10 = 6.8
  Array res_alpha1 = blend_power_law(f1, f2, 1.0f);
  EXPECT_NEAR(res_alpha1(0, 0), 6.8f, 1e-5f);

  // alpha large -> approaches max(f_i) = 8.0
  // for alpha = 10: (2^11 + 8^11) / (2^10 + 8^10) ~= 8.0
  Array res_alpha10 = blend_power_law(f1, f2, 10.0f);
  EXPECT_NEAR(res_alpha10(0, 0), 8.0f, 1e-2f);

  // Test with vector of 3 arrays
  Array                      f3(shape, 5.0f);
  std::vector<const Array *> arrays = {&f1, &f2, &f3};

  // alpha = 0 -> (2 + 8 + 5) / 3 = 5.0
  Array res_vec0 = blend_power_law(arrays, 0.0f);
  EXPECT_NEAR(res_vec0(0, 0), 5.0f, 1e-5f);

  // alpha = 1 -> (2^2 + 8^2 + 5^2) / (2 + 8 + 5) = (4 + 64 + 25) / 15 = 93 / 15
  // = 6.2
  Array res_vec1 = blend_power_law(arrays, 1.0f);
  EXPECT_NEAR(res_vec1(0, 0), 6.2f, 1e-5f);
}
