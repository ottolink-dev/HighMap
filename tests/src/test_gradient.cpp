#include "highmap/gradient.hpp"
#include "highmap/internal/assert.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(Gradient, GradientX_LinearRamp)
{
  // f(x,y) = x → df/dx = 1 everywhere
  Array input = Array({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}});
  Array dx = gradient_x(input);

  for (int j = 0; j < input.shape.y; ++j)
    for (int i = 0; i < input.shape.x; ++i)
      EXPECT_NEAR(dx(i, j), 1.f, 1e-5);
}

TEST(Gradient, GradientY_LinearRamp)
{
  // f(x,y) = y → df/dy = 1 everywhere
  Array input = Array({{0, 1, 2}, {0, 1, 2}, {0, 1, 2}});
  Array dy = gradient_y(input);

  for (int j = 0; j < input.shape.y; ++j)
    for (int i = 0; i < input.shape.x; ++i)
      EXPECT_NEAR(dy(i, j), 1.f, 1e-5);
}

TEST(Gradient, GradientNorm_ConstantField)
{
  Array input = Array({{5, 5}, {5, 5}});
  Array g = gradient_norm(input);

  EXPECT_TRUE(assert_almost_equal(g, Array(input.shape)));
}

TEST(Gradient, GradientNorm_SimpleSlope)
{
  // f(x,y) = x → |grad| = 1
  Array input = Array({{0, 1, 2}, {0, 1, 2}, {0, 1, 2}});
  Array g = gradient_norm(input);

  for (int j = 0; j < input.shape.y; ++j)
    for (int i = 0; i < input.shape.x; ++i)
      EXPECT_NEAR(g(i, j), 1.f, 1e-5);
}

TEST(Gradient, GradientNorm_OperatorsConsistency)
{
  Array input = Array({{0, 1, 2}, {1, 2, 3}, {2, 3, 4}});

  Array p = gradient_norm_prewitt(input);
  Array s = gradient_norm_sobel(input);
  Array c = gradient_norm_scharr(input);

  // same trend → values should be close
  EXPECT_TRUE(assert_almost_equal(p, s, 1e-5f));
  EXPECT_TRUE(assert_almost_equal(s, c, 1e-5f));
}

TEST(Gradient, GradientAngle_SlopeConsistency)
{
  Array input(glm::ivec2(3, 3));
  int   n_samples = 32;

  for (int k = 0; k < n_samples; ++k)
  {
    // avoid M_PI and -M_PI
    float alpha = 2.f * M_PI * float(k + 1) / float(n_samples + 1) - M_PI;

    for (int j = 0; j < input.shape.y; ++j)
      for (int i = 0; i < input.shape.x; ++i)
        input(i, j) = std::cos(alpha) * i + std::sin(alpha) * j;

    Array g = gradient_angle(input);

    for (int j = 0; j < input.shape.y; ++j)
      for (int i = 0; i < input.shape.x; ++i)
        EXPECT_NEAR(g(i, j), alpha, 1e-5);
  }
}

TEST(Gradient, GradientAngle_SlopeX)
{
  // f(x,y)=x → angle = 0
  Array input = Array({{0, 0, 0}, {1, 1, 1}, {2, 2, 2}});
  Array angle = gradient_angle(input, false);

  for (float v : angle.vector)
    EXPECT_NEAR(v, 0.f, 1e-5);
}

TEST(Gradient, GradientAngle_SlopeY)
{
  // f(x,y)=y → angle = 0
  Array input = Array({{0, 1, 2}, {0, 1, 2}, {0, 1, 2}});
  Array angle = gradient_angle(input, false);

  for (float v : angle.vector)
    EXPECT_NEAR(v, M_PI_2, 1e-5);
}

TEST(Gradient, Divergence_ZeroField)
{
  Array dx = Array({{0, 0}, {0, 0}});
  Array dy = Array({{0, 0}, {0, 0}});
  Array div = divergence_from_gradients(dx, dy);

  EXPECT_TRUE(assert_almost_equal(div, Array(dx.shape)));
}

TEST(Gradient, Laplacian_ConstantField)
{
  Array input = Array({{3, 3}, {3, 3}});
  Array lap = laplacian(input);

  EXPECT_TRUE(assert_almost_equal(lap, Array(input.shape)));
}

TEST(Gradient, Laplacian_Quadratic)
{
  // f(x,y) = x^2 + y^2 → Laplacian = 4
  Array input(glm::ivec2(3, 3));

  for (int j = 0; j < 3; ++j)
    for (int i = 0; i < 3; ++i)
      input(i, j) = i * i + j * j;

  Array lap = laplacian(input);

  EXPECT_NEAR(lap(1, 1), 4.f, 1e-5);
}

TEST(Gradient, GradientTalus_Step)
{
  Array input = Array({{0, 0, 1, 1}});
  Array talus = gradient_talus(input);

  // jump of 1 should be detected
  EXPECT_GE(talus(1, 0), 1.f);
  EXPECT_GE(talus(2, 0), 1.f);
}

TEST(Gradient, TalusJumpMask_BinaryResponse)
{
  Array input = Array({{0, 0, 10, 10}});
  Array mask = talus_jump_mask(input, 1.f, 1.f);

  // high gradient → mask near 1
  EXPECT_GT(mask(1, 0), 0.5f);
  EXPECT_GT(mask(2, 0), 0.5f);
}
