#include "highmap/array.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/core.hpp"

#include <gtest/gtest.h>

using namespace hmap;

// ------------------------------------------------------------
// ABS
// ------------------------------------------------------------

TEST(MathTest, AbsBasic)
{
  Array a = Array({{-1.f, 0.f, 2.f}});
  Array out = abs(a);

  EXPECT_FLOAT_EQ(out(0, 0), 1.f);
  EXPECT_FLOAT_EQ(out(0, 1), 0.f);
  EXPECT_FLOAT_EQ(out(0, 2), 2.f);
}

TEST(MathTest, AbsSmoothBasic)
{
  Array a = Array({{-3.f}});
  float k = 4.f;

  Array out = abs_smooth(a, k);

  EXPECT_NEAR(out(0, 0), std::sqrt(9.f + 16.f), 1e-6f);
}

// ------------------------------------------------------------
// ALMOST UNIT IDENTITY
// ------------------------------------------------------------

TEST(MathTest, AlmostUnitIdentity)
{
  EXPECT_FLOAT_EQ(almost_unit_identity(0.f), 0.f);
  EXPECT_FLOAT_EQ(almost_unit_identity(1.f), 1.f);
}

// ------------------------------------------------------------
// TRIG / EXP
// ------------------------------------------------------------

TEST(MathTest, CosZero)
{
  Array a = Array({{0.f}});
  Array out = cos(a);

  EXPECT_NEAR(out(0, 0), 1.f, 1e-6f);
}

TEST(MathTest, ExpZero)
{
  Array a = Array({{0.f}});
  Array out = exp(a);

  EXPECT_NEAR(out(0, 0), 1.f, 1e-6f);
}

// ------------------------------------------------------------
// FAST APPROX
// ------------------------------------------------------------

TEST(MathTest, FastExpApprox)
{
  float x = 1.f;
  float approx = fast_exp(x);
  float exact = std::exp(x);

  EXPECT_NEAR(approx, exact, 0.2f); // coarse tolerance
}

TEST(MathTest, FastLogApprox)
{
  float x = 2.f;
  float approx = fast_log(x);
  float exact = std::log(x);

  EXPECT_NEAR(approx, exact, 0.2f);
}

// ------------------------------------------------------------
// GAUSSIAN
// ------------------------------------------------------------

TEST(MathTest, GaussianDecayAtZero)
{
  Array a = Array({{0.f}});
  Array out = gaussian_decay(a, 1.f);

  EXPECT_NEAR(out(0, 0), 1.f, 1e-3f);
}

// ------------------------------------------------------------
// HYPOT
// ------------------------------------------------------------

TEST(MathTest, HypotBasic)
{
  Array a = Array({{3.f}});
  Array b = Array({{4.f}});

  Array out = hypot(a, b);

  EXPECT_NEAR(out(0, 0), 5.f, 1e-6f);
}

// ------------------------------------------------------------
// LOG10
// ------------------------------------------------------------

TEST(MathTest, Log10Basic)
{
  Array a = Array({{10.f}});
  Array out = log10(a);

  EXPECT_NEAR(out(0, 0), 1.f, 1e-6f);
}

// ------------------------------------------------------------
// POW (optimized paths)
// ------------------------------------------------------------

TEST(MathTest, PowSquare)
{
  Array a = Array({{3.f}});
  Array out = pow(a, 2.f);

  EXPECT_FLOAT_EQ(out(0, 0), 9.f);
}

TEST(MathTest, PowSqrt)
{
  Array a = Array({{4.f}});
  Array out = pow(a, 0.5f);

  EXPECT_FLOAT_EQ(out(0, 0), 2.f);
}

// ------------------------------------------------------------
// LERP
// ------------------------------------------------------------

TEST(MathTest, LerpScalar)
{
  EXPECT_FLOAT_EQ(lerp(0.f, 10.f, 0.5f), 5.f);
}

TEST(MathTest, LerpArray)
{
  Array a = Array({{0.f}});
  Array b = Array({{10.f}});
  Array t = Array({{0.5f}});

  Array out = lerp(a, b, t);

  EXPECT_FLOAT_EQ(out(0, 0), 5.f);
}

// ------------------------------------------------------------
// SIGMOID
// ------------------------------------------------------------

TEST(MathTest, SigmoidCenter)
{
  float v = sigmoid(0.f, 1.f, 0.f, 1.f, 0.f);

  EXPECT_NEAR(v, 0.5f, 1e-6f);
}

// ------------------------------------------------------------
// SQRT SAFE
// ------------------------------------------------------------

TEST(MathTest, SqrtSafe)
{
  Array a = Array({{-1.f, 4.f}});
  Array out = sqrt_safe(a);

  EXPECT_FLOAT_EQ(out(0, 0), 0.f);
  EXPECT_FLOAT_EQ(out(0, 1), 2.f);
}

// ------------------------------------------------------------
// THRESHOLD
// ------------------------------------------------------------

TEST(MathTest, ThresholdBasic)
{
  EXPECT_FLOAT_EQ(threshold(0.f, 1.f, 2.f), 0.f);
  EXPECT_FLOAT_EQ(threshold(3.f, 1.f, 2.f), 1.f);
  EXPECT_FLOAT_EQ(threshold(1.5f, 1.f, 2.f), 0.5f);
}

TEST(MathTest, ThresholdArray)
{
  Array a = Array({{0.f, 1.5f, 3.f}});
  Array out = threshold(a, 1.f, 2.f);

  EXPECT_FLOAT_EQ(out(0, 0), 0.f);
  EXPECT_FLOAT_EQ(out(0, 1), 0.5f);
  EXPECT_FLOAT_EQ(out(0, 2), 1.f);
}

// ------------------------------------------------------------
// SMOOTHSTEP
// ------------------------------------------------------------

TEST(MathTest, Smoothstep3Basic)
{
  EXPECT_FLOAT_EQ(smoothstep3(0.f), 0.f);
  EXPECT_FLOAT_EQ(smoothstep3(1.f), 1.f);
}

// ------------------------------------------------------------
// IS_ZERO / IS_NON_ZERO
// ------------------------------------------------------------

TEST(MathTest, IsZero)
{
  Array a = Array({{0.f, 1.f}});
  Array out = is_zero(a);

  EXPECT_FLOAT_EQ(out(0, 0), 1.f);
  EXPECT_FLOAT_EQ(out(0, 1), 0.f);
}

TEST(MathTest, IsNonZero)
{
  Array a = Array({{0.f, 1.f}});
  Array out = is_non_zero(a);

  EXPECT_FLOAT_EQ(out(0, 0), 0.f);
  EXPECT_FLOAT_EQ(out(0, 1), 1.f);
}

// ------------------------------------------------------------
// CEIL DIV
// ------------------------------------------------------------

TEST(MathTest, CeilDiv)
{
  EXPECT_EQ(ceil_div(10, 3), 4);
  EXPECT_EQ(ceil_div(9, 3), 3);
}

// ------------------------------------------------------------
// TRIANGLE
// ------------------------------------------------------------

TEST(MathTest, TriangleBasic)
{
  EXPECT_FLOAT_EQ(triangle(0.f, 1.f, 2.f), 0.f);
  EXPECT_FLOAT_EQ(triangle(1.5f, 1.f, 2.f), 1.f);
  EXPECT_FLOAT_EQ(triangle(2.f, 1.f, 2.f), 0.f);
}

// ------------------------------------------------------------
// RVACHEV R-FUNCTIONS (R_MIN / R_MAX)
// ------------------------------------------------------------

TEST(MathTest, RMinScalar)
{
  EXPECT_FLOAT_EQ(r_min(0.f, 5.f), 0.f);
  EXPECT_FLOAT_EQ(r_min(5.f, 0.f), 0.f);
  EXPECT_FLOAT_EQ(r_min(0.f, 0.f), 0.f);
  EXPECT_NEAR(r_min(3.f, 4.f), 3.f + 4.f - 5.f, 1e-6f); // alpha = 0 -> 2.0f

  // alpha = 1 -> exact min
  EXPECT_FLOAT_EQ(r_min(3.f, 4.f, 1.f), 3.f);
  EXPECT_FLOAT_EQ(r_min(7.f, 2.f, 1.f), 2.f);

  // zero-set preserved for all alpha in (-1, 1]
  EXPECT_FLOAT_EQ(r_min(0.f, 8.f, 0.5f), 0.f);
  EXPECT_FLOAT_EQ(r_min(8.f, 0.f, -0.5f), 0.f);
}

TEST(MathTest, RMaxScalar)
{
  EXPECT_FLOAT_EQ(r_max(0.f, 5.f), 10.f);
  EXPECT_FLOAT_EQ(r_max(5.f, 0.f), 10.f);
  EXPECT_FLOAT_EQ(r_max(0.f, 0.f), 0.f);
  EXPECT_NEAR(r_max(3.f, 4.f), 3.f + 4.f + 5.f, 1e-6f); // alpha = 0 -> 12.0f

  // alpha = 1 -> exact max
  EXPECT_FLOAT_EQ(r_max(3.f, 4.f, 1.f), 4.f);
  EXPECT_FLOAT_EQ(r_max(7.f, 2.f, 1.f), 7.f);
}

TEST(MathTest, RMinRMaxArray)
{
  Array a = Array({{0.f, 3.f}});
  Array b = Array({{5.f, 4.f}});

  Array out_min = r_min(a, b);
  EXPECT_FLOAT_EQ(out_min(0, 0), 0.f);
  EXPECT_NEAR(out_min(0, 1), 2.f, 1e-6f);

  Array out_min_exact = r_min(a, b, 1.f);
  EXPECT_FLOAT_EQ(out_min_exact(0, 0), 0.f);
  EXPECT_FLOAT_EQ(out_min_exact(0, 1), 3.f);

  Array out_max = r_max(a, b);
  EXPECT_FLOAT_EQ(out_max(0, 0), 10.f);
  EXPECT_NEAR(out_max(0, 1), 12.f, 1e-6f);

  Array out_max_exact = r_max(a, b, 1.f);
  EXPECT_FLOAT_EQ(out_max_exact(0, 0), 5.f);
  EXPECT_FLOAT_EQ(out_max_exact(0, 1), 4.f);
}
