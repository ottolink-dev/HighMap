#include "highmap/filters.hpp"
#include "highmap/internal/assert.hpp"
#include "highmap/primitives.hpp"
#include "highmap/statistics.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(AutoCorrLengthScale, FlatArrayReturnsMaxDimension)
{
  Array input = Array({{1, 1, 1}, {1, 1, 1}});
  float sc = autocorr_length_scale(input);

  EXPECT_EQ(sc, std::max(input.shape.x, input.shape.y));
}

TEST(AutoCorrLengthScale, RandomNoiseHasShortLengthScale)
{
  Array input = white({64, 64}, 0.f, 1.f, 42);
  float sc = autocorr_length_scale(input);

  // should be very small (almost uncorrelated), minimum gap is 2
  // cells
  EXPECT_LE(sc, 2.f);
}

TEST(AutoCorrLengthScale, SmoothFieldHasLargeLengthScale)
{
  Array input = white({64, 64}, 0.f, 1.f, 42);
  smooth_cpulse(input, 16);

  float sc = autocorr_length_scale(input);

  // scale larger than the half-width of the filter
  EXPECT_GT(sc, 8.f);
}

TEST(AutoCorrLengthScale, IncreasesWithSmoothing)
{
  Array input = white({64, 64}, 0.f, 1.f, 42);
  Array a = input;
  Array b = input;

  smooth_cpulse(a, 1);
  smooth_cpulse(b, 5);

  float sc_a = autocorr_length_scale(a);
  float sc_b = autocorr_length_scale(b);

  EXPECT_GT(sc_b, sc_a);
}

TEST(AutoCorrLengthScale, LinearGradientHighCorrelation)
{
  Array input = Array(glm::ivec2(32, 32));

  for (int j = 0; j < input.shape.y; ++j)
    for (int i = 0; i < input.shape.x; ++i)
      input(i, j) = float(i); // strong correlation along x

  float sc = autocorr_length_scale(input);

  EXPECT_GT(sc, 10.f);
}

TEST(AutoCorrLengthScale, PeriodicSignal)
{
  const int n = 64;
  Array     input = Array(glm::ivec2(n, n));

  float period = 8.f;

  for (int j = 0; j < n; ++j)
    for (int i = 0; i < n; ++i)
      input(i, j) = std::sin(2.f * M_PI * i / period);

  float sc = autocorr_length_scale(input);

  // expected scale ~ period
  EXPECT_GT(sc, 0.6f * period);
  EXPECT_LT(sc, 1.4f * period);
}

TEST(AutoCorrLengthScale, InvariantToOffset)
{
  Array input = white({64, 64}, 0.f, 1.f, 42);
  float sc1 = autocorr_length_scale(input);
  float sc2 = autocorr_length_scale(10.f + input);

  EXPECT_NEAR(sc1, sc2, 1e-4f);
}

TEST(AutoCorrLengthScale, InvariantToScaling)
{
  Array input = white({64, 64}, 0.f, 1.f, 42);
  float sc1 = autocorr_length_scale(input);
  float sc2 = autocorr_length_scale(5.f * input);

  EXPECT_NEAR(sc1, sc2, 1e-4f);
}

TEST(AutoCorrLengthScale, MaxLagFractionLimitsResult)
{
  Array input = 10.f + white({64, 64}, 0.f, 1.f, 42);
  smooth_cpulse(input, 10);

  float sc1 = autocorr_length_scale(input, 0.2f);
  float sc2 = autocorr_length_scale(input, 0.5f);

  // should be the same
  EXPECT_NEAR(sc1, sc2, 1e-4f);
}

TEST(AutoCorrLengthScaleAxial, DetectsAnisotropy)
{
  const int nx = 64;
  const int ny = 64;
  float     lambda = 16.f;
  Array     input = Array(glm::ivec2(nx, ny));

  // x-direction
  {
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        input(i, j) = std::sin(2.f * M_PI * i / nx / lambda);

    glm::vec2 sc = autocorr_length_scale_axial(input, 0.5f);
    EXPECT_NEAR(sc.x, lambda, 2.f);
    EXPECT_NEAR(sc.y, 32.f, 2.f); // cap by search radius
  }

  // y-direction
  {
    float lambda = 16.f;

    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        input(i, j) = std::sin(2.f * M_PI * j / ny / lambda);

    glm::vec2 sc = autocorr_length_scale_axial(input, 0.5f);
    EXPECT_NEAR(sc.x, 32.f, 2.f);
    EXPECT_NEAR(sc.y, lambda, 2.f);
  }
}

TEST(AutoCorrLengthScaleAxial, MatchesScalarForIsotropicField)
{
  Array input = white({64, 64}, 0.f, 1.f, 42);

  // make it isotropic via smoothing
  smooth_cpulse(input, 3);

  float     scalar_sc = autocorr_length_scale(input, 0.4f);
  glm::vec2 axial_sc = autocorr_length_scale_axial(input, 0.4f);

  // both directions should be close to scalar estimate
  EXPECT_NEAR(axial_sc.x, scalar_sc, 2.f);
  EXPECT_NEAR(axial_sc.y, scalar_sc, 2.f);

  // and close to each other
  EXPECT_NEAR(axial_sc.x, axial_sc.y, 2.f);
}
