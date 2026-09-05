#include "highmap.hpp"

#include <gtest/gtest.h>

namespace
{

// axis-aligned closed square contour, centre (cx, cy), half-size half
hmap::Path square(float cx, float cy, float half)
{
  std::vector<float> x = {cx - half, cx + half, cx + half, cx - half};
  std::vector<float> y = {cy - half, cy - half, cy + half, cy + half};
  hmap::Path         p(x, y);
  p.set_closed(true);
  return p;
}

bool is_empty(const hmap::Array &a)
{
  return a.size() == 0;
}

} // namespace

TEST(ElevationFromContours, InvalidInputsReturnEmpty)
{
  glm::ivec2 shape = {32, 32};

  // no contours
  EXPECT_TRUE(is_empty(hmap::elevation_from_contours(shape, {}, {})));

  // contours / elevations size mismatch
  EXPECT_TRUE(is_empty(hmap::elevation_from_contours(shape,
                                                     {square(0.5f, 0.5f, 0.2f)},
                                                     {0.f, 1.f})));

  // degenerate contour (2 points)
  hmap::Path two({0.f, 1.f}, {0.f, 1.f});
  EXPECT_TRUE(is_empty(hmap::elevation_from_contours(shape, {two}, {0.f})));

  // probability map with the wrong shape
  hmap::Array bad_p({16, 16}, 0.5f);
  EXPECT_TRUE(is_empty(hmap::elevation_from_contours(shape,
                                                     {square(0.5f, 0.5f, 0.2f)},
                                                     {0.f},
                                                     &bad_p)));
}

TEST(ElevationFromContours, ContourPixelsAreExact)
{
  glm::ivec2              shape = {64, 64};
  std::vector<hmap::Path> c = {square(0.5f, 0.5f, 0.4f),
                               square(0.5f, 0.5f, 0.2f)};
  std::vector<float>      h = {0.2f, 0.6f};

  hmap::Array z = hmap::elevation_from_contours(shape, c, h, nullptr, 0.f);
  ASSERT_EQ(z.shape, shape);

  // outline pixel of the outer square: x = 0.1 -> i = round(0.1 * 63) = 6
  EXPECT_NEAR(z(6, 32), 0.2f, 1e-6f);
  // outline pixel of the inner square: x = 0.3 -> i = round(0.3 * 63) = 19
  EXPECT_NEAR(z(19, 32), 0.6f, 1e-6f);
  // top edge of the outer square: y = 0.9 -> j = round(0.9 * 63) = 57
  EXPECT_NEAR(z(32, 57), 0.2f, 1e-6f);
}

TEST(ElevationFromContours, BetweenNestedContoursValuesAreBounded)
{
  glm::ivec2              shape = {64, 64};
  std::vector<hmap::Path> c = {square(0.5f, 0.5f, 0.4f),
                               square(0.5f, 0.5f, 0.2f)};
  std::vector<float>      h = {0.2f, 0.6f};

  for (float r : {0.f, 1.f})
  {
    hmap::Array z = hmap::elevation_from_contours(shape, c, h, nullptr, r, 3);

    // ring between the squares along the middle row (outer outline at i = 6,
    // inner outline at i = 19)
    for (int i = 7; i < 19; ++i)
    {
      EXPECT_GE(z(i, 32), 0.2f - 1e-5f) << "r=" << r << " i=" << i;
      EXPECT_LE(z(i, 32), 0.6f + 1e-5f) << "r=" << r << " i=" << i;
    }

    // strictly inside the ring the terrain must have left the contour levels
    EXPECT_GT(z(12, 32), 0.2f + 1e-3f) << "r=" << r;
    EXPECT_LT(z(12, 32), 0.6f - 1e-3f) << "r=" << r;

    // deterministic case: monotone along the row from outer to inner contour
    if (r == 0.f)
    {
      for (int i = 6; i < 19; ++i)
        EXPECT_LE(z(i, 32), z(i + 1, 32) + 1e-6f) << "i=" << i;
    }
  }
}

TEST(ElevationFromContours, DeterministicModeIgnoresSeed)
{
  glm::ivec2              shape = {48, 48};
  std::vector<hmap::Path> c = {square(0.5f, 0.5f, 0.4f),
                               square(0.5f, 0.5f, 0.2f)};
  std::vector<float>      h = {0.f, 1.f};

  hmap::Array z0 = hmap::elevation_from_contours(shape, c, h, nullptr, 0.f, 0);
  hmap::Array z1 = hmap::elevation_from_contours(shape, c, h, nullptr, 0.f, 42);

  for (int k = 0; k < z0.size(); ++k)
    EXPECT_FLOAT_EQ(z0(k), z1(k)) << "k=" << k;
}

TEST(ElevationFromContours, StochasticModeIsSeedReproducible)
{
  glm::ivec2              shape = {48, 48};
  std::vector<hmap::Path> c = {square(0.5f, 0.5f, 0.4f),
                               square(0.5f, 0.5f, 0.2f)};
  std::vector<float>      h = {0.f, 1.f};

  hmap::Array a = hmap::elevation_from_contours(shape, c, h, nullptr, 1.f, 7);
  hmap::Array b = hmap::elevation_from_contours(shape, c, h, nullptr, 1.f, 7);
  hmap::Array d = hmap::elevation_from_contours(shape, c, h, nullptr, 1.f, 8);

  for (int k = 0; k < a.size(); ++k)
    EXPECT_FLOAT_EQ(a(k), b(k)) << "k=" << k;

  int ndiff = 0;
  for (int k = 0; k < a.size(); ++k)
    ndiff += (a(k) != d(k));
  EXPECT_GT(ndiff, 0);
}

TEST(ElevationFromContours, LeafZoneRisesAboveItsContour)
{
  glm::ivec2              shape = {64, 64};
  std::vector<hmap::Path> c = {square(0.5f, 0.5f, 0.4f),
                               square(0.5f, 0.5f, 0.2f)};
  std::vector<float>      h = {0.2f, 0.6f};

  hmap::Array z = hmap::elevation_from_contours(shape,
                                                c,
                                                h,
                                                nullptr,
                                                0.f,
                                                0,
                                                0.5f,
                                                1.f);

  // spacing = 0.4, peak = 0.6 + 0.5 * 0.4 = 0.8 at the centre (farthest point)
  EXPECT_NEAR(z(32, 32), 0.8f, 0.05f);
  EXPECT_GT(z(28, 32), 0.6f);
  EXPECT_LT(z(28, 32), 0.8f);
}

TEST(ElevationFromContours, OutsideDropsBelowRootContour)
{
  glm::ivec2              shape = {64, 64};
  std::vector<hmap::Path> c = {square(0.5f, 0.5f, 0.2f)};
  std::vector<float>      h = {0.5f};

  // single contour: spacing falls back to 1 -> outside floor 0.5 - 1.0 * 1
  hmap::Array z = hmap::elevation_from_contours(shape,
                                                c,
                                                h,
                                                nullptr,
                                                0.f,
                                                0,
                                                0.5f,
                                                1.f);

  EXPECT_NEAR(z(0, 0), -0.5f, 0.05f); // farthest corner
  EXPECT_LT(z(10, 32), 0.5f);         // outside, near the contour
  EXPECT_GT(z(10, 32), -0.5f);
  EXPECT_NEAR(z(32, 32), 1.0f, 0.05f); // leaf peak: 0.5 + 0.5 * 1
}

TEST(ElevationFromContours, BasinLeafSinksBelowItsContour)
{
  glm::ivec2              shape = {64, 64};
  std::vector<hmap::Path> c = {square(0.5f, 0.5f, 0.4f),
                               square(0.5f, 0.5f, 0.2f)};
  std::vector<float>      h = {0.6f, 0.2f};

  hmap::Array z = hmap::elevation_from_contours(shape,
                                                c,
                                                h,
                                                nullptr,
                                                0.f,
                                                0,
                                                0.5f,
                                                1.f);

  // inner contour lower than its parent: interior sinks to 0.2 - 0.5 * 0.4
  EXPECT_NEAR(z(32, 32), 0.0f, 0.05f);
  EXPECT_LT(z(28, 32), 0.2f);

  // ring between the contours descends from 0.6 to 0.2
  EXPECT_GT(z(12, 32), 0.2f);
  EXPECT_LT(z(12, 32), 0.6f);
}
