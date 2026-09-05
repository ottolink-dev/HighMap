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
