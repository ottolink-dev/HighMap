#include <map>

#include "highmap/internal/assert.hpp"
#include "highmap/morphology.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(AreaRemove, NoRemovalWhenBelowThreshold)
{
  // threshold smaller than component → nothing removed
  Array input = Array({{1, 1}, {1, 1}});
  Array out = area_remove(input, 2.f, 0.f, -1.f);

  EXPECT_TRUE(assert_almost_equal(out, input));
}

TEST(AreaRemove, RemoveSinglePixel)
{
  Array input = Array({{0, 1, 0}, {0, 0, 0}});
  Array expected = Array({{0, -1, 0}, {0, 0, 0}});
  Array out = area_remove(input, 2.f, 0.f, -1.f);

  EXPECT_TRUE(assert_almost_equal(out, expected));
}

TEST(AreaRemove, KeepLargeComponentRemoveSmallOnes)
{
  // Components:
  // - value 1 → size 2 (keep)
  // - value 2 → size 1 (remove)
  // - value 3 → size 2 (keep)

  Array input = Array({{1, 0, 2}, {1, 0, 0}, {0, 3, 3}});
  Array expected = Array({{1, 0, -1}, {1, 0, 0}, {0, 3, 3}});
  Array out = area_remove(input, 2.f, 0.f, -1.f);

  EXPECT_TRUE(assert_almost_equal(out, expected));
}

TEST(AreaRemove, BackgroundIsIgnored)
{
  // background = 5 → only center pixel is foreground
  Array input = Array({{5, 5, 5}, {5, 1, 5}, {5, 5, 5}});
  Array expected = Array({{5, 5, 5}, {5, -1, 5}, {5, 5, 5}});
  Array out = area_remove(input, 2.f, 5.f, -1.f);

  EXPECT_TRUE(assert_almost_equal(out, expected));
}

TEST(AreaRemove, MultipleDisconnectedComponents)
{
  // all components size 1 → all removed
  Array input = Array({{1, 0, 1}, {0, 0, 0}, {1, 0, 1}});
  Array expected = Array({{-1, 0, -1}, {0, 0, 0}, {-1, 0, -1}});
  Array out = area_remove(input, 2.f, 0.f, -1.f);

  EXPECT_TRUE(assert_almost_equal(out, expected));
}

TEST(AreaRemove, DiagonalConnectivityMatters)
{
  // with this connected_components → diagonal is connected → size = 2
  // so nothing should be removed (threshold is <, not <=)
  Array input = Array({{1, 0}, {0, 1}});
  Array out = area_remove(input, 2.f, 0.f, -1.f);

  EXPECT_TRUE(assert_almost_equal(out, input));
}

TEST(AreaRemove, EmptyAreaMapReturnsOriginal)
{
  Array input = Array({{0, 0}, {0, 0}});
  Array out = area_remove(input, 2.f, 0.f, -1.f);

  EXPECT_TRUE(assert_almost_equal(out, input));
}
