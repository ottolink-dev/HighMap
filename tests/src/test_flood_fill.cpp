#include "highmap/internal/assert.hpp"
#include "highmap/morphology.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(FloodFill, SimpleRegion)
{
  Array input{
      {0.f, 0.f, 1.f},
      {0.f, 0.f, 1.f},
      {1.f, 1.f, 1.f},
  };

  Array expected{
      {9.f, 9.f, 1.f},
      {9.f, 9.f, 1.f},
      {1.f, 1.f, 1.f},
  };

  flood_fill(input, 0, 0, 9.f, 0.f);

  bool ret = assert_almost_equal(input, expected, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(FloodFill, SingleCell)
{
  Array input{
      {0.f, 1.f},
      {1.f, 1.f},
  };

  Array expected{
      {7.f, 1.f},
      {1.f, 1.f},
  };

  flood_fill(input, 0, 0, 7.f, 0.f);

  bool ret = assert_almost_equal(input, expected, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(FloodFill, NoOpOnNonBackground)
{
  Array input{
      {2.f, 2.f},
      {2.f, 2.f},
  };

  Array expected = input;

  flood_fill(input, 0, 0, 9.f, 0.f);

  bool ret = assert_almost_equal(input, expected, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(FloodFill, DisconnectedRegions)
{
  Array input{
      {0.f, 1.f, 0.f},
      {1.f, 1.f, 1.f},
      {0.f, 1.f, 0.f},
  };

  Array expected{
      {8.f, 1.f, 0.f},
      {1.f, 1.f, 1.f},
      {0.f, 1.f, 0.f},
  };

  flood_fill(input, 0, 0, 8.f, 0.f);

  bool ret = assert_almost_equal(input, expected, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(FloodFill, BoundaryClipping)
{
  Array input{
      {0.f, 0.f},
      {0.f, 0.f},
  };

  Array expected{
      {5.f, 5.f},
      {5.f, 5.f},
  };

  flood_fill(input, 1, 1, 5.f, 0.f);

  bool ret = assert_almost_equal(input, expected, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(FloodFill, DoesNotOverwriteFillValue)
{
  Array input{
      {0.f, 0.f},
      {0.f, 9.f},
  };

  Array expected{
      {3.f, 3.f},
      {3.f, 9.f},
  };

  flood_fill(input, 0, 0, 3.f, 0.f);

  bool ret = assert_almost_equal(input, expected, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(FloodFill, AlreadyFilledSeed)
{
  Array input{
      {1.f, 1.f},
      {1.f, 1.f},
  };

  Array expected = input;

  flood_fill(input, 0, 0, 7.f, 0.f);

  bool ret = assert_almost_equal(input, expected, 1e-6f);
  EXPECT_TRUE(ret);
}
