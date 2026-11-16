#include "highmap/internal/assert.hpp"
#include "highmap/morphology.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(RemoveEndpoints, SingleEndpoint)
{
  Array input{
      {0.f, 0.f, 0.f, 0.f, 0.f},
      {0.f, 1.f, 1.f, 1.f, 0.f},
      {0.f, 0.f, 0.f, 1.f, 0.f},
      {0.f, 0.f, 0.f, 0.f, 0.f},
  };

  Array expected{
      {0.f, 0.f, 0.f, 0.f, 0.f},
      {0.f, 0.f, 1.f, 1.f, 0.f},
      {0.f, 0.f, 0.f, 1.f, 0.f},
      {0.f, 0.f, 0.f, 0.f, 0.f},
  };

  Array output = remove_endpoints(input, 1, 0.f);

  bool ret = assert_almost_equal(output, expected, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(RemoveEndpoints, TwoIterations)
{
  Array input{
      {0.f, 0.f, 0.f, 0.f, 0.f},
      {0.f, 1.f, 1.f, 1.f, 1.f},
      {0.f, 0.f, 0.f, 0.f, 0.f},
  };

  Array expected{
      {0.f, 0.f, 0.f, 0.f, 0.f},
      {0.f, 0.f, 1.f, 0.f, 0.f},
      {0.f, 0.f, 0.f, 0.f, 0.f},
  };

  Array output = remove_endpoints(input, 2, 0.f);

  bool ret = assert_almost_equal(output, expected, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(RemoveEndpoints, RemoveLonelyPixel)
{
  Array input{
      {0.f, 0.f, 0.f},
      {0.f, 1.f, 0.f},
      {0.f, 0.f, 0.f},
  };

  Array expected{
      {0.f, 0.f, 0.f},
      {0.f, 0.f, 0.f},
      {0.f, 0.f, 0.f},
  };

  Array output = remove_endpoints(input, 1, 0.f);

  bool ret = assert_almost_equal(output, expected, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(RemoveEndpoints, KeepBranch)
{
  Array input{
      {0.f, 1.f, 0.f},
      {1.f, 1.f, 1.f},
      {0.f, 1.f, 0.f},
  };

  Array output = remove_endpoints(input, 1, 0.f);

  bool ret = assert_almost_equal(output, input, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(RemoveEndpoints, CustomBackground)
{
  Array input{
      {-1.f, -1.f, -1.f},
      {-1.f, 2.f, -1.f},
      {-1.f, -1.f, -1.f},
  };

  Array expected{
      {-1.f, -1.f, -1.f},
      {-1.f, -1.f, -1.f},
      {-1.f, -1.f, -1.f},
  };

  Array output = remove_endpoints(input, 1, -1.f);

  bool ret = assert_almost_equal(output, expected, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(RemoveEndpoints, ZeroIterations)
{
  Array input{
      {0.f, 1.f},
      {1.f, 0.f},
  };

  Array output = remove_endpoints(input, 0, 0.f);

  bool ret = assert_almost_equal(output, input, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(RemoveEndpoints, DiagonalConnectivity)
{
  Array input{
      {1.f, 0.f, 0.f},
      {0.f, 1.f, 0.f},
      {0.f, 0.f, 1.f},
  };

  Array expected{
      {0.f, 0.f, 0.f},
      {0.f, 1.f, 0.f},
      {0.f, 0.f, 0.f},
  };

  Array output = remove_endpoints(input, 1, 0.f);

  bool ret = assert_almost_equal(output, expected, 1e-6f);
  EXPECT_TRUE(ret);
}
