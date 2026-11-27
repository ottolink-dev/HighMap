#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/assert.hpp"
#include "highmap/primitives.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(MatchHistogram, Identity)
{
  Array input{
      {0.f, 1.f},
      {2.f, 3.f},
  };

  Array reference = input;
  Array expected = input;

  match_histogram(input, reference);

  bool ret = assert_almost_equal(input, expected, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(MatchHistogram, ConstantReference)
{
  Array input{
      {0.f, 1.f},
      {2.f, 3.f},
  };

  Array reference{
      {5.f, 5.f},
      {5.f, 5.f},
  };

  Array expected{
      {5.f, 5.f},
      {5.f, 5.f},
  };

  match_histogram(input, reference);

  bool ret = assert_almost_equal(input, expected, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(MatchHistogram, SortedRanks)
{
  Array input{
      {0.f, 1.f},
      {2.f, 3.f},
  };

  Array reference{
      {10.f, 20.f},
      {30.f, 40.f},
  };

  Array expected{
      {10.f, 20.f},
      {30.f, 40.f},
  };

  match_histogram(input, reference);

  bool ret = assert_almost_equal(input, expected, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(MatchHistogram, ReverseOrdering)
{
  Array input{
      {3.f, 2.f},
      {1.f, 0.f},
  };

  Array reference{
      {10.f, 20.f},
      {30.f, 40.f},
  };

  Array expected{
      {40.f, 30.f},
      {20.f, 10.f},
  };

  match_histogram(input, reference);

  bool ret = assert_almost_equal(input, expected, 1e-6f);
  EXPECT_TRUE(ret);
}

TEST(MatchHistogram, Interpolation)
{
  Array input{
      {0.f, 1.f, 2.f},
  };

  Array reference{
      {0.f, 10.f},
  };

  Array expected{
      {0.f, 5.f, 10.f},
  };

  match_histogram(input, reference);

  bool ret = assert_almost_equal(input, expected, 1e-5f);
  EXPECT_TRUE(ret);
}

TEST(MatchHistogram, RandomMonotonicity)
{
  Array input = white({32, 32}, 0.f, 1.f, 42);
  Array reference = white({32, 32}, 10.f, 20.f, 123);

  match_histogram(input, reference);

  float vmin = input.min();
  float vmax = input.max();

  EXPECT_GE(vmin, 10.f - 1e-4f);
  EXPECT_LE(vmax, 20.f + 1e-4f);
}

TEST(MatchHistogram, SingleValueInput)
{
  Array input{
      {42.f},
  };

  Array reference{
      {7.f},
  };

  Array expected{
      {7.f},
  };

  match_histogram(input, reference);

  ASSERT_EQ(input.vector.size(), 1u);
  EXPECT_FLOAT_EQ(input.vector.front(), expected.vector.front());
}

TEST(MatchHistogram, SingleInputUsesOddReferenceMedian)
{
  Array input{{42.f}};
  Array reference{{9.f, 1.f, 5.f}};
  match_histogram(input, reference);
  EXPECT_FLOAT_EQ(input.vector.front(), 5.f);
}

TEST(MatchHistogram, SingleInputInterpolatesEvenReferenceMedian)
{
  Array input{{42.f}};
  Array reference{{30.f, 0.f, 20.f, 10.f}};
  match_histogram(input, reference);
  EXPECT_FLOAT_EQ(input.vector.front(), 15.f);
}

TEST(MatchHistogram, SingleReferenceMapsEveryInput)
{
  Array input{{4.f, -2.f, 9.f}};
  Array reference{{7.f}};
  match_histogram(input, reference);
  for (float value : input.vector)
    EXPECT_FLOAT_EQ(value, 7.f);
}

TEST(MatchHistogram, EmptyInputAndReferenceAreNoOps)
{
  Array empty;
  Array reference{{7.f}};
  match_histogram(empty, reference);
  EXPECT_TRUE(empty.vector.empty());

  Array input{{4.f, -2.f, 9.f}};
  const auto before = input.vector;
  match_histogram(input, empty);
  EXPECT_EQ(input.vector, before);
}
