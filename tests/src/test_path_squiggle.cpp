#include "highmap.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(PathSquiggle, BasicEdgeSubdivision)
{
  Path input_path({Point(0.f, 0.f), Point(10.f, 0.f)});

  Path result = squiggle(input_path, 4, 42);

  EXPECT_GT(result.size(), input_path.size());
  EXPECT_NEAR(result.front().x, 0.f, 1e-5f);
  EXPECT_NEAR(result.front().y, 0.f, 1e-5f);
  EXPECT_NEAR(result.back().x, 10.f, 1e-5f);
  EXPECT_NEAR(result.back().y, 0.f, 1e-5f);
}

TEST(PathSquiggle, FirstAndLastEdgeSubdivision)
{
  Path input_path({Point(0.f, 0.f), Point(10.f, 0.f)});

  Path p1 = squiggle(input_path, 1, 42);
  Path p2 = squiggle(input_path, 2, 42);
  Path p3 = squiggle(input_path, 3, 42);

  EXPECT_GT(p2.size(), p1.size());
  EXPECT_GT(p3.size(), p2.size());

  // Check that the first segment length decreases with each iteration
  float first_seg_len1 = distance(p1[0], p1[1]);
  float first_seg_len2 = distance(p2[0], p2[1]);
  float first_seg_len3 = distance(p3[0], p3[1]);

  EXPECT_LT(first_seg_len2, first_seg_len1);
  EXPECT_LT(first_seg_len3, first_seg_len2);

  // Check that the last segment length decreases with each iteration
  float last_seg_len1 = distance(p1[p1.size() - 2], p1[p1.size() - 1]);
  float last_seg_len2 = distance(p2[p2.size() - 2], p2[p2.size() - 1]);
  float last_seg_len3 = distance(p3[p3.size() - 2], p3[p3.size() - 1]);

  EXPECT_LT(last_seg_len2, last_seg_len1);
  EXPECT_LT(last_seg_len3, last_seg_len2);
}

TEST(PathSquiggle, ClosedPathPreservesClosure)
{
  Path input_path(
      {Point(0.f, 0.f), Point(10.f, 0.f), Point(10.f, 10.f), Point(0.f, 10.f)});
  input_path.set_closed(true);

  Path result = squiggle(input_path, 3, 123);

  EXPECT_TRUE(result.is_closed());
  EXPECT_GT(result.size(), input_path.size());
}

TEST(PathSquiggle, Determinism)
{
  Path input_path({Point(0.f, 0.f), Point(5.f, 5.f), Point(10.f, 0.f)});

  Path p1 = squiggle(input_path, 4, 999);
  Path p2 = squiggle(input_path, 4, 999);
  Path p3 = squiggle(input_path, 4, 111);

  EXPECT_EQ(p1.size(), p2.size());
  for (size_t i = 0; i < p1.size(); ++i)
  {
    EXPECT_FLOAT_EQ(p1[i].x, p2[i].x);
    EXPECT_FLOAT_EQ(p1[i].y, p2[i].y);
  }

  bool differs = (p1.size() != p3.size());
  if (!differs)
  {
    for (size_t i = 0; i < p1.size(); ++i)
    {
      if (std::abs(p1[i].x - p3[i].x) > 1e-4f ||
          std::abs(p1[i].y - p3[i].y) > 1e-4f)
      {
        differs = true;
        break;
      }
    }
  }
  EXPECT_TRUE(differs);
}

TEST(PathSquiggle, SpatialWeights)
{
  Path input_path({Point(0.f, 0.f), Point(10.f, 0.f)});

  glm::ivec2 shape = {64, 64};
  Array      weights(shape, 1.f);
  // Heavy weights in top half y > 0
  for (int j = shape.y / 2; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      weights(i, j) = 100.f;

  glm::vec4 bbox = {0.f, 10.f, -5.f, 5.f};
  Path      result = squiggle(input_path,
                         4,
                         42,
                         0.5f,
                         1, // positive normal (y > 0)
                         &weights,
                         nullptr,
                         bbox);

  EXPECT_GT(result.size(), 2u);
}

TEST(PathSquiggle, MaskAvoidance)
{
  Path input_path({Point(0.f, 0.f), Point(10.f, 0.f)});

  glm::ivec2 shape = {64, 64};
  Array      mask(shape, 1.f);
  // Forbid upper band
  for (int j = shape.y / 2; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      mask(i, j) = 0.f;

  glm::vec4 bbox = {0.f, 10.f, -5.f, 5.f};
  Path      result = squiggle(input_path, 4, 42, 0.5f, 1, nullptr, &mask, bbox);

  EXPECT_GT(result.size(), 2u);
}

TEST(PathSquiggle, Branching)
{
  Path input_path({Point(0.f, 0.f), Point(10.f, 0.f)});

  std::vector<Path> paths = squiggle_branches(input_path, 4, 42, 0.8f, 4);

  EXPECT_GE(paths.size(), 1u);
  EXPECT_GT(paths[0].size(), 2u);
}
