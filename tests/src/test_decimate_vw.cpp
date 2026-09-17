#include "highmap.hpp"

#include <gtest/gtest.h>

using namespace hmap;

std::vector<Point> generate_line(int n)
{
  std::vector<Point> pts;
  for (int i = 0; i < n; ++i)
    pts.push_back({float(i), 0.0f});
  return pts;
}

std::vector<Point> generate_circle(int n, float r)
{
  std::vector<Point> pts;
  for (int i = 0; i < n; ++i)
  {
    float t = 2.f * M_PI * i / n;
    pts.push_back({r * std::cos(t), r * std::sin(t)});
  }
  return pts;
}

TEST(PathDecimateVw, ReducesPointCount)
{
  Path path(generate_line(100));
  path = decimate_vw(path, 10);

  EXPECT_EQ(path.size(), 10);
}

TEST(PathDecimateVw, PreservesStraightLine)
{
  Path path(generate_line(100));
  path = decimate_vw(path, 5);

  // all points should remain colinear
  for (size_t i = 1; i + 1 < path.size(); ++i)
  {
    float k = curvature(path[i - 1], path[i], path[i + 1]);
    EXPECT_NEAR(k, 0.0f, 1e-5f);
  }
}

TEST(PathDecimateVw, PreservesCircle)
{
  Path path(generate_circle(100, 1.f));
  path = decimate_vw(path, 20);

  float avg_radius = 1.f;

  for (auto &p : path)
  {
    float r = std::hypot(p.x, p.y);
    EXPECT_NEAR(r, avg_radius, 0.05f);
  }
}
