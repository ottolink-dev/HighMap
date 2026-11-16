#include <map>

#include "highmap/features.hpp"
#include "highmap/internal/assert.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(ConnectedComponents, SingleComponent)
{
  Array input = Array({{1, 1}, {1, 1}});

  std::map<float, float>                surfaces;
  std::map<float, std::array<float, 2>> centroids;

  Array labels = connected_components(input, 0.f, 0.f, &surfaces, &centroids);

  // Expect only one component
  EXPECT_EQ(surfaces.size(), 1);

  // Surface should be 4 pixels
  EXPECT_FLOAT_EQ(surfaces.begin()->second, 4.f);

  // All labels should be identical and > 0
  float lbl = labels(0, 0);
  EXPECT_GT(lbl, 0);

  for (int i = 0; i < labels.shape.x; ++i)
    for (int j = 0; j < labels.shape.y; ++j)
      EXPECT_EQ(labels(i, j), lbl);
}

TEST(ConnectedComponents, TwoSeparateComponents)
{
  Array input = Array({{1, 0, 1}, {1, 0, 1}});

  std::map<float, float> surfaces;

  Array labels = connected_components(input, 0.f, 0.f, &surfaces, nullptr);

  // Expect 2 components
  EXPECT_EQ(surfaces.size(), 2);

  // Each should have surface 2
  for (auto &[k, v] : surfaces)
    EXPECT_FLOAT_EQ(v, 2.f);
}

TEST(ConnectedComponents, DiagonalConnectivity)
{
  Array input = Array({{1, 0}, {0, 1}});

  std::map<float, float> surfaces;

  Array labels = connected_components(input, 0.f, 0.f, &surfaces, nullptr);

  // Because of your neighbor pattern (diagonal included),
  // this SHOULD be one component
  EXPECT_EQ(surfaces.size(), 1);
}

TEST(ConnectedComponents, SurfaceThresholdFiltering)
{
  Array input = Array({{1, 0, 1}, {0, 0, 0}, {1, 1, 0}});

  std::map<float, float> surfaces;

  // threshold = 2 → remove single pixels
  Array labels = connected_components(input, 2.f, 0.f, &surfaces, nullptr);

  // Only bottom component (size 2) should remain
  EXPECT_EQ(surfaces.size(), 1);
  EXPECT_FLOAT_EQ(surfaces.begin()->second, 2.f);

  EXPECT_EQ(count_non_zero(labels), 2);
}

TEST(ConnectedComponents, LabelRemappingIsCompact)
{
  Array input = Array({{1, 0, 1}, {0, 0, 0}, {1, 0, 1}});
  Array labels = connected_components(input, 0.f, 0.f, nullptr, nullptr);

  auto values = labels.unique_values();

  // Labels should be contiguous: 0,1,2,3,...
  for (size_t i = 0; i < values.size(); ++i)
    EXPECT_EQ(values[i], float(i));
}

TEST(ConnectedComponents, CentroidComputation)
{
  Array input = Array({{1, 1}, {0, 0}});

  std::map<float, float>                surfaces;
  std::map<float, std::array<float, 2>> centroids;

  connected_components(input, 0.f, 0.f, &surfaces, &centroids);

  ASSERT_EQ(centroids.size(), 1);

  auto c = centroids.begin()->second;

  // centroid sum
  EXPECT_FLOAT_EQ(c[0], 0.f + 0.f); // i indices
  EXPECT_FLOAT_EQ(c[1], 0.f + 1.f); // j indices
}
