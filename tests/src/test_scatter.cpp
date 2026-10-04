#include <cmath>
#include <fstream>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/scatter.hpp"

#include <gtest/gtest.h>

using namespace hmap;

// --- Helper Functions

static constexpr float eps = 1e-5f;

static bool float_eq(float a, float b, float tol = eps)
{
  return std::abs(a - b) < tol;
}

// --- ScatterItem Unit Tests

TEST(ScatterItemTest, DefaultConstructor)
{
  ScatterItem item;
  EXPECT_TRUE(float_eq(item.position.x, 0.f));
  EXPECT_TRUE(float_eq(item.position.y, 0.f));
  EXPECT_TRUE(float_eq(item.position.z, 0.f));
  EXPECT_EQ(item.class_id, 0u);
  EXPECT_TRUE(float_eq(item.radius, HMAP_DEFAULT_SCATTER_RADIUS));
}

TEST(ScatterItemTest, ParameterizedConstructors)
{
  ScatterItem it1(glm::vec3(1.f, 2.f, 3.f), 2u, 4.5f);
  EXPECT_TRUE(float_eq(it1.position.x, 1.f));
  EXPECT_TRUE(float_eq(it1.position.y, 2.f));
  EXPECT_TRUE(float_eq(it1.position.z, 3.f));
  EXPECT_EQ(it1.class_id, 2u);
  EXPECT_TRUE(float_eq(it1.radius, 4.5f));

  ScatterItem it2(glm::vec2(5.f, 6.f), 1u, 2.0f);
  EXPECT_TRUE(float_eq(it2.position.x, 5.f));
  EXPECT_TRUE(float_eq(it2.position.y, 6.f));
  EXPECT_TRUE(float_eq(it2.position.z, 0.f));
  EXPECT_EQ(it2.class_id, 1u);

  ScatterItem it3(7.f, 8.f, 9.f, 3u, 1.5f);
  EXPECT_TRUE(float_eq(it3.position.x, 7.f));
  EXPECT_TRUE(float_eq(it3.position.y, 8.f));
  EXPECT_TRUE(float_eq(it3.position.z, 9.f));
  EXPECT_EQ(it3.class_id, 3u);
}

TEST(ScatterItemTest, ConversionsAndOperators)
{
  ScatterItem it(1.f, 2.f, 3.f, 1u, 2.5f);

  Point p = it.to_point();
  EXPECT_TRUE(float_eq(p.x, 1.f));
  EXPECT_TRUE(float_eq(p.y, 2.f));
  EXPECT_TRUE(float_eq(p.v, 2.5f));

  glm::vec2 v2 = it.to_vec2();
  EXPECT_TRUE(float_eq(v2.x, 1.f));
  EXPECT_TRUE(float_eq(v2.y, 2.f));

  glm::vec3 v3 = it.to_vec3();
  EXPECT_TRUE(float_eq(v3.x, 1.f));
  EXPECT_TRUE(float_eq(v3.y, 2.f));
  EXPECT_TRUE(float_eq(v3.z, 3.f));

  std::string s = it.to_string();
  EXPECT_FALSE(s.empty());

  ScatterItem same(1.f, 2.f, 3.f, 1u, 2.5f);
  ScatterItem diff(1.f, 2.f, 3.f, 2u, 2.5f);

  EXPECT_TRUE(it == same);
  EXPECT_FALSE(it == diff);
  EXPECT_TRUE(it != diff);
}

// --- ScatterField Unit Tests

TEST(ScatterFieldTest, ContainerOperations)
{
  ScatterField field;
  EXPECT_TRUE(field.empty());
  EXPECT_EQ(field.size(), 0u);

  field.push_back(ScatterItem(1.f, 2.f, 0.f, 1u, 3.f));
  field.emplace_back(4.f, 5.f, 0.f, 2u, 4.f);

  EXPECT_EQ(field.size(), 2u);
  EXPECT_FALSE(field.empty());
  EXPECT_TRUE(float_eq(field[0].position.x, 1.f));
  EXPECT_TRUE(float_eq(field.at(1).position.x, 4.f));
  EXPECT_TRUE(float_eq(field.front().position.y, 2.f));
  EXPECT_TRUE(float_eq(field.back().position.y, 5.f));

  field.clear();
  EXPECT_TRUE(field.empty());
}

TEST(ScatterFieldTest, FromAndToCloud)
{
  Cloud cloud;
  cloud.push_back(Point(10.f, 20.f, 3.5f));
  cloud.push_back(Point(30.f, 40.f, 0.0f));

  ScatterField field(cloud, 4u, 2.0f);
  ASSERT_EQ(field.size(), 2u);
  EXPECT_TRUE(float_eq(field[0].radius, 3.5f));
  EXPECT_TRUE(float_eq(field[1].radius, 2.0f));

  Cloud exported = field.to_cloud();
  ASSERT_EQ(exported.size(), 2u);
  EXPECT_TRUE(float_eq(exported[0].x, 10.f));
  EXPECT_TRUE(float_eq(exported[0].v, 3.5f));
}

TEST(ScatterFieldTest, FilterByClassAndBbox)
{
  ScatterField field;
  field.push_back(ScatterItem(1.f, 1.f, 0.f, 1u, 2.f));
  field.push_back(ScatterItem(2.f, 2.f, 0.f, 2u, 3.f));
  field.push_back(ScatterItem(3.f, 3.f, 0.f, 1u, 2.f));

  auto class_ids = field.get_class_ids();
  ASSERT_EQ(class_ids.size(), 2u);
  EXPECT_EQ(class_ids[0], 1u);
  EXPECT_EQ(class_ids[1], 2u);

  ScatterField filtered = field.filter_by_class(1u);
  EXPECT_EQ(filtered.size(), 2u);

  glm::vec4 bbox = field.get_bbox();
  EXPECT_TRUE(float_eq(bbox.x, 1.f));
  EXPECT_TRUE(float_eq(bbox.y, 3.f));
  EXPECT_TRUE(float_eq(bbox.z, 1.f));
  EXPECT_TRUE(float_eq(bbox.w, 3.f));
}

TEST(ScatterFieldTest, PruneCollisions)
{
  ScatterField field;
  // Two overlapping items (dist = 1.0, radii = 2.0 and 1.0)
  field.push_back(ScatterItem(0.0f, 0.0f, 0.0f, 0u, 2.0f));
  field.push_back(ScatterItem(1.0f, 0.0f, 0.0f, 0u, 1.0f));
  // One distant item
  field.push_back(ScatterItem(10.0f, 0.0f, 0.0f, 0u, 1.0f));

  field.prune_collisions();
  EXPECT_EQ(field.size(), 2u);
  EXPECT_TRUE(float_eq(field[0].radius, 2.0f));
  EXPECT_TRUE(float_eq(field[1].position.x, 10.0f));
}

TEST(ScatterFieldTest, SetElevationFromTerrain)
{
  Array elevation(glm::ivec2(3, 3), 0.f);
  for (int j = 0; j < 3; ++j)
    for (int i = 0; i < 3; ++i)
      elevation(i, j) = static_cast<float>(i) * 50.f;

  glm::vec4 bbox = {0.f, 10.f, 0.f, 10.f};

  ScatterField field;
  field.push_back(ScatterItem(0.f, 5.f, 0.f));
  field.push_back(ScatterItem(5.f, 5.f, 0.f));
  field.push_back(ScatterItem(10.f, 5.f, 0.f));

  field.set_elevation_from_terrain(elevation, bbox);

  EXPECT_TRUE(float_eq(field[0].position.z, 0.f));
  EXPECT_TRUE(float_eq(field[1].position.z, 50.f));
  EXPECT_TRUE(float_eq(field[2].position.z, 100.f));
}

TEST(ScatterFieldTest, Densify)
{
  ScatterField field;
  field.push_back(ScatterItem(0.0f, 0.0f, 0.0f, 1u, 0.1f));
  field.push_back(ScatterItem(1.0f, 0.0f, 0.0f, 1u, 0.1f));
  field.push_back(ScatterItem(0.0f, 1.0f, 0.0f, 2u, 0.1f));

  size_t original_size = field.size();
  field.densify(0.05f);
  EXPECT_GT(field.size(), original_size);
}

TEST(ScatterFieldTest, ReinforceClassClusters)
{
  // A center item surrounded by 4 neighbor items with class 1
  // Center item with class 2 should be converted to class 1
  ScatterField field;
  field.push_back(ScatterItem(0.0f, 0.0f, 0.0f, 2u, 0.1f)); // center
  field.push_back(ScatterItem(1.0f, 0.0f, 0.0f, 1u, 0.1f));
  field.push_back(ScatterItem(-1.0f, 0.0f, 0.0f, 1u, 0.1f));
  field.push_back(ScatterItem(0.0f, 1.0f, 0.0f, 1u, 0.1f));
  field.push_back(ScatterItem(0.0f, -1.0f, 0.0f, 1u, 0.1f));

  field.reinforce_class_clusters(1, 4, true);

  // Center item should now have class 1 (4 votes vs 1 vote)
  EXPECT_EQ(field[0].class_id, 1u);
}

TEST(ScatterFieldTest, ShuffleClasses)
{
  ScatterField field;
  field.push_back(ScatterItem(0.0f, 0.0f, 0.0f, 0u, 1.0f));
  field.push_back(ScatterItem(0.1f, 0.0f, 0.0f, 1u, 2.0f));
  field.push_back(ScatterItem(0.2f, 0.0f, 0.0f, 0u, 1.0f));
  field.push_back(ScatterItem(0.3f, 0.0f, 0.0f, 1u, 2.0f));

  field.shuffle_classes(1.0f, 2, 42);

  // Total class counts should remain invariant
  auto c0 = field.filter_by_class(0u);
  auto c1 = field.filter_by_class(1u);
  EXPECT_EQ(c0.size(), 2u);
  EXPECT_EQ(c1.size(), 2u);
}

TEST(ScatterSeedingTest, SeedScatterClustersAndKMeans)
{
  Array density({32, 32}, 1.0f);
  Array exclusion({32, 32}, 0.0f);

  ScatterSeedingOptions opts;
  opts.seed = 42;

  ScatterField clustered =
      seed_scatter_clusters(2, 50, density, exclusion, 0.05f, 5, opts);
  EXPECT_EQ(clustered.size(), 50u);

  ScatterField kmeans =
      seed_scatter_kmeans(2, 50, density, exclusion, 0.0f, 4, opts);
  EXPECT_EQ(kmeans.size(), 50u);
}
