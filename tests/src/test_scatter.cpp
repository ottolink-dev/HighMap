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

TEST(ScatterFieldTest, PruneCollisionsAgainstOtherField)
{
  // field1 (e.g. rocks) - unmodified reference
  ScatterField rocks;
  rocks.push_back(ScatterItem(0.0f, 0.0f, 0.0f, 0u, 2.0f));
  rocks.push_back(ScatterItem(10.0f, 10.0f, 0.0f, 0u, 1.5f));

  // field2 (e.g. forest trees)
  ScatterField trees;
  trees.push_back(
      ScatterItem(1.0f,
                  0.0f,
                  0.0f,
                  1u,
                  1.0f)); // inside rock radius (dist = 1.0 < 2.0 + 1.0)
  trees.push_back(ScatterItem(5.0f,
                              0.0f,
                              0.0f,
                              1u,
                              1.0f)); // far from rocks (dist = 5.0 > 2.0 + 1.0)
  trees.push_back(
      ScatterItem(10.5f,
                  10.0f,
                  0.0f,
                  1u,
                  1.0f)); // inside second rock radius (dist = 0.5 < 1.5 + 1.0)

  trees.prune_collisions(rocks);

  // Only the middle tree (at x=5.0) should remain
  EXPECT_EQ(trees.size(), 1u);
  EXPECT_TRUE(float_eq(trees[0].position.x, 5.0f));

  // Reference rocks field remains unchanged
  EXPECT_EQ(rocks.size(), 2u);
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

TEST(ScatterFieldTest, PerturbPositions)
{
  ScatterField field;
  field.push_back(ScatterItem(0.5f, 0.5f, 0.0f, 0u, 1.0f));
  field.push_back(ScatterItem(1.0f, 1.0f, 0.0f, 1u, 1.0f));

  field.perturb_positions(0.1f, 0.1f, 42);

  // Positions should have moved within +/- dx, dy
  EXPECT_FALSE(float_eq(field[0].position.x, 0.5f));
  EXPECT_FALSE(float_eq(field[0].position.y, 0.5f));
  EXPECT_NEAR(field[0].position.x, 0.5f, 0.1001f);
  EXPECT_NEAR(field[0].position.y, 0.5f, 0.1001f);

  EXPECT_FALSE(float_eq(field[1].position.x, 1.0f));
  EXPECT_FALSE(float_eq(field[1].position.y, 1.0f));
  EXPECT_NEAR(field[1].position.x, 1.0f, 0.1001f);
  EXPECT_NEAR(field[1].position.y, 1.0f, 0.1001f);
}

TEST(ScatterFieldTest, RegularizePositions)
{
  ScatterField field;
  // Two points placed very close to each other
  field.push_back(ScatterItem(0.5f, 0.5f, 0.0f, 0u, 1.0f));
  field.push_back(ScatterItem(0.501f, 0.5f, 0.0f, 0u, 1.0f));

  float initial_dist = std::abs(field[1].position.x - field[0].position.x);

  field.regularize_positions(1, 0.05f, 5);

  float new_dist = std::abs(field[1].position.x - field[0].position.x);
  EXPECT_GT(new_dist, initial_dist);
}

TEST(ScatterFieldTest, ResolveCollisionsSelf)
{
  ScatterField field;
  // 3 overlapping items with radii 1.0 each (so required min dist is 2.0)
  field.push_back(ScatterItem(0.0f, 0.0f, 0.0f, 0u, 1.0f));
  field.push_back(ScatterItem(0.5f, 0.0f, 0.0f, 0u, 1.0f));
  field.push_back(ScatterItem(0.25f, 0.4f, 0.0f, 0u, 1.0f));

  size_t initial_count = field.size();

  // Resolve collisions over 30 iterations with triangulation recomputed every 5
  // steps
  field.resolve_collisions(40, 0.0f, 0.5f, 5);

  // Object count must remain unchanged
  EXPECT_EQ(field.size(), initial_count);

  // Check distances between all pairs are >= required sum of radii (within
  // small numerical margin)
  for (size_t i = 0; i < field.size(); ++i)
  {
    for (size_t j = i + 1; j < field.size(); ++j)
    {
      float dx = field[i].position.x - field[j].position.x;
      float dy = field[i].position.y - field[j].position.y;
      float dist = std::sqrt(dx * dx + dy * dy);
      float min_dist = field[i].radius + field[j].radius;
      EXPECT_GE(dist, min_dist - 1e-2f);
    }
  }
}

TEST(ScatterFieldTest, ResolveCollisionsAgainstOther)
{
  ScatterField obstacles;
  // Fixed obstacle at (0, 0) with radius 2.0
  obstacles.push_back(ScatterItem(0.0f, 0.0f, 0.0f, 0u, 2.0f));

  ScatterField moving;
  // Item overlapping the obstacle at (0.5, 0.0) with radius 1.0 -> required
  // dist = 3.0
  moving.push_back(ScatterItem(0.5f, 0.0f, 0.0f, 1u, 1.0f));

  moving.resolve_collisions(obstacles, 40, 0.0f, 0.5f);

  EXPECT_EQ(moving.size(), 1u);
  EXPECT_EQ(obstacles.size(), 1u);

  // Obstacle should remain unmoved at (0, 0)
  EXPECT_TRUE(float_eq(obstacles[0].position.x, 0.0f));
  EXPECT_TRUE(float_eq(obstacles[0].position.y, 0.0f));

  // Moving item should now be pushed to distance >= 3.0
  float dx = moving[0].position.x - obstacles[0].position.x;
  float dy = moving[0].position.y - obstacles[0].position.y;
  float dist = std::sqrt(dx * dx + dy * dy);
  EXPECT_GE(dist, 3.0f - 1e-2f);
}

TEST(ScatterFieldTest, ResolveCollisionsLargeRandom)
{
  Array density({64, 64}, 1.0f);
  Array exclusion({64, 64}, 0.0f);

  ScatterSeedingOptions opts;
  opts.seed = 12345;
  ScatterField field =
      seed_scatter_clusters(1, 100, density, exclusion, 0.05f, 1, opts);

  // Set uniform radius
  for (size_t i = 0; i < field.size(); ++i)
  {
    field[i].radius = 0.03f;
  }

  auto count_collisions = [](const ScatterField &sf)
  {
    size_t collisions = 0;
    for (size_t i = 0; i < sf.size(); ++i)
    {
      for (size_t j = i + 1; j < sf.size(); ++j)
      {
        float dx = sf[i].position.x - sf[j].position.x;
        float dy = sf[i].position.y - sf[j].position.y;
        float dist_sq = dx * dx + dy * dy;
        float min_dist = sf[i].radius + sf[j].radius;
        if (dist_sq < min_dist * min_dist)
        {
          collisions++;
        }
      }
    }
    return collisions;
  };

  size_t before = count_collisions(field);

  field.resolve_collisions(50, 0.0f, 0.5f, 1);

  size_t after = count_collisions(field);

  EXPECT_LE(after, before);
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

// --- Merge Functions Unit Tests

TEST(ScatterFieldTest, MergeScatterFieldsPreserveClasses)
{
  ScatterField f1;
  f1.push_back(ScatterItem(1.f, 2.f, 3.f, 10u, 0.5f));
  f1.push_back(ScatterItem(4.f, 5.f, 6.f, 20u, 1.0f));

  ScatterField f2;
  f2.push_back(ScatterItem(7.f, 8.f, 9.f, 10u, 1.5f));
  f2.push_back(ScatterItem(10.f, 11.f, 12.f, 30u, 2.0f));

  ScatterField merged = merge_scatter_fields({f1, f2}, true);

  EXPECT_EQ(merged.size(), 4u);
  EXPECT_EQ(merged[0].class_id, 10u);
  EXPECT_TRUE(float_eq(merged[0].position.x, 1.f));
  EXPECT_TRUE(float_eq(merged[0].radius, 0.5f));

  EXPECT_EQ(merged[1].class_id, 20u);
  EXPECT_EQ(merged[2].class_id, 10u);
  EXPECT_TRUE(float_eq(merged[2].position.x, 7.f));
  EXPECT_TRUE(float_eq(merged[2].radius, 1.5f));

  EXPECT_EQ(merged[3].class_id, 30u);
}

TEST(ScatterFieldTest, MergeScatterFieldsRemapClasses)
{
  ScatterField f1;
  f1.push_back(ScatterItem(1.f, 2.f, 3.f, 10u, 0.5f));
  f1.push_back(ScatterItem(4.f, 5.f, 6.f, 20u, 1.0f));

  ScatterField f2;
  f2.push_back(ScatterItem(7.f, 8.f, 9.f, 10u, 1.5f));
  f2.push_back(ScatterItem(10.f, 11.f, 12.f, 30u, 2.0f));

  ScatterField merged = merge_scatter_fields({f1, f2}, false);

  EXPECT_EQ(merged.size(), 4u);
  // f1 items should have class 0
  EXPECT_EQ(merged[0].class_id, 0u);
  EXPECT_EQ(merged[1].class_id, 0u);
  // f2 items should have class 1
  EXPECT_EQ(merged[2].class_id, 1u);
  EXPECT_EQ(merged[3].class_id, 1u);

  // Check that other properties were preserved
  EXPECT_TRUE(float_eq(merged[0].position.x, 1.f));
  EXPECT_TRUE(float_eq(merged[0].position.y, 2.f));
  EXPECT_TRUE(float_eq(merged[0].position.z, 3.f));
  EXPECT_TRUE(float_eq(merged[0].radius, 0.5f));

  EXPECT_TRUE(float_eq(merged[2].position.x, 7.f));
  EXPECT_TRUE(float_eq(merged[2].position.y, 8.f));
  EXPECT_TRUE(float_eq(merged[2].position.z, 9.f));
  EXPECT_TRUE(float_eq(merged[2].radius, 1.5f));
}

TEST(ScatterFieldTest, MergeScatterFieldsMultiple)
{
  ScatterField f1;
  f1.push_back(ScatterItem(1.f, 1.f, 0.f, 5u, 0.1f));

  ScatterField f2;
  f2.push_back(ScatterItem(2.f, 2.f, 0.f, 6u, 0.2f));

  ScatterField f3;
  f3.push_back(ScatterItem(3.f, 3.f, 0.f, 7u, 0.3f));

  // Merge empty vector
  EXPECT_EQ(merge_scatter_fields({}).size(), 0u);

  // Merge with merge_by_class = true
  ScatterField merged_by_class = merge_scatter_fields({f1, f2, f3}, true);
  EXPECT_EQ(merged_by_class.size(), 3u);
  EXPECT_EQ(merged_by_class[0].class_id, 5u);
  EXPECT_EQ(merged_by_class[1].class_id, 6u);
  EXPECT_EQ(merged_by_class[2].class_id, 7u);

  // Merge with merge_by_class = false (new class per input field)
  ScatterField merged_remapped = merge_scatter_fields({f1, f2, f3}, false);
  EXPECT_EQ(merged_remapped.size(), 3u);
  EXPECT_EQ(merged_remapped[0].class_id, 0u);
  EXPECT_EQ(merged_remapped[1].class_id, 1u);
  EXPECT_EQ(merged_remapped[2].class_id, 2u);
  EXPECT_TRUE(float_eq(merged_remapped[2].position.x, 3.f));
  EXPECT_TRUE(float_eq(merged_remapped[2].radius, 0.3f));
}

TEST(ScatterFieldTest, MergeScatterFieldsWithBbox)
{
  // Field 1 covers [0, 10] x [0, 10]
  ScatterField f1;
  f1.push_back(ScatterItem(2.f,
                           2.f,
                           0.f,
                           10u,
                           0.5f)); // inside non-overlapping region of f1
  f1.push_back(
      ScatterItem(8.f, 5.f, 0.f, 10u, 0.5f)); // inside overlap region with f2
  glm::vec4 bbox1 = {0.f, 10.f, 0.f, 10.f};

  // Field 2 covers [6, 16] x [0, 10]
  ScatterField f2;
  f2.push_back(ScatterItem(8.f,
                           5.f,
                           0.f,
                           20u,
                           0.8f)); // in overlap -> should be discarded
  f2.push_back(ScatterItem(14.f,
                           5.f,
                           0.f,
                           20u,
                           0.8f)); // inside non-overlapping region of f2
  f2.push_back(ScatterItem(20.f,
                           5.f,
                           0.f,
                           20u,
                           0.8f)); // outside bbox2 -> should be discarded
  glm::vec4 bbox2 = {6.f, 16.f, 0.f, 10.f};

  // Field 3 covers [12, 22] x [0, 10]
  ScatterField f3;
  f3.push_back(ScatterItem(14.f,
                           5.f,
                           0.f,
                           30u,
                           1.0f)); // in overlap with f2 -> should be discarded
  f3.push_back(ScatterItem(18.f,
                           5.f,
                           0.f,
                           30u,
                           1.0f)); // inside non-overlapping region of f3
  glm::vec4 bbox3 = {12.f, 22.f, 0.f, 10.f};

  // Merge with merge_by_class = true
  ScatterField merged = merge_scatter_fields({f1, f2, f3},
                                             {bbox1, bbox2, bbox3},
                                             true);

  // Expected items:
  // From f1: (2, 2) and (8, 5)
  // From f2: (14, 5)
  // From f3: (18, 5)
  EXPECT_EQ(merged.size(), 4u);

  EXPECT_TRUE(float_eq(merged[0].position.x, 2.f));
  EXPECT_EQ(merged[0].class_id, 10u);

  EXPECT_TRUE(float_eq(merged[1].position.x, 8.f));
  EXPECT_EQ(merged[1].class_id, 10u);

  EXPECT_TRUE(float_eq(merged[2].position.x, 14.f));
  EXPECT_EQ(merged[2].class_id, 20u);

  EXPECT_TRUE(float_eq(merged[3].position.x, 18.f));
  EXPECT_EQ(merged[3].class_id, 30u);

  // Merge with merge_by_class = false (remapped classes)
  ScatterField merged_remap = merge_scatter_fields({f1, f2, f3},
                                                   {bbox1, bbox2, bbox3},
                                                   false);
  EXPECT_EQ(merged_remap.size(), 4u);
  EXPECT_EQ(merged_remap[0].class_id, 0u);
  EXPECT_EQ(merged_remap[1].class_id, 0u);
  EXPECT_EQ(merged_remap[2].class_id, 1u);
  EXPECT_EQ(merged_remap[3].class_id, 2u);

  // Mismatched size returns empty
  ScatterField merged_mismatch = merge_scatter_fields({f1, f2}, {bbox1});
  EXPECT_EQ(merged_mismatch.size(), 0u);
}

TEST(ScatterFieldTest, ToHeightmapEmpty)
{
  ScatterField empty_field;
  Array        hmap = empty_field.to_heightmap({64, 64});
  EXPECT_EQ(hmap.shape.x, 64);
  EXPECT_EQ(hmap.shape.y, 64);
  EXPECT_FLOAT_EQ(hmap.max(), 0.0f);
  EXPECT_FLOAT_EQ(hmap.min(), 0.0f);
}

TEST(ScatterFieldTest, ToHeightmapDisk)
{
  ScatterField field;
  field.push_back(ScatterItem(0.5f, 0.5f, 0.0f, 0u, 0.2f));

  glm::ivec2 shape = {128, 128};
  float      ratio = 1.5f;
  Array      hmap = field.to_heightmap(shape, SCATTER_SHAPE_DISK, ratio);

  EXPECT_EQ(hmap.shape.x, shape.x);
  EXPECT_EQ(hmap.shape.y, shape.y);

  // Peak should be at/near center with elevation ~ radius * ratio = 0.2 * 1.5 =
  // 0.3
  float max_val = hmap.max();
  EXPECT_NEAR(max_val, 0.3f, 0.02f);

  // Corners should have 0 elevation
  EXPECT_FLOAT_EQ(hmap(0, 0), 0.0f);
  EXPECT_FLOAT_EQ(hmap(shape.x - 1, shape.y - 1), 0.0f);
}

TEST(ScatterFieldTest, ToHeightmapPolygonAndShapes)
{
  ScatterField field;
  field.push_back(ScatterItem(0.5f, 0.5f, 0.0f, 0u, 0.2f));

  glm::ivec2 shape = {64, 64};
  Array      hmap_poly = field.to_heightmap(shape, SCATTER_SHAPE_POLYGON, 1.0f);
  EXPECT_GT(hmap_poly.max(), 0.15f);
  EXPECT_FLOAT_EQ(hmap_poly(0, 0), 0.0f);

  Array hmap_cone = field.to_heightmap(shape, SCATTER_SHAPE_CONE, 1.0f);
  EXPECT_GT(hmap_cone.max(), 0.15f);

  Array hmap_pyr = field.to_heightmap(shape, SCATTER_SHAPE_PYRAMID, 1.0f);
  EXPECT_GT(hmap_pyr.max(), 0.15f);

  Array hmap_dome = field.to_heightmap(shape, SCATTER_SHAPE_SMOOTH_DOME, 1.0f);
  EXPECT_GT(hmap_dome.max(), 0.15f);
}

TEST(ScatterFieldTest, ToHeightmapClassFilterAndFreeFunction)
{
  ScatterField field;
  field.push_back(ScatterItem(0.3f, 0.3f, 0.0f, 0u, 0.1f));
  field.push_back(ScatterItem(0.7f, 0.7f, 0.0f, 1u, 0.1f));

  glm::ivec2 shape = {64, 64};

  // Filter only class 1
  Array hmap_class1 = field.to_heightmap(shape,
                                         SCATTER_SHAPE_DISK,
                                         1.0f,
                                         std::make_optional<uint32_t>(1u));

  // Region around (0.3, 0.3) should be 0, while region around (0.7, 0.7) should
  // be > 0
  int ix0 = static_cast<int>(0.3f * shape.x);
  int iy0 = static_cast<int>(0.3f * shape.y);
  int ix1 = static_cast<int>(0.7f * shape.x);
  int iy1 = static_cast<int>(0.7f * shape.y);

  EXPECT_FLOAT_EQ(hmap_class1(ix0, iy0), 0.0f);
  EXPECT_GT(hmap_class1(ix1, iy1), 0.05f);

  // Test free function wrapper
  Array hmap_free = scatter_field_to_heightmap(field,
                                               shape,
                                               SCATTER_SHAPE_DISK,
                                               1.0f);
  EXPECT_GT(hmap_free(ix0, iy0), 0.05f);
  EXPECT_GT(hmap_free(ix1, iy1), 0.05f);
}

TEST(ScatterFieldTest, ToHeightmapRockMap)
{
  ScatterField field;
  // Two overlapping items at the same position to test collision clamping to
  // [0, 1]
  field.push_back(ScatterItem(0.5f, 0.5f, 0.0f, 0u, 0.2f));
  field.push_back(ScatterItem(0.5f, 0.5f, 0.0f, 0u, 0.2f));

  glm::ivec2 shape = {64, 64};
  Array      rock_map;
  Array      hmap = field.to_heightmap(shape,
                                  SCATTER_SHAPE_DISK,
                                  2.0f,
                                  std::nullopt,
                                  &rock_map);

  EXPECT_EQ(rock_map.shape.x, shape.x);
  EXPECT_EQ(rock_map.shape.y, shape.y);
  // Normalized rock map max amplitude must be <= 1.0f (clamped to [0, 1])
  EXPECT_NEAR(rock_map.max(), 1.0f, 0.02f);
  EXPECT_GE(rock_map.min(), 0.0f);
  EXPECT_LE(rock_map.max(), 1.0f);

  // Heightmap max should be ~ 2 * radius * ratio = 2 * 0.2 * 2.0 = 0.8
  EXPECT_NEAR(hmap.max(), 0.8f, 0.05f);

  // Test free function with rock_map
  Array rock_map_free;
  scatter_field_to_heightmap(field,
                             shape,
                             SCATTER_SHAPE_DISK,
                             2.0f,
                             std::nullopt,
                             &rock_map_free);
  EXPECT_NEAR(rock_map_free.max(), 1.0f, 0.02f);
  EXPECT_GE(rock_map_free.min(), 0.0f);
  EXPECT_LE(rock_map_free.max(), 1.0f);
}

TEST(ScatterFieldTest, ToImg8bit)
{
  ScatterField field;
  field.push_back(ScatterItem(0.5f, 0.5f, 0.0f, 0u, 0.2f));

  glm::ivec2           shape = {64, 64};
  std::vector<uint8_t> img = field.to_img_8bit(shape);

  EXPECT_EQ(img.size(), static_cast<size_t>(shape.x * shape.y * 3));

  // Background is not empty and has content
  bool has_non_bg = false;
  for (size_t i = 0; i < img.size(); i += 3)
  {
    // Default background is rgb(35, 35, 35)
    if (img[i] != 35 || img[i + 1] != 35 || img[i + 2] != 35)
    {
      has_non_bg = true;
      break;
    }
  }
  EXPECT_TRUE(has_non_bg);
}
