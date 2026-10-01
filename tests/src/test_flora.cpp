#include <cmath>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/flora.hpp"
#include "highmap/geometry/cloud.hpp"

#include <gtest/gtest.h>

using namespace hmap;

// --- Helper Functions

static constexpr float eps = 1e-5f;

static bool float_eq(float a, float b, float tol = eps)
{
  return std::abs(a - b) < tol;
}

// --- Tree Unit Tests

TEST(TreeTest, DefaultConstructor)
{
  Tree tree;
  EXPECT_TRUE(float_eq(tree.position.x, 0.f));
  EXPECT_TRUE(float_eq(tree.position.y, 0.f));
  EXPECT_TRUE(float_eq(tree.position.z, 0.f));
  EXPECT_EQ(tree.species_id, 0u);
  EXPECT_TRUE(float_eq(tree.radius, 1.0f));
}

TEST(TreeTest, ParameterizedConstructors)
{
  // 3D vector constructor
  Tree t1(glm::vec3(1.f, 2.f, 3.f), 2u, 4.5f);
  EXPECT_TRUE(float_eq(t1.position.x, 1.f));
  EXPECT_TRUE(float_eq(t1.position.y, 2.f));
  EXPECT_TRUE(float_eq(t1.position.z, 3.f));
  EXPECT_EQ(t1.species_id, 2u);
  EXPECT_TRUE(float_eq(t1.radius, 4.5f));

  // 2D vector constructor
  Tree t2(glm::vec2(5.f, 6.f), 1u, 2.0f);
  EXPECT_TRUE(float_eq(t2.position.x, 5.f));
  EXPECT_TRUE(float_eq(t2.position.y, 6.f));
  EXPECT_TRUE(float_eq(t2.position.z, 0.f));
  EXPECT_EQ(t2.species_id, 1u);

  // Component constructor
  Tree t3(7.f, 8.f, 9.f, 3u, 1.5f);
  EXPECT_TRUE(float_eq(t3.position.x, 7.f));
  EXPECT_TRUE(float_eq(t3.position.y, 8.f));
  EXPECT_TRUE(float_eq(t3.position.z, 9.f));
  EXPECT_EQ(t3.species_id, 3u);
}

TEST(TreeTest, ConversionsAndOperators)
{
  Tree t(1.f, 2.f, 3.f, 1u, 2.5f);

  Point p = t.to_point();
  EXPECT_TRUE(float_eq(p.x, 1.f));
  EXPECT_TRUE(float_eq(p.y, 2.f));
  EXPECT_TRUE(float_eq(p.v, 2.5f));

  glm::vec2 v2 = t.to_vec2();
  EXPECT_TRUE(float_eq(v2.x, 1.f));
  EXPECT_TRUE(float_eq(v2.y, 2.f));

  glm::vec3 v3 = t.to_vec3();
  EXPECT_TRUE(float_eq(v3.x, 1.f));
  EXPECT_TRUE(float_eq(v3.y, 2.f));
  EXPECT_TRUE(float_eq(v3.z, 3.f));

  Tree same(1.f, 2.f, 3.f, 1u, 2.5f);
  Tree diff(1.f, 2.f, 3.f, 2u, 2.5f);

  EXPECT_TRUE(t == same);
  EXPECT_FALSE(t == diff);
  EXPECT_TRUE(t != diff);
}

// --- Forest Unit Tests

TEST(ForestTest, DefaultConstructorEmpty)
{
  Forest forest;
  EXPECT_EQ(forest.size(), 0u);
  EXPECT_TRUE(forest.empty());
}

TEST(ForestTest, ContainerOperations)
{
  Forest forest;
  forest.push_back(Tree(1.f, 2.f, 0.f, 1u, 3.f));
  forest.emplace_back(4.f, 5.f, 0.f, 2u, 4.f);

  EXPECT_EQ(forest.size(), 2u);
  EXPECT_FALSE(forest.empty());

  EXPECT_TRUE(float_eq(forest[0].position.x, 1.f));
  EXPECT_TRUE(float_eq(forest.at(1).position.x, 4.f));
  EXPECT_TRUE(float_eq(forest.front().position.y, 2.f));
  EXPECT_TRUE(float_eq(forest.back().position.y, 5.f));

  EXPECT_THROW(forest.at(10), std::out_of_range);

  forest.clear();
  EXPECT_EQ(forest.size(), 0u);
  EXPECT_TRUE(forest.empty());
}

TEST(ForestTest, FromAndToCloud)
{
  Cloud cloud;
  cloud.push_back(Point(10.f, 20.f, 3.5f));
  cloud.push_back(Point(30.f, 40.f, 0.0f));

  Forest forest(cloud, 4u, 2.0f);
  ASSERT_EQ(forest.size(), 2u);

  EXPECT_TRUE(float_eq(forest[0].position.x, 10.f));
  EXPECT_TRUE(float_eq(forest[0].position.y, 20.f));
  EXPECT_TRUE(float_eq(forest[0].radius, 3.5f));
  EXPECT_EQ(forest[0].species_id, 4u);

  // second point had v = 0.0, so it gets default_radius 2.0
  EXPECT_TRUE(float_eq(forest[1].position.x, 30.f));
  EXPECT_TRUE(float_eq(forest[1].position.y, 40.f));
  EXPECT_TRUE(float_eq(forest[1].radius, 2.0f));
  EXPECT_EQ(forest[1].species_id, 4u);

  Cloud exported = forest.to_cloud();
  ASSERT_EQ(exported.size(), 2u);
  EXPECT_TRUE(float_eq(exported[0].x, 10.f));
  EXPECT_TRUE(float_eq(exported[0].y, 20.f));
  EXPECT_TRUE(float_eq(exported[0].v, 3.5f));
}

TEST(ForestTest, FilterBySpeciesAndSpeciesIds)
{
  Forest forest;
  forest.push_back(Tree(1.f, 1.f, 0.f, 1u, 2.f));
  forest.push_back(Tree(2.f, 2.f, 0.f, 2u, 3.f));
  forest.push_back(Tree(3.f, 3.f, 0.f, 1u, 2.f));
  forest.push_back(Tree(4.f, 4.f, 0.f, 3u, 1.f));

  auto species_ids = forest.get_species_ids();
  ASSERT_EQ(species_ids.size(), 3u);
  EXPECT_EQ(species_ids[0], 1u);
  EXPECT_EQ(species_ids[1], 2u);
  EXPECT_EQ(species_ids[2], 3u);

  Forest sp1 = forest.filter_by_species(1u);
  EXPECT_EQ(sp1.size(), 2u);
  EXPECT_EQ(sp1[0].species_id, 1u);
  EXPECT_EQ(sp1[1].species_id, 1u);

  Forest sp99 = forest.filter_by_species(99u);
  EXPECT_EQ(sp99.size(), 0u);
}

TEST(ForestTest, GetBbox)
{
  Forest    empty_forest;
  glm::vec4 bbox_empty = empty_forest.get_bbox();
  EXPECT_TRUE(float_eq(bbox_empty.x, 0.f));
  EXPECT_TRUE(float_eq(bbox_empty.y, 1.f));

  Forest forest;
  forest.push_back(Tree(10.f, 20.f, 0.f));
  forest.push_back(Tree(50.f, 5.f, 0.f));
  forest.push_back(Tree(30.f, 60.f, 0.f));

  glm::vec4 bbox = forest.get_bbox();
  EXPECT_TRUE(float_eq(bbox.x, 10.f));
  EXPECT_TRUE(float_eq(bbox.y, 50.f));
  EXPECT_TRUE(float_eq(bbox.z, 5.f));
  EXPECT_TRUE(float_eq(bbox.w, 60.f));
}

TEST(ForestTest, SetElevationFromTerrain)
{
  // 3x3 grid with linear slope along x
  Array elevation(glm::ivec2(3, 3), 0.f);
  for (int j = 0; j < 3; ++j)
  {
    for (int i = 0; i < 3; ++i)
    {
      elevation(i, j) = static_cast<float>(i) * 100.f; // 0, 100, 200
    }
  }

  glm::vec4 bbox = {0.f, 10.f, 0.f, 10.f};

  Forest forest;
  forest.push_back(Tree(0.f, 5.f, 0.f)); // at x = 0 -> elevation = 0
  forest.push_back(
      Tree(5.f, 5.f, 0.f)); // at x = 5 (midpoint) -> elevation = 100
  forest.push_back(Tree(10.f, 5.f, 0.f)); // at x = 10 -> elevation = 200

  forest.set_elevation_from_terrain(elevation, bbox);

  EXPECT_TRUE(float_eq(forest[0].position.z, 0.f));
  EXPECT_TRUE(float_eq(forest[1].position.z, 100.f));
  EXPECT_TRUE(float_eq(forest[2].position.z, 200.f));
}
