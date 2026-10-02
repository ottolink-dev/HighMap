#include <cmath>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/flora.hpp"

#include <gtest/gtest.h>

using namespace hmap;

// --- Helper Functions

static constexpr float eps = 1e-4f;

static bool float_eq(float a, float b, float tol = eps)
{
  return std::abs(a - b) < tol;
}

// --- InteractionMatrix Tests

TEST(InteractionMatrixTest, BasicOperations)
{
  InteractionMatrix mat(3, 1.5f);
  EXPECT_EQ(mat.size, 3u);
  EXPECT_TRUE(float_eq(mat.get(0, 0), 1.5f));
  EXPECT_TRUE(float_eq(mat.get(1, 2), 1.5f));

  mat.set(0, 1, 3.0f);
  EXPECT_TRUE(float_eq(mat.get(0, 1), 3.0f));
  EXPECT_TRUE(float_eq(mat.get(1, 0), 1.5f));

  mat.set_symmetric(1, 2, 4.2f);
  EXPECT_TRUE(float_eq(mat.get(1, 2), 4.2f));
  EXPECT_TRUE(float_eq(mat.get(2, 1), 4.2f));

  mat.fill(0.5f);
  for (size_t i = 0; i < 3; ++i)
    for (size_t j = 0; j < 3; ++j)
      EXPECT_TRUE(float_eq(mat.get(i, j), 0.5f));
}

TEST(InteractionMatrixTest, FactoryMethods)
{
  std::vector<float> diag = {1.0f, 2.0f, 3.0f};
  InteractionMatrix  mat = InteractionMatrix::diagonal(diag, 0.1f);

  EXPECT_EQ(mat.size, 3u);
  EXPECT_TRUE(float_eq(mat.get(0, 0), 1.0f));
  EXPECT_TRUE(float_eq(mat.get(1, 1), 2.0f));
  EXPECT_TRUE(float_eq(mat.get(2, 2), 3.0f));
  EXPECT_TRUE(float_eq(mat.get(0, 1), 0.1f));
  EXPECT_TRUE(float_eq(mat.get(2, 0), 0.1f));

  InteractionMatrix uni = InteractionMatrix::uniform(4, 2.5f);
  EXPECT_EQ(uni.size, 4u);
  EXPECT_TRUE(float_eq(uni.get(3, 3), 2.5f));

  // random matrix with variation around 1.0
  float             offset = 0.2f;
  InteractionMatrix rnd = InteractionMatrix::random(3, 42, offset, true);
  EXPECT_EQ(rnd.size, 3u);
  for (size_t i = 0; i < 3; ++i)
  {
    for (size_t j = 0; j < 3; ++j)
    {
      float val = rnd.get(i, j);
      EXPECT_GE(val, 1.0f - offset);
      EXPECT_LE(val, 1.0f + offset);
      EXPECT_TRUE(float_eq(rnd.get(i, j), rnd.get(j, i))); // symmetric
    }
  }

  // from_radii
  std::vector<float> radii = {1.0f, 2.0f, 0.5f};
  InteractionMatrix  from_r = InteractionMatrix::from_radii(radii, 2.0f);
  EXPECT_EQ(from_r.size, 3u);
  EXPECT_TRUE(float_eq(from_r.get(0, 0), 4.0f)); // 2.0 * (1.0 + 1.0)
  EXPECT_TRUE(float_eq(from_r.get(0, 1), 6.0f)); // 2.0 * (1.0 + 2.0)
  EXPECT_TRUE(float_eq(from_r.get(1, 2), 5.0f)); // 2.0 * (2.0 + 0.5)

  // from_species
  std::vector<Species> spec_list = {Species(0, 1.0f), Species(1, 2.0f)};
  InteractionMatrix from_sp = InteractionMatrix::from_species(spec_list, 1.0f);
  EXPECT_EQ(from_sp.size, 2u);
  EXPECT_TRUE(float_eq(from_sp.get(0, 0), 2.0f));
  EXPECT_TRUE(float_eq(from_sp.get(0, 1), 3.0f));
}

// --- Species Class Tests

TEST(SpeciesTest, DefaultAndCustomConstructors)
{
  Species s0;
  EXPECT_EQ(s0.id, 0u);
  EXPECT_TRUE(float_eq(s0.radius, HMAP_DEFAULT_TREE_RADIUS));
  EXPECT_TRUE(float_eq(s0.weight, 1.0f));
  EXPECT_TRUE(float_eq(s0.radius_min, 0.5f * HMAP_DEFAULT_TREE_RADIUS));
  EXPECT_TRUE(float_eq(s0.radius_max, 1.5f * HMAP_DEFAULT_TREE_RADIUS));
  EXPECT_TRUE(float_eq(s0.competition_factor, 0.4f));

  Species s1(2u, 0.01f, 2.5f, 0.005f, 0.02f, 0.35f, "Oak");
  EXPECT_EQ(s1.id, 2u);
  EXPECT_EQ(s1.name, "Oak");
  EXPECT_TRUE(float_eq(s1.radius, 0.01f));
  EXPECT_TRUE(float_eq(s1.weight, 2.5f));
  EXPECT_TRUE(float_eq(s1.radius_min, 0.005f));
  EXPECT_TRUE(float_eq(s1.radius_max, 0.02f));
  EXPECT_TRUE(float_eq(s1.competition_factor, 0.35f));

  // auto default radius_min / radius_max when negative
  Species s2(1u, 0.02f, 1.0f, -1.0f, -1.0f);
  EXPECT_TRUE(float_eq(s2.radius_min, 0.01f));
  EXPECT_TRUE(float_eq(s2.radius_max, 0.03f));
}

// --- Forest Seeding Tests

TEST(ForestSeedingTest, EmptyAndInvalidInputs)
{
  Array density(glm::ivec2(16, 16), 1.0f);
  Array exclusion(glm::ivec2(16, 16), 0.0f);
  Array empty_array;

  // 0 species
  Forest f1 = seed_forest_clusters(0, 50, density, exclusion);
  EXPECT_EQ(f1.size(), 0u);

  // 0 trees
  Forest f2 = seed_forest_clusters(2, 0, density, exclusion);
  EXPECT_EQ(f2.size(), 0u);

  // empty density
  Forest f3 = seed_forest_clusters(2, 50, empty_array, exclusion);
  EXPECT_EQ(f3.size(), 0u);

  // kmeans seeding 0 species
  Forest f4 = seed_forest_kmeans(0, 50, density, exclusion);
  EXPECT_EQ(f4.size(), 0u);

  // kmeans seeding 0 trees
  Forest f5 = seed_forest_kmeans(2, 0, density, exclusion);
  EXPECT_EQ(f5.size(), 0u);

  // kmeans seeding empty density
  Forest f6 = seed_forest_kmeans(2, 50, empty_array, exclusion);
  EXPECT_EQ(f6.size(), 0u);

  // shape mismatch between density and exclusion
  Array  mismatched_exclusion(glm::ivec2(8, 8), 0.0f);
  Forest f7 = seed_forest_clusters(2, 50, density, mismatched_exclusion);
  EXPECT_EQ(f7.size(), 0u);

  Forest f8 = seed_forest_kmeans(2, 50, density, mismatched_exclusion);
  EXPECT_EQ(f8.size(), 0u);
}

TEST(ForestSeedingTest, BasicClusterSeeding)
{
  Array density(glm::ivec2(32, 32), 1.0f);
  Array exclusion(glm::ivec2(32, 32), 0.0f);

  ForestSeedingOptions options;
  options.seed = 42;
  options.bbox = {0.0f, 100.0f, 0.0f, 100.0f};

  size_t species_count = 3;
  size_t tree_count = 60;

  Forest forest = seed_forest_clusters(species_count,
                                       tree_count,
                                       density,
                                       exclusion,
                                       0.05f,
                                       16,
                                       options);
  EXPECT_GT(forest.size(), 0u);
  EXPECT_LE(forest.size(), tree_count);

  for (const auto &tree : forest)
  {
    EXPECT_GE(tree.position.x, 0.0f);
    EXPECT_LE(tree.position.x, 100.0f);
    EXPECT_GE(tree.position.y, 0.0f);
    EXPECT_LE(tree.position.y, 100.0f);
    EXPECT_LT(tree.species_id, species_count);
  }

  auto species_ids = forest.get_species_ids();
  EXPECT_GT(species_ids.size(), 0u);
}

TEST(ForestSeedingTest, BasicKMeansSeeding)
{
  Array density(glm::ivec2(32, 32), 1.0f);
  Array exclusion(glm::ivec2(32, 32), 0.0f);

  ForestSeedingOptions options;
  options.seed = 42;
  options.bbox = {0.0f, 100.0f, 0.0f, 100.0f};
  options.species = {Species(0, 0.5f), Species(1, 1.0f), Species(2, 1.5f)};

  size_t species_count = 3;
  size_t tree_count = 100;

  Forest forest = seed_forest_kmeans(species_count,
                                     tree_count,
                                     density,
                                     exclusion,
                                     0.0f,
                                     4,
                                     options);
  EXPECT_EQ(forest.size(), tree_count);

  for (const auto &tree : forest)
  {
    EXPECT_GE(tree.position.x, 0.0f);
    EXPECT_LE(tree.position.x, 100.0f);
    EXPECT_GE(tree.position.y, 0.0f);
    EXPECT_LE(tree.position.y, 100.0f);
    EXPECT_LT(tree.species_id, species_count);
    EXPECT_FLOAT_EQ(tree.radius, options.species[tree.species_id].radius);
  }

  auto species_ids = forest.get_species_ids();
  EXPECT_EQ(species_ids.size(), species_count);

  // test with cluster_randomness > 0
  Forest forest_rand = seed_forest_kmeans(species_count,
                                          tree_count,
                                          density,
                                          exclusion,
                                          0.5f,
                                          4,
                                          options);
  EXPECT_EQ(forest_rand.size(), tree_count);
  EXPECT_EQ(forest_rand.get_species_ids().size(), species_count);

  Forest f0 = seed_forest_kmeans(3, 500, density, exclusion, 0.0f, 4, options);
  Forest f1 = seed_forest_kmeans(3, 500, density, exclusion, 0.5f, 4, options);
  Forest f2 = seed_forest_kmeans(3, 500, density, exclusion, 1.0f, 4, options);

  // test with empty species definitions
  ForestSeedingOptions options_empty_species = options;
  options_empty_species.species = {};
  Forest forest_empty_r = seed_forest_kmeans(species_count,
                                             tree_count,
                                             density,
                                             exclusion,
                                             0.2f,
                                             4,
                                             options_empty_species);
  EXPECT_EQ(forest_empty_r.size(), tree_count);
  for (const auto &tree : forest_empty_r)
  {
    EXPECT_FALSE(std::isnan(tree.position.x));
    EXPECT_FALSE(std::isnan(tree.position.y));
    EXPECT_FALSE(std::isnan(tree.position.z));
    EXPECT_FALSE(std::isnan(tree.radius));
    EXPECT_FLOAT_EQ(tree.radius, 1e-3f);
  }

  // test with custom k_neighbors
  Forest forest_k8 = seed_forest_kmeans(species_count,
                                        tree_count,
                                        density,
                                        exclusion,
                                        0.0f,
                                        8,
                                        options);
  EXPECT_EQ(forest_k8.size(), tree_count);
  EXPECT_EQ(forest_k8.get_species_ids().size(), species_count);
}

TEST(ForestSeedingTest, DensityAdherence)
{
  // Density only on right half (x in [0.5, 1.0])
  Array density(glm::ivec2(32, 32), 0.0f);
  for (int j = 0; j < 32; ++j)
    for (int i = 16; i < 32; ++i)
      density(i, j) = 1.0f;

  Array exclusion(glm::ivec2(32, 32), 0.0f);

  ForestSeedingOptions options;
  options.seed = 123;
  options.bbox = {0.0f, 1.0f, 0.0f, 1.0f};

  Forest forest =
      seed_forest_clusters(2, 40, density, exclusion, 0.03f, 8, options);
  EXPECT_GT(forest.size(), 0u);

  // Most cluster centers and offspring should be on the right half (allowing
  // small cluster spread around 0.5)
  size_t right_count = 0;
  for (const auto &tree : forest)
  {
    if (tree.position.x >= 0.4f) right_count++;
  }

  EXPECT_GT(right_count, forest.size() * 8 / 10);
}

TEST(ForestSeedingTest, ExclusionMapMasking)
{
  Array density(glm::ivec2(32, 32), 1.0f);
  Array exclusion(glm::ivec2(32, 32), 0.0f);

  // Exclude entire top half (y in [0.5, 1.0])
  for (int j = 16; j < 32; ++j)
    for (int i = 0; i < 32; ++i)
      exclusion(i, j) = 1.0f;

  ForestSeedingOptions options;
  options.seed = 999;
  options.bbox = {0.0f, 1.0f, 0.0f, 1.0f};
  options.exclusion_threshold = 0.5f;

  Forest forest =
      seed_forest_clusters(2, 50, density, exclusion, 0.05f, 16, options);
  EXPECT_GT(forest.size(), 0u);

  for (const auto &tree : forest)
  {
    EXPECT_LT(tree.position.y, 0.55f);
  }

  // Also test KMeans seeding with exclusion
  Forest forest_km =
      seed_forest_kmeans(2, 50, density, exclusion, 0.0f, 4, options);
  EXPECT_GT(forest_km.size(), 0u);

  for (const auto &tree : forest_km)
  {
    EXPECT_LT(tree.position.y, 0.55f);
  }
}

TEST(ForestSeedingTest, SpeciesWeights)
{
  Array density(glm::ivec2(32, 32), 1.0f);
  Array exclusion(glm::ivec2(32, 32), 0.0f);

  ForestSeedingOptions options;
  options.seed = 777;
  options.bbox = {0.0f, 10.0f, 0.0f, 10.0f};
  // 90% species 0, 10% species 1
  options.species = {Species(0, 0.001f, 9.0f), Species(1, 0.001f, 1.0f)};

  Forest forest =
      seed_forest_clusters(2, 50, density, exclusion, 0.05f, 16, options);
  EXPECT_GT(forest.size(), 0u);

  Forest sp0 = forest.filter_by_species(0);
  Forest sp1 = forest.filter_by_species(1);

  EXPECT_GT(sp0.size(), sp1.size());
}

TEST(ForestSeedingTest, SoftCoreThinning)
{
  // Create a tight cluster of overlapping trees
  Forest dense_forest;
  for (int i = 0; i < 20; ++i)
  {
    float offset = static_cast<float>(i) * 0.01f;
    dense_forest.push_back(Tree(5.0f + offset, 5.0f + offset, 0.0f, 0u, 1.0f));
  }

  InteractionMatrix dist = InteractionMatrix::uniform(1, 0.5f);
  InteractionMatrix strength = InteractionMatrix::uniform(1, 1.0f);

  Forest thinned = thin_forest_soft_core(dense_forest,
                                         dist,
                                         strength,
                                         20,
                                         42,
                                         {0.f, 10.f, 0.f, 10.f});

  // Since all 20 trees are within 0.2 units of each other and repulsion
  // distance is 0.5 with strength 1.0, thinning should significantly reduce the
  // tree count
  EXPECT_LT(thinned.size(), dense_forest.size());
  EXPECT_GT(thinned.size(), 0u);
}

// --- Nearest-Neighbor Competition Growth Tests

TEST(ForestGrowthTest, CompetitionNN)
{
  // Create two trees separated by distance = 2.0
  Forest forest;
  forest.push_back(Tree(0.0f, 0.0f, 0.0f, 0u, 0.01f));
  forest.push_back(Tree(2.0f, 0.0f, 0.0f, 1u, 0.01f));

  // Species 0: alpha = 0.4, r_min = 0.1, r_max = 2.0 -> r_est = 0.4 * 2.0 = 0.8
  // Species 1: alpha = 0.3, r_min = 0.1, r_max = 0.5 -> r_est = 0.3 * 2.0 = 0.6
  // -> clamped to 0.5
  std::vector<Species> species = {
      Species(0u, 0.8f, 1.0f, 0.1f, 2.0f, 0.4f),
      Species(1u, 0.5f, 1.0f, 0.1f, 0.5f, 0.3f),
  };

  Forest grown = grow_forest_competition_nn(forest, species);
  ASSERT_EQ(grown.size(), 2u);
  EXPECT_TRUE(float_eq(grown[0].radius, 0.8f));
  EXPECT_TRUE(float_eq(grown[1].radius, 0.5f));
}

TEST(ForestGrowthTest, CompetitionNNPruneUnviable)
{
  // 3 trees: tree 0 and tree 1 are very close (dist = 0.02), tree 2 is far
  // (dist = 5.0)
  Forest forest;
  forest.push_back(Tree(0.0f, 0.0f, 0.0f, 0u, 0.01f));
  forest.push_back(Tree(0.02f, 0.0f, 0.0f, 0u, 0.01f));
  forest.push_back(Tree(5.0f, 0.0f, 0.0f, 0u, 0.01f));

  // Species 0: alpha = 0.5, r_min = 0.1, r_max = 1.0
  // For trees 0 & 1: r_est = 0.5 * 0.02 = 0.01 < r_min (0.1) -> should be
  // culled if prune_unviable = true For tree 2: nearest neighbor is tree 1 at
  // dist ~4.98 -> r_est = 0.5 * 4.98 = 2.49 -> clamped to 1.0
  std::vector<Species> species = {
      Species(0u, 0.5f, 1.0f, 0.1f, 1.0f, 0.5f),
  };

  Forest grown = grow_forest_competition_nn(forest,
                                            species,
                                            {},
                                            {},
                                            1.0f,
                                            true);
  ASSERT_EQ(grown.size(), 1u);
  EXPECT_TRUE(float_eq(grown[0].position.x, 5.0f));
  EXPECT_TRUE(float_eq(grown[0].radius, 1.0f));
}

TEST(ForestGrowthTest, CompetitionNNMaxRadiusScale)
{
  // 2 widely separated trees (dist = 10.0)
  Forest forest;
  forest.push_back(Tree(0.0f, 0.5f, 0.0f, 0u, 0.01f)); // at x = 0.0
  forest.push_back(Tree(1.0f, 0.5f, 0.0f, 0u, 0.01f)); // at x = 1.0

  // Species 0: r_min = 0.1, r_max = 1.0, alpha = 0.5 -> r_est = 0.5 * 1.0 = 0.5
  std::vector<Species> species = {
      Species(0u, 0.5f, 1.0f, 0.1f, 1.0f, 0.5f),
  };

  // Scale map: left (x=0) has value 0.0 => r_max_eff = r_min = 0.1
  //            right (x=1) has value 1.0 => r_max_eff = r_max = 1.0
  Array scale_map({2, 1}, 0.0f);
  scale_map(0, 0) = 0.0f;
  scale_map(1, 0) = 1.0f;

  glm::vec4 bbox = {0.0f, 1.0f, 0.0f, 1.0f};

  // Full strength (1.0)
  Forest grown_full = grow_forest_competition_nn(forest,
                                                 species,
                                                 {},
                                                 scale_map,
                                                 1.0f,
                                                 false,
                                                 bbox);

  ASSERT_EQ(grown_full.size(), 2u);
  // Tree at x = 0.0: r_max_eff = 0.1 + 0.0 * (1.0 - 0.1) = 0.1 -> clamped to
  // 0.1
  EXPECT_TRUE(float_eq(grown_full[0].radius, 0.1f));
  // Tree at x = 1.0: r_max_eff = 0.1 + 1.0 * (1.0 - 0.1) = 1.0 -> r_est = 0.5
  // fits
  EXPECT_TRUE(float_eq(grown_full[1].radius, 0.5f));

  // Zero strength (0.0): scale array ignored -> r_max_eff = r_max = 1.0
  // everywhere
  Forest grown_zero = grow_forest_competition_nn(forest,
                                                 species,
                                                 {},
                                                 scale_map,
                                                 0.0f,
                                                 false,
                                                 bbox);

  ASSERT_EQ(grown_zero.size(), 2u);
  EXPECT_TRUE(float_eq(grown_zero[0].radius, 0.5f));
  EXPECT_TRUE(float_eq(grown_zero[1].radius, 0.5f));

  // Half strength (0.5): at x = 0.0, s_raw = 0.0 -> s_eff = 1.0 + 0.5 * (0 - 1)
  // = 0.5 r_max_eff = 0.1 + 0.5 * (1.0 - 0.1) = 0.55 -> r_est = 0.5 fits
  Forest grown_half = grow_forest_competition_nn(forest,
                                                 species,
                                                 {},
                                                 scale_map,
                                                 0.5f,
                                                 false,
                                                 bbox);
  ASSERT_EQ(grown_half.size(), 2u);
  EXPECT_TRUE(float_eq(grown_half[0].radius, 0.5f));
}

TEST(ForestGrowthTest, CompetitionVoronoiBasic)
{
  // 4 points forming a square [0, 10] x [0, 10] -> total area = 100
  Forest forest;
  forest.push_back(Tree(0.0f, 0.0f, 0.0f, 0u, 0.01f));
  forest.push_back(Tree(10.0f, 0.0f, 0.0f, 0u, 0.01f));
  forest.push_back(Tree(10.0f, 10.0f, 0.0f, 0u, 0.01f));
  forest.push_back(Tree(0.0f, 10.0f, 0.0f, 0u, 0.01f));

  // Species 0: alpha = 0.5, r_min = 0.1, r_max = 5.0
  std::vector<Species> species = {
      Species(0u, 1.0f, 1.0f, 0.1f, 5.0f, 0.5f),
  };

  Forest grown = grow_forest_competition_voronoi(forest, species);
  ASSERT_EQ(grown.size(), 4u);

  // Each corner vertex receives non-zero vertex area and positive radius within
  // [r_min, r_max]
  for (const auto &tree : grown)
  {
    EXPECT_GE(tree.radius, 0.1f);
    EXPECT_LE(tree.radius, 5.0f);
  }
}

TEST(ForestGrowthTest, CompetitionVoronoiPrune)
{
  // Dense cluster of 3 points with very small triangle area, plus 1 distant
  // point
  Forest forest;
  forest.push_back(Tree(0.0f, 0.0f, 0.0f, 0u, 0.01f));
  forest.push_back(Tree(0.01f, 0.0f, 0.0f, 0u, 0.01f));
  forest.push_back(Tree(0.0f, 0.01f, 0.0f, 0u, 0.01f));
  forest.push_back(Tree(10.0f, 10.0f, 0.0f, 0u, 0.01f));

  // Species 0: alpha = 1.0, r_min = 0.5, r_max = 5.0
  // Dense cluster points have area ~ 1e-4 / 3, r_raw ~ sqrt(1e-4 / (3*pi)) ~
  // 0.003 < r_min (0.5)
  std::vector<Species> species = {
      Species(0u, 1.0f, 1.0f, 0.5f, 5.0f, 1.0f),
  };

  Forest grown = grow_forest_competition_voronoi(forest,
                                                 species,
                                                 {},
                                                 {},
                                                 1.0f,
                                                 true);
  // Dense points should be pruned
  EXPECT_LT(grown.size(), forest.size());
}

TEST(ForestGrowthTest, CompetitionVoronoiMaxRadiusScale)
{
  // 4 corners of a large square
  Forest forest;
  forest.push_back(Tree(0.0f, 0.0f, 0.0f, 0u, 0.01f));
  forest.push_back(Tree(1.0f, 0.0f, 0.0f, 0u, 0.01f));
  forest.push_back(Tree(1.0f, 1.0f, 0.0f, 0u, 0.01f));
  forest.push_back(Tree(0.0f, 1.0f, 0.0f, 0u, 0.01f));

  std::vector<Species> species = {
      Species(0u, 0.5f, 1.0f, 0.1f, 2.0f, 1.0f),
  };

  Array scale_map({2, 1}, 0.0f);
  scale_map(0, 0) = 0.0f; // left side scale = 0 -> r_max_eff = r_min = 0.1
  scale_map(1, 0) = 1.0f; // right side scale = 1 -> r_max_eff = r_max = 2.0

  glm::vec4 bbox = {0.0f, 1.0f, 0.0f, 1.0f};

  Forest grown_full = grow_forest_competition_voronoi(forest,
                                                      species,
                                                      {},
                                                      scale_map,
                                                      1.0f,
                                                      false,
                                                      bbox);
  ASSERT_EQ(grown_full.size(), 4u);

  // Left trees (x=0) must be clamped to r_min = 0.1
  for (const auto &t : grown_full)
  {
    if (t.position.x == 0.0f)
    {
      EXPECT_TRUE(float_eq(t.radius, 0.1f));
    }
  }
}

TEST(ForestGrowthTest, CompetitionVoronoiCollisionPruning)
{
  // 4 trees in a rectangle where two trees are close and their grown radii will
  // overlap
  Forest forest;
  forest.push_back(Tree(0.0f, 0.0f, 0.0f, 0u, 0.01f));
  forest.push_back(
      Tree(1.0f, 0.0f, 0.0f, 1u, 0.01f)); // species 1 will have smaller radius
  forest.push_back(Tree(10.0f, 10.0f, 0.0f, 0u, 0.01f));
  forest.push_back(Tree(0.0f, 10.0f, 0.0f, 0u, 0.01f));

  // Species 0: alpha = 1.0, r_min = 0.5, r_max = 5.0
  // Species 1: alpha = 0.1, r_min = 0.1, r_max = 0.2
  // Distance between (0,0) and (1,0) is 1.0. Tree 0 will grow larger, e.g.
  // radius > 0.9, causing r0 + r1 > 1.0 -> tree 1 (smaller) should be removed
  // due to collision.
  std::vector<Species> species = {
      Species(0u, 1.0f, 1.0f, 0.5f, 5.0f, 1.0f),
      Species(1u, 0.2f, 1.0f, 0.1f, 0.2f, 0.1f),
  };

  Forest grown = grow_forest_competition_voronoi(forest,
                                                 species,
                                                 {},
                                                 {},
                                                 1.0f,
                                                 true);
  // Ensure no colliding pairs remain in grown forest
  for (size_t i = 0; i < grown.size(); ++i)
  {
    for (size_t j = i + 1; j < grown.size(); ++j)
    {
      float dx = grown[i].position.x - grown[j].position.x;
      float dy = grown[i].position.y - grown[j].position.y;
      float dist = std::sqrt(dx * dx + dy * dy);
      EXPECT_GE(dist, grown[i].radius + grown[j].radius - 1e-5f);
    }
  }
}

// --- Forest::densify Tests

TEST(ForestTest, DensifyBasic)
{
  // A single right triangle with vertices at (0, 0), (2, 0), (0, 2)
  // Circumcenter is at (1, 1)
  Forest forest;
  forest.push_back(Tree(0.0f, 0.0f, 0.0f, 1u, 0.5f));
  forest.push_back(Tree(2.0f, 0.0f, 0.0f, 1u, 0.5f));
  forest.push_back(Tree(0.0f, 2.0f, 0.0f, 2u, 0.5f));

  forest.densify(0.05f);

  // Expect 3 original trees + 1 new tree at (1, 1)
  ASSERT_EQ(forest.size(), 4u);

  const Tree &new_tree = forest.back();
  EXPECT_TRUE(float_eq(new_tree.position.x, 1.0f));
  EXPECT_TRUE(float_eq(new_tree.position.y, 1.0f));
  EXPECT_EQ(new_tree.species_id,
            1u); // species 1 has count 2 vs species 2 count 1
  EXPECT_TRUE(float_eq(new_tree.radius, 0.05f));
}

TEST(ForestTest, DensifySmallForest)
{
  Forest forest;
  forest.push_back(Tree(0.0f, 0.0f, 0.0f, 0u, 0.1f));
  forest.push_back(Tree(1.0f, 1.0f, 0.0f, 0u, 0.1f));

  // Fewer than 3 trees: densify should be a no-op
  forest.densify();
  EXPECT_EQ(forest.size(), 2u);
}

TEST(ForestTest, DensifyGridSpeciesMajority)
{
  // 4 trees forming a square with species 0, 0, 1, 2
  Forest forest;
  forest.push_back(Tree(0.0f, 0.0f, 0.0f, 0u, 0.1f));
  forest.push_back(Tree(1.0f, 0.0f, 0.0f, 0u, 0.1f));
  forest.push_back(Tree(1.0f, 1.0f, 0.0f, 1u, 0.1f));
  forest.push_back(Tree(0.0f, 1.0f, 0.0f, 2u, 0.1f));

  size_t initial_count = forest.size();
  forest.densify(0.02f);

  // Delaunay of 4 points forms 2 triangles -> 2 new trees added
  EXPECT_EQ(forest.size(), initial_count + 2);

  for (size_t i = initial_count; i < forest.size(); ++i)
  {
    EXPECT_TRUE(float_eq(forest[i].radius, 0.02f));
  }
}

// --- Forest::reinforce_species_clusters Tests

TEST(ForestTest, ReinforceSpeciesClustersBasic)
{
  // A center tree surrounded by 4 neighbor trees with species 1
  // Center tree with species 2 should be converted to species 1
  Forest forest;
  forest.push_back(Tree(0.0f, 0.0f, 0.0f, 2u, 0.1f)); // center
  forest.push_back(Tree(1.0f, 0.0f, 0.0f, 1u, 0.1f));
  forest.push_back(Tree(-1.0f, 0.0f, 0.0f, 1u, 0.1f));
  forest.push_back(Tree(0.0f, 1.0f, 0.0f, 1u, 0.1f));
  forest.push_back(Tree(0.0f, -1.0f, 0.0f, 1u, 0.1f));

  forest.reinforce_species_clusters(1, 4, true);

  // Center tree should now have species 1 (4 votes vs 1 vote)
  EXPECT_EQ(forest[0].species_id, 1u);
}

TEST(ForestTest, ReinforceSpeciesClustersMultiIteration)
{
  // Line of trees: [1, 1, 2, 2, 2]
  Forest forest;
  forest.push_back(Tree(0.0f, 0.0f, 0.0f, 1u, 0.1f));
  forest.push_back(Tree(1.0f, 0.0f, 0.0f, 1u, 0.1f));
  forest.push_back(Tree(2.0f, 0.0f, 0.0f, 2u, 0.1f));
  forest.push_back(Tree(3.0f, 0.0f, 0.0f, 2u, 0.1f));
  forest.push_back(Tree(4.0f, 0.0f, 0.0f, 2u, 0.1f));

  // Trivial edge cases
  forest.reinforce_species_clusters(0, 2);
  EXPECT_EQ(forest[0].species_id, 1u);
  EXPECT_EQ(forest[2].species_id, 2u);

  // With k=2 and multiple iterations
  forest.reinforce_species_clusters(3, 2, true);
  // All should have valid species
  for (const auto &t : forest)
  {
    EXPECT_TRUE(t.species_id == 1u || t.species_id == 2u);
  }
}
