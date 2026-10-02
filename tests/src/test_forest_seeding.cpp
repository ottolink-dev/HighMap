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
  options.species_radii = {0.5f, 1.0f, 1.5f};

  size_t species_count = 3;
  size_t tree_count = 100;

  Forest forest = seed_forest_kmeans(species_count,
                                     tree_count,
                                     density,
                                     exclusion,
                                     0.0f,
                                     options);
  EXPECT_EQ(forest.size(), tree_count);

  for (const auto &tree : forest)
  {
    EXPECT_GE(tree.position.x, 0.0f);
    EXPECT_LE(tree.position.x, 100.0f);
    EXPECT_GE(tree.position.y, 0.0f);
    EXPECT_LE(tree.position.y, 100.0f);
    EXPECT_LT(tree.species_id, species_count);
    EXPECT_FLOAT_EQ(tree.radius, options.species_radii[tree.species_id]);
  }

  auto species_ids = forest.get_species_ids();
  EXPECT_EQ(species_ids.size(), species_count);

  // test with cluster_randomness > 0
  Forest forest_rand = seed_forest_kmeans(species_count,
                                          tree_count,
                                          density,
                                          exclusion,
                                          0.5f,
                                          options);
  EXPECT_EQ(forest_rand.size(), tree_count);
  EXPECT_EQ(forest_rand.get_species_ids().size(), species_count);

  Forest f0 = seed_forest_kmeans(3, 500, density, exclusion, 0.0f, options);
  Forest f1 = seed_forest_kmeans(3, 500, density, exclusion, 0.5f, options);
  Forest f2 = seed_forest_kmeans(3, 500, density, exclusion, 1.0f, options);

  // test with empty species_radii
  ForestSeedingOptions options_empty_radii = options;
  options_empty_radii.species_radii = {};
  Forest forest_empty_r = seed_forest_kmeans(species_count,
                                             tree_count,
                                             density,
                                             exclusion,
                                             0.2f,
                                             options_empty_radii);
  EXPECT_EQ(forest_empty_r.size(), tree_count);
  for (const auto &tree : forest_empty_r)
  {
    EXPECT_FALSE(std::isnan(tree.position.x));
    EXPECT_FALSE(std::isnan(tree.position.y));
    EXPECT_FALSE(std::isnan(tree.position.z));
    EXPECT_FALSE(std::isnan(tree.radius));
    EXPECT_FLOAT_EQ(tree.radius, 1e-3f);
  }
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
  Forest forest_km = seed_forest_kmeans(2,
                                        50,
                                        density,
                                        exclusion,
                                        0.0f,
                                        options);
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
  options.species_weights = {9.0f, 1.0f};

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
