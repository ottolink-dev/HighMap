#include <iostream>

#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {1024, 1024};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 1;

  // --- Terrain Elevation

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::SIMPLEX2, shape, kw, seed);
  z = hmap::bulkify(z, hmap::PrimitiveType::PRIM_CUBIC_PULSE, 2.f);
  hmap::remap(z);

  // --- Environmental Suitability Criteria

  hmap::Array density = hmap::build_tree_density(z,
                                                 0.0f, // min_elev
                                                 0.9f, // max_elev
                                                 0.1f, // elev_transition_width
                                                 3.f / shape.x, // min_talus
                                                 5.f / shape.x, // max_talus
                                                 30.f,          // angle
                                                 90.f,          // angle_width
                                                 1.f,           // weight_elev
                                                 1.f,           // weight_talus
                                                 0.8f,          // weight_angle
                                                 0.5f);         // weight_twi

  // exclusion map - exclude steep slopes and high mountain crests
  hmap::Array gn = hmap::gradient_norm(z);
  hmap::Array exclusion = hmap::threshold_smooth(gn,
                                                 4.f / shape.x,
                                                 6.f / shape.x);

  // --- Multi-Species Forest Seeding
  // (Method A: k-means spatial partitioning)

  size_t species_count = 4;
  size_t tree_count = 10000;

  hmap::ForestSeedingOptions options;
  options.seed = static_cast<uint32_t>(seed);
  options.species = {
      hmap::Species(0, 0.003f, 1.0f), // id | radius | weight
      hmap::Species(1, 0.002f, 0.8f),
      hmap::Species(2, 0.001f, 0.5f),
      hmap::Species(3, 0.001f, 0.5f),
  };
  options.exclusion_threshold = 0.5f;

  float cluster_randomness = 0.f;

  // generate forest distribution using 2D inverse sampling and k-means
  // clustering
  hmap::Forest forest_kmeans = hmap::seed_forest_kmeans(species_count,
                                                        tree_count,
                                                        density,
                                                        exclusion,
                                                        cluster_randomness,
                                                        4,
                                                        options);

  std::cout << "--- K-Means Forest ---\n"
            << forest_kmeans.to_string() << "\n\n";

  // sample z elevation from terrain
  forest_kmeans.set_elevation_from_terrain(z);

  // (Method B: parent-child cluster hierarchy)
  float  cluster_spread = 0.05f;
  size_t points_per_cluster = 8;

  hmap::Forest forest_clusters = hmap::seed_forest_clusters(species_count,
                                                            tree_count,
                                                            density,
                                                            exclusion,
                                                            cluster_spread,
                                                            points_per_cluster,
                                                            options);

  // forest_clusters.densify();
  forest_clusters.reinforce_class_clusters();

  // sample z elevation from terrain
  forest_clusters.set_elevation_from_terrain(z);

  std::cout << "--- Clustered Forest ---\n"
            << forest_clusters.to_string() << "\n";

  // --- Grow forest
  // (Competition NN)
  hmap::Forest forest_grown_nn = hmap::grow_forest_competition_nn(
      forest_clusters,
      options.species,
      hmap::InteractionMatrix{},
      density,
      0.8f);

  std::cout << "--- Clustered Forest Grown NN ---\n"
            << forest_grown_nn.to_string() << "\n\n";

  // (Competition Voronoi)
  hmap::Forest forest_grown_voronoi = hmap::grow_forest_competition_voronoi(
      forest_clusters,
      options.species,
      hmap::InteractionMatrix{},
      density,
      0.8f);

  std::cout << "--- Clustered Forest Grown Voronoi ---\n"
            << forest_grown_voronoi.to_string() << "\n\n";

  // (Iterative Growth)
  hmap::Forest forest_grown_iterative = hmap::grow_forest_iterative(
      forest_clusters,
      options.species,
      15,
      0.1f,
      hmap::InteractionMatrix{},
      density,
      0.8f);

  std::cout << "--- Clustered Forest Grown Iterative ---\n"
            << forest_grown_iterative.to_string() << "\n\n";

  // (Strauss Soft-Core Growth)
  hmap::InteractionMatrix repulsion_distances =
      hmap::InteractionMatrix::from_species(options.species, 2.0f);
  hmap::InteractionMatrix repulsion_strengths =
      hmap::InteractionMatrix::uniform(species_count, 0.85f);

  hmap::Forest forest_thinned_soft_core = hmap::grow_forest_soft_core(
      forest_clusters,
      options.species,
      repulsion_distances,
      repulsion_strengths,
      0,
      density,
      0.8f);

  std::cout << "--- Clustered Forest Grown Soft Core ---\n"
            << forest_thinned_soft_core.to_string() << "\n";

  // --- Resulting densities

  auto density_soft_core = forest_thinned_soft_core.to_density_map(shape);

  // --- Export & Visualization

  z.dump("terrain.png");
  density.dump("density_linear.png");
  exclusion.dump("exclusion.png");

  density_soft_core.dump("density_soft_core.png");

  forest_kmeans.to_png("forest_kmeans.png", shape, density);
  forest_clusters.to_png("forest_clusters.png", shape, density);
  forest_grown_nn.to_png("forest_grown_nn.png", shape, z);
  forest_grown_voronoi.to_png("forest_grown_voronoi.png", shape, z);
  forest_grown_iterative.to_png("forest_grown_iterative.png", shape, z);
  forest_thinned_soft_core.to_png("forest_thinned_soft_core.png", shape, z);

  hmap::export_usd("scene_forest.usdc",
                   z,
                   forest_thinned_soft_core,
                   {},
                   {},
                   hmap::MeshType::TRI_OPTIMIZED,
                   0.25f);

  hmap::export_banner_png("ex_forest_seeding.png",
                          {z, density, exclusion},
                          hmap::Cmap::JET);

  return 0;
}
