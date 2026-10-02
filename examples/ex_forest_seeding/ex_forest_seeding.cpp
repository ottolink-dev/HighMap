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

  hmap::Array density;
  hmap::Array exclusion;

  {
    auto gn = hmap::gradient_norm(z);

    auto cz = 1.f - hmap::threshold_smooth(z, 0.f, 0.9f);
    auto cg = 1.f - hmap::threshold_smooth(gn, 3.f / shape.x, 5.f / shape.x);
    auto ca = hmap::select_angle(z, 30.f, 90.f);

    auto cw = hmap::topographic_wetness_index(z);
    hmap::remap(cw);
    hmap::saturate_percentile(cw, 0.f, 0.95f);

    // combine linear density
    std::vector<const hmap::Array *> vec = {&cz, &cg, &ca, &cw};
    std::vector<float>               w = {1.f, 1.f, 0.8f, 0.5f};

    density = hmap::build_density_linear(vec, w);

    // exclusion map - exclude steep slopes and high mountain crests
    exclusion = hmap::threshold_smooth(gn, 4.f / shape.x, 6.f / shape.x);
  }

  // --- Multi-Species Forest Seeding
  // (Method A: k-means spatial partitioning)

  size_t species_count = 3;
  size_t tree_count = 5000;

  hmap::ForestSeedingOptions options;
  options.seed = static_cast<uint32_t>(seed);
  options.species_weights = {1.f, 1.f, 1.f};
  options.species_radii = {0.005f, 0.003f, 0.002f};
  options.exclusion_threshold = 0.5f;

  float cluster_randomness = 0.2f;

  // generate forest distribution using 2D inverse sampling and k-means
  // clustering
  hmap::Forest forest_kmeans = hmap::seed_forest_kmeans(species_count,
                                                        tree_count,
                                                        density,
                                                        exclusion,
                                                        cluster_randomness,
                                                        options);

  std::cout << "--- K-Means Forest ---\n"
            << forest_kmeans.to_string() << "\n\n";

  // sample z elevation from terrain
  forest_kmeans.set_elevation_from_terrain(z);

  // (Method B: parent-child cluster hierarchy)
  float  cluster_spread = 0.05f;
  size_t points_per_cluster = 32;

  hmap::Forest forest_clusters = hmap::seed_forest_clusters(species_count,
                                                            tree_count,
                                                            density,
                                                            exclusion,
                                                            cluster_spread,
                                                            points_per_cluster,
                                                            options);

  std::cout << "--- Clustered Forest ---\n"
            << forest_clusters.to_string() << "\n";

  // sample z elevation from terrain
  forest_clusters.set_elevation_from_terrain(z);

  // --- Export & Visualization

  z.dump("terrain.png");
  density.dump("density_linear.png");
  exclusion.dump("exclusion.png");

  forest_kmeans.to_png("forest_kmeans.png", shape, density);
  forest_clusters.to_png("forest_clusters.png", shape, density);

  hmap::export_banner_png("ex_forest_seeding.png",
                          {z, density, exclusion},
                          hmap::Cmap::JET);

  return 0;
}
