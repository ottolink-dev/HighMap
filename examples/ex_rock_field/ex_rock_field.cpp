#include <iostream>

#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {512, 512};
  glm::vec2  kw = {4.f, 4.f};
  int        seed = 42;

  // --- 1. Generate Terrain Elevation and Slope

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::SIMPLEX2, shape, kw, seed);
  z = hmap::bulkify(z, hmap::PrimitiveType::PRIM_CUBIC_PULSE, 2.f);
  hmap::remap(z);

  // --- 2. Build Scree / Talus Rock Density

  // compute slope / gradient magnitude for geomorphological sorting
  hmap::Array slope = hmap::gradient_norm(z);

  // scree fields naturally accumulate in steep-to-moderate chute
  // regions
  hmap::Array density = hmap::select_range(slope,
                                           1.f / shape.x,
                                           2.f / shape.x,
                                           0.5f / shape.x);
  density.infos();

  // exclusion map (e.g., lakes or very flat valleys)
  hmap::Array exclusion(
      shape); //  1.0f - hmap::threshold_smooth(z, 0.0f, 0.2f);

  // --- 3. Seed Rock Field with Pareto Power-Law Distribution

  hmap::RockSeedingOptions options;
  options.seed = static_cast<uint32_t>(seed);
  options.bbox = {0.f, 1.f, 0.f, 1.f};
  options.distribution = hmap::RockDistribution(0u,
                                                0.001f,
                                                0.005f,
                                                2.f); // Pareto alpha = 2.0

  size_t          target_count = 1500;
  hmap::RockField rocks = hmap::seed_rock_field(target_count,
                                                density,
                                                exclusion,
                                                options);

  std::cout << "=== Initial Rock Field ===\n";
  std::cout << rocks.to_string() << "\n\n";

  // --- 4. Apply Gravitational Slope Sorting

  // larger boulders roll down steeper slopes and accumulate at
  // gentler bases
  rocks.apply_slope_sorting(slope, 0.8f, options.bbox);

  // --- 5. Interstitial Pebble Packing

  // pack fine gravel / interstitial stones (new class_id = 1) in the voids
  rocks.pack_interstitial_rocks(500,
                                density,
                                0.0005f,
                                0.0015f,
                                1u,
                                1337,
                                options.bbox);

  // --- 6. Prune Collisions and Project Elevation

  // resolve physical overlaps between boulders
  rocks.prune_collisions();

  // project rock vertical elevation directly onto the terrain heightmap
  rocks.set_elevation_from_terrain(z);

  std::cout << "=== Final Rock Field (Sorted & Packed) ===\n";
  std::cout << rocks.to_string() << "\n\n";

  // --- 7. Query and Export

  // filter by class ID
  auto class_ids = rocks.get_class_ids();
  std::cout << "Found " << class_ids.size() << " rock classes.\n";

  hmap::RockField boulders = rocks.filter_by_class(0u);
  hmap::RockField pebbles = rocks.filter_by_class(1u);
  std::cout << "  - Boulders (class 0): " << boulders.size() << " items\n";
  std::cout << "  - Pebbles  (class 1): " << pebbles.size() << " items\n";

  // --- 8. Export

  rocks.to_png("rock_density.png", shape, density);

  std::cout << "\nExported rock_density.png and rock_field.csv successfully.\n";
  return 0;
}
