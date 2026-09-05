#include <chrono>
#include <iostream>

#include "highmap.hpp"

int main(void)
{
  const glm::ivec2 shape = {512, 512};
  const uint32_t   seed = 1;

  const float water_level = 0.5f;
  const float coast_margin = 0.1f;
  const float road_level = 0.65f;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  shape,
                                  {4.f, 4.f},
                                  seed);
  hmap::remap(z);

  // 1. Shape constraint: a five-pointed star island. Vertices inside the
  // mask must be above the water level, vertices outside must be below (with
  // a small margin on each side to get a clear coastline).

  hmap::Path star;
  for (int k = 0; k < 10; ++k)
  {
    float angle = 0.5f * M_PI + 2.f * M_PI * static_cast<float>(k) / 10.f;
    float r = (k % 2 == 0) ? 0.42f : 0.18f;
    star.add_point(
        {0.5f + r * std::cos(angle), 0.5f + r * std::sin(angle), 1.f});
  }
  star.set_closed(true);

  hmap::Array mask(shape, 0.f);
  star.to_array(mask, {0.f, 1.f, 0.f, 1.f}, true);

  hmap::Array sea(shape, 1.f);
  sea -= mask;

  hmap::DeformationConstraint land;
  land.target = hmap::Array(shape, water_level + coast_margin);
  land.weight = mask;
  land.type = hmap::DeformationConstraintType::ABOVE;

  hmap::DeformationConstraint ocean;
  ocean.target = hmap::Array(shape, water_level - coast_margin);
  ocean.weight = sea;
  ocean.type = hmap::DeformationConstraintType::BELOW;

  // 2. Path constraint: a flat S-shaped road at a fixed elevation, the rest
  // of the terrain being preserved as much as possible (to avoid a global
  // flattening of the terrain)

  hmap::Path road;
  road.add_point({0.08f, 0.15f, 1.f});
  road.add_point({0.35f, 0.25f, 1.f});
  road.add_point({0.50f, 0.50f, 1.f});
  road.add_point({0.65f, 0.75f, 1.f});
  road.add_point({0.92f, 0.85f, 1.f});
  road = hmap::catmullrom(road);

  hmap::Array path(shape, 0.f);
  road.to_array(path);
  path = hmap::dilation(path, 4);

  hmap::Array not_path(shape, 1.f);
  not_path -= path;

  hmap::DeformationConstraint flat_road;
  flat_road.target = hmap::Array(shape, road_level);
  flat_road.weight = path;
  flat_road.type = hmap::DeformationConstraintType::MATCH;

  hmap::DeformationConstraint preserve;
  preserve.target = z;
  preserve.weight = not_path;
  preserve.type = hmap::DeformationConstraintType::MATCH;

  // 3. Run the search for each set of constraints, and for the combination
  // (the road is allowed to cross the sea as a causeway: sea vertices under
  // the road are released from the ocean constraint)

  auto run = [&](const std::string                              &label,
                 const std::vector<hmap::DeformationConstraint> &constraints,
                 int                                             iterations)
  {
    std::vector<hmap::GaussianPush> pushes;

    auto        t0 = std::chrono::steady_clock::now();
    hmap::Array out = hmap::sls_deformation(z,
                                            constraints,
                                            seed,
                                            iterations,
                                            4,
                                            64,
                                            7,
                                            0.f,
                                            512,
                                            0.65f,
                                            0.1f,
                                            1e-3f,
                                            &pushes);
    auto        t1 = std::chrono::steady_clock::now();

    std::cout << label << ": " << pushes.size() << " pushes in "
              << std::chrono::duration<float>(t1 - t0).count() << " s\n";
    return out;
  };

  hmap::Array z_star = run("star island", {land, ocean}, 600);
  hmap::Array z_road = run("flat road", {flat_road, preserve}, 800);

  hmap::DeformationConstraint ocean_causeway = ocean;
  ocean_causeway.weight *= not_path;

  hmap::Array z_both = run("star island + road",
                           {land, ocean_causeway, flat_road},
                           1000);

  // render with the sea flattened at the water level
  for (hmap::Array *p : {&z_star, &z_road, &z_both})
    hmap::clamp_min(*p, water_level);

  hmap::export_banner_png("ex_sls_deformation.png",
                          {z, z_star, z_road, z_both},
                          hmap::Cmap::TERRAIN,
                          true);

  return 0;
}
