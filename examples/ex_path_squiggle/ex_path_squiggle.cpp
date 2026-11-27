#include "highmap.hpp"

int main(void)
{
  glm::ivec2    shape = {256, 256};
  hmap::BBox2D  bbox = {0.f, 1.f, 0.f, 1.f};
  std::uint32_t seed = 42;

  // 1. Initial base path (open polyline traversing domain)
  hmap::Path base_path({hmap::Point(0.1f, 0.3f, 0.2f),
                        hmap::Point(0.4f, 0.55f, 0.5f),
                        hmap::Point(0.65f, 0.45f, 0.8f),
                        hmap::Point(0.9f, 0.8f, 1.0f)});

  hmap::Array z_base = base_path.to_array(shape, bbox);

  // 2. Standard procedural squiggle path (alternating orientation)
  hmap::Path p_squiggle =
      hmap::squiggle(base_path, 4, seed, 0.45f, 0, nullptr, nullptr, bbox);

  hmap::Array z_squiggle = p_squiggle.to_array(shape, bbox);

  // 3. Squiggle with spatial weighting (attracted towards top half)
  hmap::Array weights(shape, 1.f);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float y_norm = (float)j / (float)shape.y;
      weights(i, j) = std::exp(4.f * y_norm);
    }

  hmap::Path p_weighted =
      hmap::squiggle(base_path, 4, seed, 0.45f, 0, &weights, nullptr, bbox);

  hmap::Array z_weighted = p_weighted.to_array(shape, bbox);

  // 4. Squiggle with branching tributaries
  std::vector<hmap::Path> branched_paths = hmap::squiggle_branches(base_path,
                                                                   4,
                                                                   seed,
                                                                   0.9f,
                                                                   6,
                                                                   0.4f,
                                                                   0,
                                                                   nullptr,
                                                                   nullptr,
                                                                   bbox);

  hmap::Array z_branches(shape, 0.f);

  for (const auto &p : branched_paths)
  {
    hmap::Array z_p = p.to_array(shape, bbox);
    z_branches = hmap::maximum(z_branches, z_p);
  }

  // Export visual comparison banner
  hmap::export_banner_png("ex_path_squiggle.png",
                          {z_base, z_squiggle, z_weighted, z_branches},
                          hmap::Cmap::INFERNO);

  return 0;
}
