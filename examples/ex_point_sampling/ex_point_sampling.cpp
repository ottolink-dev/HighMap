#include <iostream>

#include "highmap.hpp"

int main(void)
{
  // for render only
  glm::ivec2               shape = {256, 256};
  hmap::Array              raster(shape);
  std::vector<hmap::Array> zs = {};
  std::uint32_t            seed = 0;

  // density field
  glm::vec2   kw = {2.f, 2.f};
  hmap::Array density = hmap::noise(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(density);   // /!\ NEEDS TO BE IN [0, 1]
  hmap::gain(density, 2); // sharper transition

  // base param
  size_t count = 1000;

  // --- random

  for (auto &type : {hmap::PointSamplingMethod::RND_RANDOM,
                     hmap::PointSamplingMethod::RND_HALTON,
                     hmap::PointSamplingMethod::RND_HAMMERSLEY,
                     hmap::PointSamplingMethod::RND_LHS})
  {
    hmap::Cloud cloud = hmap::random_cloud(count, seed, type);

    cloud.to_array(raster);
    zs.push_back(raster);
  }

  hmap::export_banner_png("ex_point_sampling0.png", zs, hmap::Cmap::BONE);

  // --- density

  {
    zs.clear();
    zs.push_back(density);
    raster = 0.f;

    // auto        xy = hmap::random_points_density(count, density, seed);
    // hmap::Cloud cloud(xy[0], xy[1], 1.f /* value */);

    // rejection sampling
    hmap::Cloud cloud = hmap::random_cloud_density(count, density, seed);

    cloud.to_array(raster);
    zs.push_back(raster);

    std::cout << "cloud_density count: " << cloud.size() << "\n";

    // inverse transform sampling (batch)
    raster = 0.f;
    hmap::Cloud cloud_inv = hmap::random_cloud_inverse_sampling(count,
                                                                density,
                                                                seed);
    cloud_inv.to_array(raster);
    zs.push_back(raster);

    std::cout << "cloud_inv count: " << cloud_inv.size() << "\n";

    // sequential inverse transform sampling with dynamic local updates (e.g.
    // repelling radius)
    raster = 0.f;
    hmap::InverseSampler2D seq_sampler(density, seed);
    hmap::Cloud            cloud_seq;
    int inhibit_r = 4; // grid radius to depress after each point
    for (size_t k = 0; k < count; ++k)
    {
      glm::vec2 p = seq_sampler.sample();
      cloud_seq.push_back(hmap::Point(p.x, p.y, 1.f));

      // Depress density around the sampled point on the grid
      int cx = static_cast<int>(p.x * (shape.x - 1));
      int cy = static_cast<int>(p.y * (shape.y - 1));
      seq_sampler.update_region(
          [cx, cy, inhibit_r](int i, int j, float cur)
          {
            float d2 = static_cast<float>((i - cx) * (i - cx) +
                                          (j - cy) * (j - cy));
            float r2 = static_cast<float>(inhibit_r * inhibit_r);
            if (d2 < r2)
            {
              float factor = d2 / r2;
              return cur * factor;
            }
            return cur;
          },
          {cx - inhibit_r, cx + inhibit_r},
          {cy - inhibit_r, cy + inhibit_r});
    }

    cloud_seq.to_array(raster);
    zs.push_back(raster);

    std::cout << "cloud_seq count: " << cloud_seq.size() << "\n";

    // filter
    cloud = hmap::random_cloud(count,
                               seed,
                               hmap::PointSamplingMethod::RND_HALTON);

    hmap::rejection_filter_density(cloud, density, seed);

    raster = 0.f;
    cloud.to_array(raster);
    zs.push_back(raster);

    std::cout << "cloud_rejection count: " << cloud.size() << "\n";

    hmap::export_banner_png("ex_point_sampling1.png", zs, hmap::Cmap::BONE);
  }

  // --- distance based

  {
    zs.clear();
    zs.push_back(density);

    float min_dist = 0.02f;
    float max_dist = 0.08f;

    hmap::Cloud cloud = hmap::random_cloud_distance(min_dist, seed);

    raster = 0.f;
    cloud.to_array(raster);
    zs.push_back(raster);

    cloud = hmap::random_cloud_distance(min_dist, max_dist, density, seed);

    std::cout << "count: " << cloud.size() << "\n";

    raster = 0.f;
    cloud.to_array(raster);
    zs.push_back(raster);

    cloud = hmap::random_cloud_distance_power_law(0.01f /* dist_min */,
                                                  0.2f /* dist_max */,
                                                  1.2f /* alpha */,
                                                  seed);

    std::cout << "count: " << cloud.size() << "\n";

    raster = 0.f;
    cloud.to_array(raster);
    zs.push_back(raster);

    cloud = hmap::random_cloud_distance_weibull(0.01f /* dist_min */,
                                                0.05f /* lambda */,
                                                1.f /* k */,
                                                seed);

    std::cout << "count: " << cloud.size() << "\n";

    raster = 0.f;
    cloud.to_array(raster);
    zs.push_back(raster);

    hmap::export_banner_png("ex_point_sampling2.png", zs, hmap::Cmap::BONE);
  }

  // --- grid jittered

  {
    zs.clear();

    glm::vec2 jitter_amount = {0.3f, 0.3f};
    glm::vec2 stagger_ratio = {0.5f, 0.f};

    hmap::Cloud cloud = hmap::random_cloud_jittered(count,
                                                    jitter_amount,
                                                    stagger_ratio,
                                                    seed);
    cloud.set_values_from_min_distance();

    raster = 0.f;
    cloud.to_array(raster);
    zs.push_back(raster);

    hmap::export_banner_png("ex_point_sampling3.png", zs, hmap::Cmap::BONE);
  }
}
