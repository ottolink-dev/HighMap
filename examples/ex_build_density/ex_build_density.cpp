#include "highmap.hpp"

int main(void)
{
  glm::ivec2 shape = {256, 256};
  shape = {1024, 1024};
  glm::vec2 kw = {4.f, 4.f};
  int       seed = 1;

  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::SIMPLEX2, shape, kw, seed);
  z = hmap::bulkify(z, hmap::PrimitiveType::PRIM_CUBIC_PULSE, 2.f);
  hmap::remap(z);

  // build criteria
  auto gn = hmap::gradient_norm(z);

  auto cz = 1.f - hmap::threshold_smooth(z, 0.f, 0.9f);
  auto cg = 1.f - hmap::threshold_smooth(gn, 3.f / shape.x, 5.f / shape.x);
  auto ca = hmap::select_angle(z, 30.f, 90.f);

  auto cw = hmap::topographic_wetness_index(z);
  hmap::remap(cw);
  hmap::saturate_percentile(cw, 0.f, 0.95f);

  // auto cq = hmap::curvature_quadric(z, 1, hmap::CurvatureType::CT_MEAN);
  // hmap::clamp_min(cq, 0.f);
  // hmap::remap(cq);
  // hmap::saturate_percentile(cq, 0.f, 0.98f);

  // combine
  std::vector<const hmap::Array *> vec = {&cz, &cg, &ca, &cw};
  std::vector<float>               w = {1.f, 1.f, 0.5f, 0.5f};

  auto dl = hmap::build_density_linear(vec, w);
  auto dg = hmap::build_density_log(vec, w);

  z.dump("out0.png");
  dl.dump("out1.png");
  dg.dump("out2.png");

  size_t      count = 20000;
  hmap::Cloud cloud = hmap::random_cloud_inverse_sampling(count, dl, seed);
  cloud.set_values_from_array(z);
  cloud.to_png("cloud.png", hmap::Cmap::INFERNO);

  hmap::export_usd("scene.usdc",
                   z,
                   {},
                   {cloud},
                   {},
                   hmap::MeshType::TRI,
                   0.15f);

  hmap::export_banner_png("ex_build_density.png", {z, dl, dg}, hmap::Cmap::JET);
}
