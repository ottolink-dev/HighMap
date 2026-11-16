#include "highmap.hpp"
#include "highmap/internal/timer.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};

  glm::vec2 kw = {4.f, 4.f};
  int       seed = 1;
  glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f};
  size_t    count = 1000;

  // --- generate base mesh and elevation

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                   shape,
                                   kw,
                                   seed,
                                   8,
                                   0.7f);
  z0 = hmap::bulkify(z0, hmap::PrimitiveType::PRIM_CONE_SMOOTH, 2.f);
  hmap::remap(z0);

  hmap::Cloud cloud = hmap::random_cloud(count,
                                         seed,
                                         hmap::PointSamplingMethod::RND_LHS,
                                         bbox);
  cloud.snap_points_to_bounding_box(bbox);
  cloud.set_values_from_array(z0, bbox);

  // --- drainage basin

  auto db = hmap::DrainageBasin(cloud.to_vec3());

  db.update_stream_tree();
  db.to_csv("db.csv");

  // --- main channels

  auto mc = db.get_main_channels();

  if (false)
  {
    for (const auto &channel : mc)
      for (const auto &i : channel)
        db.get_mesh().get_points()[i].z = 0.f;
  }

  // --- export

  std::vector<float> area = db.compute_vertex_areas();
  std::vector<float> flow(db.size(), 1.f);
  db.accumulate_area_by_outlet(area, flow);

  hmap::Array zm = db.get_mesh().to_array(shape);
  hmap::Array a_area = db.get_mesh().to_array(shape, area);
  hmap::Array a_flow = db.get_mesh().to_array(shape, flow);

  hmap::remap(a_area);
  hmap::remap(a_flow);

  hmap::export_banner_png("ex_drainage_basin.png",
                          {z0, zm, a_area, a_flow},
                          hmap::Cmap::INFERNO);
}
