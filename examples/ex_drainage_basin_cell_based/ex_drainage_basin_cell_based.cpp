#include "highmap.hpp"
#include "highmap/internal/timer.hpp"

int main(void)
{
  hmap::gpu::init_opencl();

  glm::ivec2 shape = {256, 256};
  // shape = {1024, 1024};

  glm::vec2 kw = {4.f, 4.f};
  int       seed = 1;

  // --- generate base mesh and elevation

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                   shape,
                                   kw,
                                   seed,
                                   8,
                                   0.7f);
  z0 = hmap::bulkify(z0, hmap::PrimitiveType::PRIM_CONE_SMOOTH, 2.f);
  hmap::remap(z0);

  // --- drainage basin

  float noise_strength = 0.2f;

  auto db = hmap::DrainageBasinCellBased(z0);

  // db.update_stream_tree(seed, noise_strength);

  hmap::Array erodibility = 1.f - z0; //(shape, 1.f);
  hmap::Array max_slope = z0;
  hmap::remap(max_slope, 1.f / shape.x, 4.f / shape.x);

  float m_exp = 0.5f;
  float uplift_rate = 1.f;

  for (int it = 0; it < 100; ++it)
  {
    db.update_stream_tree(seed, noise_strength);

    hmap::Array acc(shape);
    db.accumulate_area_by_outlet(acc);

    auto  response_times = db.compute_response_times(acc, erodibility, m_exp);
    float diff = db.update_elevations(response_times, uplift_rate, max_slope);

    hmap::log::trace("{}, {}", it, diff);
  }

  db.get_z().dump("out1.png");

  hmap::Array acc(shape, 0.f);
  db.accumulate_area_by_outlet(acc);
  acc.dump("acc.png");

  // hmap::Cloud cloud = hmap::random_cloud(count,
  //                                        seed,
  //                                        hmap::PointSamplingMethod::RND_LHS,
  //                                        bbox);
  // cloud.snap_points_to_bounding_box(bbox);
  // cloud.set_values_from_array(z0, bbox);

  // // --- drainage basin

  // auto db = hmap::DrainageBasin(cloud.to_vec3());

  // db.update_stream_tree();
  // db.to_csv("db.csv");

  // // --- export

  // std::vector<float> area = db.compute_vertex_areas();
  // std::vector<float> flow(db.size(), 1.f);
  // db.accumulate_area_by_outlet(area, flow);

  // hmap::Array zm = db.get_mesh().to_array(shape);
  // hmap::Array a_area = db.get_mesh().to_array(shape, area);
  // hmap::Array a_flow = db.get_mesh().to_array(shape, flow);

  // hmap::remap(a_area);
  // hmap::remap(a_flow);

  hmap::export_banner_png("ex_drainage_basin_cell_based.png",
                          {z0},
                          hmap::Cmap::INFERNO);
}
