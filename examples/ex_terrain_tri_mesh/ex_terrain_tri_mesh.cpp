#include "highmap.hpp"

int main(void)
{
  // --- parameters

  glm::ivec2   shape = {256, 256};
  glm::vec2    kw = {4.f, 4.f};
  int          seed = 0;
  hmap::BBox2D bbox = {0.f, 1.f, 0.f, 1.f};
  size_t       count = 1000;

  // --- generate sample points

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::clamp_min(z0, 0.f);

  hmap::Cloud cloud = hmap::random_cloud(count,
                                         seed,
                                         hmap::PointSamplingMethod::RND_LHS,
                                         bbox);
  cloud.snap_points_to_bounding_box(bbox);
  cloud.set_values_from_array(z0, bbox);

  std::vector<glm::vec3> points = cloud.to_vec3();

  // --- tests TerrainTriMesh

  {
    auto mesh = hmap::TerrainTriMesh(points);
    mesh.relax_xy(0.5f, 10);
    mesh.export_obj("mesh.obj");

    auto areas = mesh.get_vertex_areas(true);

    mesh.to_array(shape).to_png("mesh.png", hmap::Cmap::JET);
    mesh.to_array(shape, areas).to_png("mesh_areas.png", hmap::Cmap::JET);
  }

  // --- optimized remeshing

  {
    float max_error = 1e-2f;
    auto  mesh = hmap::generate_terrain_tri_mesh_from_heightmap(z0, max_error);
    mesh.to_array(shape).to_png("mesh_opt.png", hmap::Cmap::JET);
    mesh.print_info();
    mesh.export_obj("mesh_opt.obj");
  }

  // --- subdivise

  {
    auto mesh = hmap::TerrainTriMesh(points);
    mesh.print_info();
    mesh.subdivise();
    mesh.print_info();
    mesh.to_array(shape).to_png("mesh_subdivise.png", hmap::Cmap::JET);
  }

  // --- slope limiter

  {
    float max_slope = 0.5f;
    int   iterations = 10;
    float sigma = 0.2f;

    auto mesh = hmap::TerrainTriMesh(points);
    mesh.slope_limiter(max_slope, iterations, sigma);
    mesh.to_array(shape).to_png("mesh_slope_limiter.png", hmap::Cmap::JET);
    mesh.export_obj("mesh_slope_limiter.obj");
  }
}
