#include <iostream>

#include "highmap.hpp"

int main(void)
{
  // --- parameters

  glm::ivec2   shape = {1024, 1024};
  glm::vec2    kw = {4.f, 4.f};
  int          seed = 0;
  hmap::BBox2D bbox = {0.f, 1.f, 0.f, 1.f};
  size_t       count = 50000;

  // --- generate sample points

  hmap::Array z0 = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, kw, seed);
  hmap::remap(z0);

  hmap::Cloud cloud = hmap::random_cloud(count,
                                         seed,
                                         hmap::PointSamplingMethod::RND_LHS,
                                         bbox);
  cloud.snap_points_to_bounding_box(bbox);
  cloud.set_values_from_array(z0, bbox);

  std::vector<glm::vec3> points = cloud.to_vec3();

  // --- tests TerrainTriMesh

  auto mesh = hmap::TerrainTriMesh(points);
  mesh.relax_xy(0.5f, 10);
  mesh.to_array(shape).to_png("mesh0.png", hmap::Cmap::JET, true);

  if (false)
  {
    std::vector<size_t> path = mesh.shortest_path(0, count - 1, true);
    auto               &pts = mesh.get_points();

    for (const auto &p : path)
      pts[p].z -= 0.1f;

    mesh.to_array(shape).to_png("mesh1.png", hmap::Cmap::JET);
  }

  if (false)
  {
    hmap::TerrainTriMesh::ShortestPathResult results =
        mesh.compute_shortest_paths_to_hull(true);

    auto &pts = mesh.get_points();

    for (const auto &p : std::vector<size_t>{0, 10, 50, 100, 150})
    {
      std::vector<size_t> path = mesh.path_to_hull(p, results);

      for (const auto &p : path)
        pts[p].z -= 0.1f;
    }

    mesh.to_array(shape).to_png("mesh2.png", hmap::Cmap::JET);
  }
}
