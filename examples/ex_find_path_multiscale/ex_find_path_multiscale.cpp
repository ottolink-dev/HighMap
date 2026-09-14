#include <vector>

#include "highmap.hpp"
#include "highmap/dbg/timer.hpp"

int main(void)
{
  glm::ivec2 shape = {512, 512};
  glm::vec2  res = {2.f, 2.f};
  int        seed = 42;

  // generate terrain cost map
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);

  glm::ivec2 ij_start = {50, 50};
  glm::ivec2 ij_end = {460, 460};
  glm::vec4  bbox = {0.f, 1.f, 0.f, 1.f};

  // --- 1. find_path_dijkstra (full-resolution Dijkstra)

  hmap::Timer::Start("find_path_dijkstra");
  std::vector<int> i_dijk, j_dijk;
  hmap::find_path_dijkstra(z,
                           ij_start,
                           ij_end,
                           i_dijk,
                           j_dijk,
                           0.1f,
                           2.f,
                           10.f);
  hmap::Timer::Stop("find_path_dijkstra");

  std::vector<hmap::Point> pts_dijk;
  pts_dijk.reserve(i_dijk.size());

  for (size_t k = 0; k < i_dijk.size(); ++k)
  {
    pts_dijk.emplace_back(float(i_dijk[k]) / float(shape.x - 1),
                          float(j_dijk[k]) / float(shape.y - 1),
                          z(i_dijk[k], j_dijk[k]));
  }

  hmap::Path  path_dijk(pts_dijk);
  hmap::Array w_dijk = path_dijk.to_array(shape);

  // --- 2. find_path_midpoint (midpoint displacement heuristic)

  hmap::Timer::Start("find_path_midpoint");
  std::vector<glm::ivec2> idx_midp = hmap::find_path_midpoint(z,
                                                              ij_start,
                                                              ij_end,
                                                              0.2f);
  hmap::Timer::Stop("find_path_midpoint");

  std::vector<hmap::Point> pts_midp;
  pts_midp.reserve(idx_midp.size());

  for (const auto &p : idx_midp)
  {
    pts_midp.emplace_back(float(p.x) / float(shape.x - 1),
                          float(p.y) / float(shape.y - 1),
                          z(p.x, p.y));
  }

  hmap::Path  path_midp(pts_midp);
  hmap::Array w_midp = path_midp.to_array(shape);

  // --- 3. find_path_multiscale (coarse-to-fine shortest path)

  hmap::Timer::Start("find_path_multiscale");
  hmap::Path path_multi = hmap::find_path_multiscale(
      z,
      ij_start,
      ij_end,
      bbox,
      4,    // n_levels = 4
      12,   // corridor_radius
      1.f,  // corridor_decay
      0.1f, // elevation_ratio
      2.f,  // distance_exponent
      10.f, // upward_penalization
      nullptr,
      false,  // use_astar
      false); // smooth_path
  hmap::Timer::Stop("find_path_multiscale");

  hmap::Array w_multi = path_multi.to_array(shape);

  // --- Export banner

  hmap::export_banner_png("ex_find_path_multiscale.png",
                          {z, w_dijk, w_midp, w_multi},
                          hmap::Cmap::TERRAIN);

  hmap::Timer::Dump();

  return 0;
}
