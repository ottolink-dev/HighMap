#include <vector>

#include "highmap.hpp"
#include "highmap/dbg/timer.hpp"

int main(void)
{
  const glm::ivec2 shape = {512, 512};
  const glm::vec2  res = {2.f, 2.f};
  const int        seed = 42;

  // generate terrain cost map
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN, shape, res, seed);

  const glm::ivec2 ij_start = {50, 50};
  const glm::ivec2 ij_end = {460, 460};
  const glm::vec4  bbox = {0.f, 1.f, 0.f, 1.f};

  const float elevation_ratio = 0.1f;
  const float distance_exponent = 2.f;
  const float upward_penalization = 0.1f;

  // --- 1. find_path_dijkstra (full-resolution Dijkstra)

  hmap::Timer::Start("find_path_dijkstra");
  hmap::Path path_dijk = hmap::find_path_dijkstra(z,
                                                  ij_start,
                                                  ij_end,
                                                  bbox,
                                                  elevation_ratio,
                                                  distance_exponent,
                                                  upward_penalization);
  hmap::Timer::Stop("find_path_dijkstra");

  hmap::Array w_dijk = path_dijk.to_array(shape);

  // --- 2. find_path_midpoint (midpoint displacement heuristic)

  const float offset_ratio = 0.2f;

  hmap::Timer::Start("find_path_midpoint");
  hmap::Path path_midp = hmap::find_path_midpoint(z,
                                                  ij_start,
                                                  ij_end,
                                                  bbox,
                                                  offset_ratio);
  hmap::Timer::Stop("find_path_midpoint");

  hmap::Array w_midp = path_midp.to_array(shape);

  // --- 3. find_path_multiscale (coarse-to-fine shortest path)

  const int   n_levels = 4;
  const int   corridor_radius = 12;
  const float corridor_decay = 1.f;
  const bool  use_astar = false;
  const bool  smooth_path = false;

  hmap::Timer::Start("find_path_multiscale");
  hmap::Path path_multi = hmap::find_path_multiscale(z,
                                                     ij_start,
                                                     ij_end,
                                                     bbox,
                                                     n_levels,
                                                     corridor_radius,
                                                     corridor_decay,
                                                     elevation_ratio,
                                                     distance_exponent,
                                                     upward_penalization,
                                                     nullptr,
                                                     use_astar,
                                                     smooth_path);
  hmap::Timer::Stop("find_path_multiscale");

  hmap::Array w_multi = path_multi.to_array(shape);

  // --- Export banner

  hmap::export_banner_png("ex_find_path_multiscale.png",
                          {z, w_dijk, w_midp, w_multi},
                          hmap::Cmap::TERRAIN);

  hmap::Timer::Dump();

  return 0;
}
