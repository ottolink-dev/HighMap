/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <queue>
#include <random>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/filters.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate2d.hpp"
#include "highmap/primitives/coherent_noise.hpp"
#include "highmap/terrain_tri_mesh.hpp"

namespace hmap
{

Array diffusion_limited_aggregation(glm::ivec2    shape,
                                    float         scale,
                                    std::uint32_t seed,
                                    float         seeding_radius,
                                    float         seeding_outer_radius_ratio,
                                    float         slope,
                                    float         noise_ratio)
{
  if (!validate_shape(shape)) return Array();

  std::mt19937                          gen(seed);
  std::uniform_real_distribution<float> dis(0.f, 1.f);

  // neighbor search
  std::vector<int> di = {-1, 0, 0, 1, -1, -1, 1, 1};
  std::vector<int> dj = {0, 1, -1, 0, -1, 1, -1, 1};

  // --- work on a grid with a resolution defined by the 'scale'
  int ncells = std::max(1, (int)(1.f / scale));

  glm::ivec2 shape_wrk = {ncells, ncells};
  Array      wrk = Array(shape_wrk);

  int   nwalkers = ncells * ncells;
  float ratio = std::pow(0.01f, 1.f / ncells);

  // seed the diffusion process (center of domain)
  int ic = (int)(0.5f * ncells);
  int jc = (int)(0.5f * ncells);
  wrk(ic, jc) = 1.f;

  for (int k = 0; k < nwalkers; k++)
  {
    // pick a random cell on a circle
    float theta = 2.f * M_PI * dis(gen);
    int   i = (int)(0.5f * shape_wrk.x +
                  seeding_radius *
                      (1.f + seeding_outer_radius_ratio * dis(gen)) *
                      (shape_wrk.x - 1.0) * std::cos(theta));
    int   j = (int)(0.5f * shape_wrk.y +
                  seeding_radius *
                      (1.f + seeding_outer_radius_ratio * dis(gen)) *
                      (shape_wrk.y - 1.0) * std::sin(theta));

    // random_walk
    bool keep_walking = true;

    while (i > 0 && j > 0 && i < shape_wrk.x - 1 && j < shape_wrk.y - 1 &&
           keep_walking)
    {
      // check neighbors for encounter with an already cell
      // touched by diffusion
      for (size_t p = 0; p < di.size(); p++)
        if (wrk(i + di[p], j + dj[p]) > 0.f)
        {
          wrk(i, j) = ratio * wrk(i + di[p], j + dj[p]);
          keep_walking = false;
          break;
        }

      // next step in random direction
      int p = (int)(std::floor(8.f * dis(gen)));

      i += di[p];
      j += dj[p];
    }
  }

  // clean-up, remove spurious values outward the seeding radius
  for (int j = 0; j < shape_wrk.y; j++)
    for (int i = 0; i < shape_wrk.x; i++)
    {
      float dx = (float)(i - ic) / (shape_wrk.x - 1);
      float dy = (float)(j - jc) / (shape_wrk.y - 1);
      float r = std::hypot(dx, dy);

      if (r > 0.95f * seeding_radius) wrk(i, j) = 0.f;
    }

  fill_borders(wrk);
  fill_talus(wrk, slope / (float)wrk.shape.x, seed, noise_ratio);
  extrapolate_borders(wrk, 2, 0.75f);

  // --- generate the output array within the requested shape

  Array z = wrk.resample_to_shape(shape); // output

  return z;
}

Array diffusion_limited_aggregation_trimesh(
    glm::ivec2            shape,
    std::uint32_t         seed,
    size_t                control_points_count,
    glm::vec2             seed_position,
    float                 ratio,
    float                 stop_proba,
    float                 slope,
    InterpolationMethod2D interpolation_method,
    const Array          *p_noise_x,
    const Array          *p_noise_y)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  // --- Generate triangle mesh

  const glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f};

  Cloud cloud = random_cloud_jittered(control_points_count,
                                      {0.5f, 0.5f},
                                      {0.f, 0.f},
                                      seed,
                                      bbox);
  cloud.snap_points_to_bounding_box(bbox);
  cloud.set_values(0.f);

  // initialize seed
  size_t kc = cloud.nearest_point(seed_position);
  cloud[kc].v = 1.f;

  auto        mesh = TerrainTriMesh(cloud.to_vec3());
  const auto &nbrs_data = mesh.get_neighbors();
  auto       &points = mesh.get_points();

  // --- DLA process

  std::vector<size_t> outlets = mesh.get_convex_hull();
  size_t              nwalkers = mesh.size();

  std::mt19937 gen(seed);

  for (size_t k = 0; k < nwalkers; k++)
  {
    // pick a random outlet cell
    int    random_index = std::uniform_int_distribution<int>(0,
                                                          outlets.size() -
                                                              1)(gen);
    size_t i = outlets[random_index];

    // random_walk
    bool keep_walking = true;

    while (keep_walking)
    {
      // check neighbors for encounter with an already cell touched
      // by diffusion
      for (const auto &nb : nbrs_data.adjacency[i])
      {
        size_t n = nb.index;
        float  rd = std::uniform_real_distribution<float>(0.f, 1.f)(gen);

        if (points[n].z > 0.f && rd < stop_proba)
        {
          points[i].z = ratio * points[n].z;
          keep_walking = false;
          break;
        }
      }

      // next step in random direction
      size_t n_neighbors = nbrs_data.adjacency[i].size();
      int    random_nb = std::uniform_int_distribution<int>(0,
                                                         n_neighbors - 1)(gen);
      i = nbrs_data.adjacency[i][random_nb].index;
    }
  }

  // normalize non-zero values to [0, 1]
  float zmin = std::numeric_limits<float>::max();
  for (const auto &p : points)
    if (p.z > 0.f && p.z < zmin) zmin = p.z;

  for (auto &p : points)
    if (p.z > 0.f) p.z = (p.z - zmin) / (1.f - zmin);

  // --- Fill zero-value with a talus

  if (true)
  {
    struct HeapNode
    {
      float  h;
      size_t node;

      bool operator<(const HeapNode &o) const noexcept
      {
        return h < o.h; // max-heap
      }
    };

    std::vector<HeapNode> heap_storage;
    heap_storage.reserve(mesh.size());

    std::priority_queue<HeapNode> heap; // max-heap

    for (size_t k = 0; k < mesh.size(); ++k)
    {
      if (points[k].z > 0.f)
      {
        heap.push({points[k].z, k});
      }
    }

    while (!heap.empty())
    {
      HeapNode top = heap.top();
      heap.pop();

      size_t i = top.node;

      for (const auto &nb : nbrs_data.adjacency[i])
      {
        size_t j = nb.index;

        if (points[j].z == 0.f)
        {
          glm::vec2 pi(points[i].x, points[i].y);
          glm::vec2 pj(points[j].x, points[j].y);
          float     dist = glm::length(pi - pj);

          points[j].z = std::max(0.001f,
                                 points[i].z - dist * slope * points[i].z);

          heap.push({points[j].z, j});
        }
      }
    }
  }

  // --- Interpolate back to an heightmap

  // interpolate
  std::vector<float> xc, yc, zc;
  for (const auto &p : mesh.get_points())
  {
    xc.push_back(p.x);
    yc.push_back(p.y);
    zc.push_back(p.z);
  }

  Array ze = interpolate2d(shape,
                           xc,
                           yc,
                           zc,
                           interpolation_method,
                           p_noise_x,
                           p_noise_y);

  return ze;
}

} // namespace hmap
