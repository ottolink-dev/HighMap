/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/carving.hpp"
#include "highmap/filters.hpp"
#include "highmap/geometry/bvh.hpp"
#include "highmap/geometry/graph.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/math/profiles.hpp"

namespace hmap
{

void carve_river(Array        &z,
                 const Graph  &graph,
                 float         talus,
                 float         outlet_width,
                 float         headwater_width,
                 float         width_exponent,
                 float         depth,
                 float         depth_exponent,
                 std::uint32_t seed,
                 float         talus_noise_ratio,
                 RadialProfile radial_profile,
                 float         radial_profile_parameter,
                 Array        *p_river_mask,
                 glm::vec4     bbox)
{
  if (!validate_non_empty(z)) return;
  if (graph.empty() || graph.num_edges() == 0) return;

  const size_t n = graph.size();
  const auto  &edges = graph.get_edges();

  // compute flow accumulation along the graph
  std::vector<float> accum = graph.compute_flow_accumulation();

  float max_accum = 0.f;
  for (float a : accum)
    if (a > max_accum) max_accum = a;

  // compute per-node bed widths and depths
  float min_w = std::max(0.f, headwater_width);
  float max_w = std::max(min_w, outlet_width);

  if (max_w <= 0.f) return;

  std::vector<float> node_widths(n, min_w);
  std::vector<float> node_depths(n, depth);

  for (size_t i = 0; i < n; ++i)
  {
    float norm = (max_accum > 1e-6f) ? (accum[i] / max_accum) : 0.f;
    norm = std::clamp(norm, 0.f, 1.f);

    float w_factor = std::pow(norm, width_exponent);
    node_widths[i] = min_w + (max_w - min_w) * w_factor;

    float d_factor = (depth_exponent > 0.f) ? std::pow(norm, depth_exponent)
                                            : 1.f;
    node_depths[i] = depth * std::max(0.1f, d_factor);
  }

  const glm::ivec2 &shape = z.shape;
  float             dx = (bbox.y - bbox.x) / float(shape.x);
  float             dy = (bbox.w - bbox.z) / float(shape.y);

  // --- 1. Rasterize river bottom elevation and seed mask
  Array z_bed(shape, 0.f);
  Array mask(shape, 0.f);

  std::vector<Segment2D> all_segments;
  all_segments.reserve(edges.size());

  float lx = bbox.y - bbox.x;
  float ly = bbox.w - bbox.z;
  float ppu = std::max(shape.x / lx, shape.y / ly);

  float ai = (shape.x - 1) / lx;
  float bi = -bbox.x * (shape.x - 1) / lx;
  float aj = (shape.y - 1) / ly;
  float bj = -bbox.z * (shape.y - 1) / ly;

  for (const auto &edge : edges)
  {
    int u = edge.u;
    int v = edge.v;

    if (u < 0 || static_cast<size_t>(u) >= n || v < 0 ||
        static_cast<size_t>(v) >= n || u == v)
      continue;

    float vu = graph[u].v;
    float vv = graph[v].v;

    int from = u;
    int to = v;

    if (vu < vv || (vu == vv && u > v))
    {
      from = v;
      to = u;
    }

    const auto &p0 = graph[from];
    const auto &p1 = graph[to];

    Segment2D seg;
    seg.p0 = glm::vec2(p0.x, p0.y);
    seg.p1 = glm::vec2(p1.x, p1.y);
    seg.v0 = p0.v - node_depths[from];
    seg.v1 = p1.v - node_depths[to];
    seg.s0 = node_widths[from];
    seg.s1 = node_widths[to];
    seg.k0 = static_cast<size_t>(from);
    seg.k1 = static_cast<size_t>(to);
    seg.path_index = 0;
    seg.id = static_cast<uint32_t>(all_segments.size());
    seg.d = seg.p1 - seg.p0;

    float len_sq = seg.d.x * seg.d.x + seg.d.y * seg.d.y;
    seg.inv_len_sq = (len_sq > 1e-12f) ? (1.f / len_sq) : 0.f;

    all_segments.push_back(seg);

    // Rasterize segment into z_bed and mask
    float seg_len = std::sqrt(len_sq);
    int   npixels = static_cast<int>(std::ceil(seg_len * ppu)) + 1;

    for (int k = 0; k < npixels; ++k)
    {
      float t = (npixels > 1)
                    ? (static_cast<float>(k) / static_cast<float>(npixels - 1))
                    : 0.f;
      float px = seg.p0.x + t * seg.d.x;
      float py = seg.p0.y + t * seg.d.y;
      float pz = (1.f - t) * seg.v0 + t * seg.v1;

      int pi = static_cast<int>(std::round(ai * px + bi));
      int pj = static_cast<int>(std::round(aj * py + bj));

      if (pi >= 0 && pi < shape.x && pj >= 0 && pj < shape.y)
      {
        if (mask(pi, pj) == 0.f)
        {
          z_bed(pi, pj) = pz;
          mask(pi, pj) = 1.f;
        }
        else
        {
          z_bed(pi, pj) = std::min(z_bed(pi, pj), pz);
        }
      }
    }
  }

  if (all_segments.empty()) return;

  // --- 2. Expand talus upward from riverbed seeds
  int ir = 2;
  expand_talus(z_bed, mask, talus, seed, ir, talus_noise_ratio);

  // --- 3. Build BVH for exact network SDF and per-segment width interpolation
  BVH2D bvh(all_segments, max_w);

  AABB2D domain_aabb(bbox.x, bbox.z, bbox.y, bbox.w);
  if (!bvh.root_bbox().overlaps(domain_aabb)) return;

  int tile_size = 32;
  int num_tiles_x = (shape.x + tile_size - 1) / tile_size;
  int num_tiles_y = (shape.y + tile_size - 1) / tile_size;
  int total_tiles = num_tiles_x * num_tiles_y;

  std::vector<std::vector<uint32_t>> tile_segments(total_tiles);

#pragma omp parallel for schedule(static)
  for (int tile_idx = 0; tile_idx < total_tiles; ++tile_idx)
  {
    int tx = tile_idx % num_tiles_x;
    int ty = tile_idx / num_tiles_x;

    float tile_xmin = bbox.x + float(tx * tile_size) * dx;
    float tile_xmax = bbox.x +
                      float(std::min((tx + 1) * tile_size, shape.x)) * dx;
    float tile_ymin = bbox.z + float(ty * tile_size) * dy;
    float tile_ymax = bbox.z +
                      float(std::min((ty + 1) * tile_size, shape.y)) * dy;

    AABB2D tile_box(tile_xmin, tile_ymin, tile_xmax, tile_ymax);
    bvh.query_box(tile_box, tile_segments[tile_idx]);
  }

  auto profile_fct = get_radial_profile_function(radial_profile,
                                                 radial_profile_parameter);

  Array river_mask;
  if (p_river_mask) river_mask = Array(shape, 0.f);

  std::vector<float> xg, yg;
  grid_xy_vector(xg, yg, shape, bbox, false);

  // --- 4. Merge z_bed with original z using SDF and local width scaling
#pragma omp parallel for schedule(dynamic, 8)
  for (int tile_idx = 0; tile_idx < total_tiles; ++tile_idx)
  {
    const auto &segs = tile_segments[tile_idx];
    if (segs.empty()) continue;

    int tx = tile_idx % num_tiles_x;
    int ty = tile_idx / num_tiles_x;

    int i_start = tx * tile_size;
    int i_end = std::min(i_start + tile_size, shape.x);
    int j_start = ty * tile_size;
    int j_end = std::min(j_start + tile_size, shape.y);

    for (int j = j_start; j < j_end; ++j)
    {
      float y = yg[j];

      for (int i = i_start; i < i_end; ++i)
      {
        float x = xg[i];
        float orig_z = z(i, j);
        float valley_z = z_bed(i, j);

        float min_carved_z = orig_z;
        float max_river_weight = 0.f;

        for (uint32_t seg_idx : segs)
        {
          const auto &seg = all_segments[seg_idx];

          float qx = x - seg.p0.x;
          float qy = y - seg.p0.y;

          float t = std::clamp((qx * seg.d.x + qy * seg.d.y) * seg.inv_len_sq,
                               0.f,
                               1.f);

          float proj_x = seg.p0.x + t * seg.d.x;
          float proj_y = seg.p0.y + t * seg.d.y;

          float d_sq = (x - proj_x) * (x - proj_x) +
                       (y - proj_y) * (y - proj_y);

          float w_local = (1.f - t) * seg.s0 + t * seg.s1;

          if (w_local <= 0.f || d_sq >= w_local * w_local) continue;

          float dist = std::sqrt(d_sq);
          float r = dist / w_local;
          float profile_t = profile_fct(r);

          float candidate_z = lerp(valley_z, orig_z, profile_t);

          if (candidate_z < min_carved_z) min_carved_z = candidate_z;

          float weight = 1.f - profile_t;
          if (weight > max_river_weight) max_river_weight = weight;
        }

        z(i, j) = min_carved_z;

        if (p_river_mask) river_mask(i, j) = max_river_weight;
      }
    }
  }

  if (p_river_mask) *p_river_mask = std::move(river_mask);
}

} // namespace hmap
