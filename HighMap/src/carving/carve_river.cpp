/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/carving.hpp"
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
                 float         outlet_width,
                 float         headwater_width,
                 float         width_exponent,
                 float         depth,
                 float         bank_width_ratio,
                 RadialProfile radial_profile,
                 float         radial_profile_parameter,
                 const Array  *p_noise_r,
                 Array        *p_river_mask,
                 glm::vec4     bbox)
{
  if (!validate_non_empty(z)) return;
  if (p_noise_r && !validate_same_shape(z, *p_noise_r)) return;
  if (graph.empty() || graph.num_edges() == 0) return;

  const size_t n = graph.size();
  const auto  &edges = graph.get_edges();

  // compute flow accumulation along the graph
  std::vector<float> accum = graph.compute_flow_accumulation();

  float max_accum = 0.f;
  for (float a : accum)
    if (a > max_accum) max_accum = a;

  // compute per-node bed widths and influence widths
  float min_w = std::max(0.f, headwater_width);
  float max_w = std::max(min_w, outlet_width);

  if (max_w <= 0.f) return;

  std::vector<float> node_bed_widths(n, min_w);
  std::vector<float> node_influence_widths(n, min_w);

  float bank_ratio = std::max(0.f, bank_width_ratio);

  for (size_t i = 0; i < n; ++i)
  {
    float norm = (max_accum > 1e-6f) ? (accum[i] / max_accum) : 0.f;
    norm = std::clamp(norm, 0.f, 1.f);

    float factor = std::pow(norm, width_exponent);
    float w_bed = min_w + (max_w - min_w) * factor;

    node_bed_widths[i] = w_bed;
    node_influence_widths[i] = w_bed * (1.f + bank_ratio);
  }

  float max_influence_width = max_w * (1.f + bank_ratio);

  if (p_noise_r)
  {
    max_influence_width = std::max(0.f,
                                   max_influence_width *
                                       (1.f + p_noise_r->max()));
  }

  if (max_influence_width <= 0.f) return;

  const glm::ivec2 &shape = z.shape;
  auto              profile_fct = get_radial_profile_function(radial_profile,
                                                 radial_profile_parameter);

  float dx = (bbox.y - bbox.x) / float(shape.x);
  float dy = (bbox.w - bbox.z) / float(shape.y);

  Array river_mask;
  if (p_river_mask) river_mask = Array(shape, 0.f);

  // build oriented downstream segments
  std::vector<Segment2D> all_segments;
  all_segments.reserve(edges.size());

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
    seg.v0 = p0.v - depth;
    seg.v1 = p1.v - depth;
    seg.s0 = node_influence_widths[from];
    seg.s1 = node_influence_widths[to];
    seg.k0 = static_cast<size_t>(from);
    seg.k1 = static_cast<size_t>(to);
    seg.path_index = 0;
    seg.id = static_cast<uint32_t>(all_segments.size());
    seg.d = seg.p1 - seg.p0;

    float len_sq = seg.d.x * seg.d.x + seg.d.y * seg.d.y;
    seg.inv_len_sq = (len_sq > 1e-12f) ? (1.f / len_sq) : 0.f;

    all_segments.push_back(seg);
  }

  if (all_segments.empty()) return;

  // build BVH over segments
  BVH2D bvh(all_segments, max_influence_width);

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

  std::vector<float> xg, yg;
  grid_xy_vector(xg, yg, shape, bbox, false);

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

        float min_carved_z = orig_z;
        float max_river_weight = 0.f;

        float dr = p_noise_r ? (*p_noise_r)(i, j) : 0.f;

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

          float w_influence = (1.f - t) * seg.s0 + t * seg.s1;
          w_influence = std::max(0.f, w_influence * (1.f + dr));

          if (w_influence <= 0.f || d_sq >= w_influence * w_influence) continue;

          float r = std::sqrt(d_sq) / w_influence;
          float profile_t = profile_fct(r);

          float z_floor = (1.f - t) * seg.v0 + t * seg.v1;
          float candidate_z = lerp(z_floor, orig_z, profile_t);

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
