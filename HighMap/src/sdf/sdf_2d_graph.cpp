/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/bvh.hpp"
#include "highmap/geometry/graph.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/sdf.hpp"

namespace hmap
{

namespace
{

inline float aabb_distance_sq(const glm::vec2 &p, const AABB2D &box)
{
  float dx = std::max(0.f, std::max(box.min.x - p.x, p.x - box.max.x));
  float dy = std::max(0.f, std::max(box.min.y - p.y, p.y - box.max.y));
  return dx * dx + dy * dy;
}

inline float segment_distance_sq(const glm::vec2 &p, const Segment2D &seg)
{
  glm::vec2 w = p - seg.p0;
  float     t = std::clamp(glm::dot(w, seg.d) * seg.inv_len_sq, 0.f, 1.f);
  glm::vec2 b = w - seg.d * t;
  return glm::dot(b, b);
}

float query_nearest_dist_sq_bvh(const BVH2D &bvh, const glm::vec2 &p)
{
  const auto &nodes = bvh.get_nodes();
  const auto &segs = bvh.get_segments();
  if (nodes.empty()) return std::numeric_limits<float>::max();

  float best_dist_sq = std::numeric_limits<float>::max();

  struct StackEntry
  {
    int32_t node_idx;
    float   dist_sq;
  };
  StackEntry stack[64];
  int        stack_ptr = 0;

  stack[stack_ptr++] = {0, aabb_distance_sq(p, nodes[0].bbox)};

  while (stack_ptr > 0)
  {
    StackEntry cur = stack[--stack_ptr];
    if (cur.dist_sq >= best_dist_sq) continue;

    const auto &node = nodes[cur.node_idx];
    if (node.is_leaf())
    {
      for (uint32_t i = 0; i < node.seg_count; ++i)
      {
        float d2 = segment_distance_sq(p, segs[node.seg_start + i]);
        if (d2 < best_dist_sq) best_dist_sq = d2;
      }
    }
    else
    {
      float d_left = aabb_distance_sq(p, nodes[node.left].bbox);
      float d_right = aabb_distance_sq(p, nodes[node.right].bbox);

      if (d_left < d_right)
      {
        if (d_right < best_dist_sq) stack[stack_ptr++] = {node.right, d_right};
        if (d_left < best_dist_sq) stack[stack_ptr++] = {node.left, d_left};
      }
      else
      {
        if (d_left < best_dist_sq) stack[stack_ptr++] = {node.left, d_left};
        if (d_right < best_dist_sq) stack[stack_ptr++] = {node.right, d_right};
      }
    }
  }

  return best_dist_sq;
}

} // namespace

Array sdf_2d_graph(const Graph &graph,
                   glm::ivec2   shape,
                   glm::vec4    bbox,
                   const Array *p_noise_x,
                   const Array *p_noise_y,
                   bool         use_bvh)
{
  if (use_bvh)
    return sdf_2d_graph_bvh(graph, shape, bbox, p_noise_x, p_noise_y);

  return sdf_2d_graph_brute_force(graph, shape, bbox, p_noise_x, p_noise_y);
}

Array sdf_2d_graph_brute_force(const Graph &graph,
                               glm::ivec2   shape,
                               glm::vec4    bbox,
                               const Array *p_noise_x,
                               const Array *p_noise_y)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  if (graph.empty() || graph.num_edges() == 0)
  {
    return Array(shape, std::numeric_limits<float>::max());
  }

  Array sdf(shape);

  std::vector<float> xp = graph.get_edge_x_pairs();
  std::vector<float> yp = graph.get_edge_y_pairs();

  std::vector<float> xg(shape.x);
  std::vector<float> yg(shape.y);
  grid_xy_vector(xg, yg, shape, bbox, false);

#pragma omp parallel for collapse(2) schedule(static)
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float dx = p_noise_x ? (*p_noise_x)(i, j) : 0.f;
      float dy = p_noise_y ? (*p_noise_y)(i, j) : 0.f;

      float x = xg[i] + dx;
      float y = yg[j] + dy;

      float d = std::numeric_limits<float>::max();

      for (size_t k = 0; k < xp.size() - 1; k += 2)
      {
        size_t    p = k + 1;
        glm::vec2 e = {xp[p] - xp[k], yp[p] - yp[k]};
        glm::vec2 w = {x - xp[k], y - yp[k]};
        float     len2 = glm::dot(e, e);
        float t = (len2 > 1e-12f) ? std::clamp(glm::dot(w, e) / len2, 0.f, 1.f)
                                  : 0.f;
        glm::vec2 b = {w.x - e.x * t, w.y - e.y * t};
        d = std::min(d, glm::dot(b, b));
      }

      sdf(i, j) = std::sqrt(d);
    }

  return sdf;
}

Array sdf_2d_graph_bvh(const Graph &graph,
                       glm::ivec2   shape,
                       glm::vec4    bbox,
                       const Array *p_noise_x,
                       const Array *p_noise_y)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  if (graph.empty() || graph.num_edges() == 0)
  {
    return Array(shape, std::numeric_limits<float>::max());
  }

  const auto            &edges = graph.get_edges();
  std::vector<Segment2D> segs;
  segs.reserve(edges.size());

  for (size_t k = 0; k < edges.size(); ++k)
  {
    int u = edges[k].u;
    int v = edges[k].v;

    if (u < 0 || static_cast<size_t>(u) >= graph.size() || v < 0 ||
        static_cast<size_t>(v) >= graph.size() || u == v)
      continue;

    Segment2D seg;
    seg.p0 = glm::vec2(graph[u].x, graph[u].y);
    seg.p1 = glm::vec2(graph[v].x, graph[v].y);
    seg.d = seg.p1 - seg.p0;
    float lsq = glm::dot(seg.d, seg.d);
    seg.inv_len_sq = (lsq > 0.f) ? (1.f / lsq) : 0.f;
    seg.id = static_cast<uint32_t>(k);
    segs.push_back(seg);
  }

  if (segs.empty())
  {
    return Array(shape, std::numeric_limits<float>::max());
  }

  BVH2D bvh(segs, 0.f);

  Array sdf(shape);

  std::vector<float> xg(shape.x);
  std::vector<float> yg(shape.y);
  grid_xy_vector(xg, yg, shape, bbox, false);

#pragma omp parallel for collapse(2) schedule(static)
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float dx = p_noise_x ? (*p_noise_x)(i, j) : 0.f;
      float dy = p_noise_y ? (*p_noise_y)(i, j) : 0.f;

      float x = xg[i] + dx;
      float y = yg[j] + dy;

      float dist_sq = query_nearest_dist_sq_bvh(bvh, glm::vec2(x, y));
      sdf(i, j) = std::sqrt(dist_sq);
    }

  return sdf;
}

} // namespace hmap
