/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/bvh.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/math/core.hpp"
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

Array sdf_2d_polyline(const Path  &path,
                      glm::ivec2   shape,
                      glm::vec4    bbox,
                      const Array *p_noise_x,
                      const Array *p_noise_y,
                      bool         use_bvh)
{
  if (use_bvh)
    return sdf_2d_polyline_bvh(path, shape, bbox, p_noise_x, p_noise_y);

  return sdf_2d_polyline_brute_force(path, shape, bbox, p_noise_x, p_noise_y);
}

Array sdf_2d_polyline_brute_force(const Path  &path,
                                  glm::ivec2   shape,
                                  glm::vec4    bbox,
                                  const Array *p_noise_x,
                                  const Array *p_noise_y)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  if (path.size() < 2)
  {
    hmap::log::error("at least 2 points needed in the Path to compute the SDF");
    return Array(shape);
  }

  Array sdf2(shape);

  // normalize path nodes coord to unit square
  std::vector<float> xp = path.get_x();
  std::vector<float> yp = path.get_y();

  // array cell coordinates
  std::vector<float> xg(shape.x);
  std::vector<float> yg(shape.y);
  grid_xy_vector(xg, yg, shape, bbox, false);

  // compute sdf
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float dx = p_noise_x ? (*p_noise_x)(i, j) : 0.f;
      float dy = p_noise_y ? (*p_noise_y)(i, j) : 0.f;

      float x = xg[i] + dx;
      float y = yg[j] + dy;

      // distance
      float d = std::numeric_limits<float>::max();

      for (size_t k = 0; k < path.size() - 1; k++)
      {
        size_t    p = k + 1;
        glm::vec2 e = {xp[p] - xp[k], yp[p] - yp[k]};
        glm::vec2 w = {x - xp[k], y - yp[k]};
        float     t = std::clamp(dot(w, e) / dot(e, e), 0.f, 1.f);
        glm::vec2 b = {w.x - e.x * t, w.y - e.y * t};
        d = std::min(d, dot(b, b));
      }
      sdf2(i, j) = std::sqrt(d);
    }

  return sdf2;
}

Array sdf_2d_polyline_bvh(const Path  &path,
                          glm::ivec2   shape,
                          glm::vec4    bbox,
                          const Array *p_noise_x,
                          const Array *p_noise_y)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  if (path.size() < 2)
  {
    hmap::log::error("at least 2 points needed in the Path to compute the SDF");
    return Array(shape);
  }

  // build segments
  std::vector<Segment2D> segs;
  segs.reserve(path.size() - 1);

  for (size_t k = 0; k < path.size() - 1; ++k)
  {
    Segment2D seg;
    seg.p0 = glm::vec2(path[k].x, path[k].y);
    seg.p1 = glm::vec2(path[k + 1].x, path[k + 1].y);
    seg.d = seg.p1 - seg.p0;
    float lsq = glm::dot(seg.d, seg.d);
    seg.inv_len_sq = (lsq > 0.f) ? (1.f / lsq) : 0.f;
    seg.id = static_cast<uint32_t>(k);
    segs.push_back(seg);
  }

  BVH2D bvh(segs, 0.f);

  Array sdf2(shape);

  // array cell coordinates
  std::vector<float> xg(shape.x);
  std::vector<float> yg(shape.y);
  grid_xy_vector(xg, yg, shape, bbox, false);

  // compute sdf
#pragma omp parallel for collapse(2) schedule(static)
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float dx = p_noise_x ? (*p_noise_x)(i, j) : 0.f;
      float dy = p_noise_y ? (*p_noise_y)(i, j) : 0.f;

      float x = xg[i] + dx;
      float y = yg[j] + dy;

      float dist_sq = query_nearest_dist_sq_bvh(bvh, glm::vec2(x, y));
      sdf2(i, j) = std::sqrt(dist_sq);
    }

  return sdf2;
}

Array sdf_2d_polyline_bezier(const Path  &path,
                             glm::ivec2   shape,
                             glm::vec4    bbox,
                             const Array *p_noise_x,
                             const Array *p_noise_y)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  if (path.size() < 3)
  {
    hmap::log::error("at least 3 points needed in the Path to compute the SDF");
    return Array(shape);
  }

  Array sdf2(shape);

  // get points, repeat every 3 points to have a continuous Bezier
  // curve and add a point if the input size is even, to ensure start
  // and end points will be taken into account
  std::vector<float> x_input = path.get_x();
  std::vector<float> y_input = path.get_y();
  std::vector<float> xp;
  std::vector<float> yp;

  if (x_input.size() % 2 == 0)
  {
    // not just in the middle to avoid Bezier degeneration
    size_t k = x_input.size();
    float  new_x = lerp(x_input[k - 1], x_input[k - 2], 0.1f);
    float  new_y = lerp(y_input[k - 1], y_input[k - 2], 0.1f);

    x_input.insert(x_input.end() - 1, new_x);
    y_input.insert(y_input.end() - 1, new_y);
  }

  for (size_t k = 0; k < x_input.size() - 2; k += 2)
  {
    xp.push_back(x_input[k]);
    xp.push_back(x_input[k + 1]);
    xp.push_back(x_input[k + 2]);

    yp.push_back(y_input[k]);
    yp.push_back(y_input[k + 1]);
    yp.push_back(y_input[k + 2]);
  }

  // array cell coordinates
  std::vector<float> xg(shape.x);
  std::vector<float> yg(shape.y);
  grid_xy_vector(xg, yg, shape, bbox, false);

  // --- compute sdf

  // https://iquilezles.org/articles/distfunctions2d/
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float dx = p_noise_x ? (*p_noise_x)(i, j) : 0.f;
      float dy = p_noise_y ? (*p_noise_y)(i, j) : 0.f;

      float x = xg[i] + dx;
      float y = yg[j] + dy;

      // --- distance

      float dmin = std::numeric_limits<float>::max();

      for (size_t i = 0; i < xp.size() - 1; i += 3)
      {
        int j = i + 1;
        int k = i + 2;

        float dist;

        glm::vec2 a = {xp[j] - xp[i], yp[j] - yp[i]};
        glm::vec2 b = {xp[i] - 2.f * xp[j] + xp[k],
                       yp[i] - 2.f * yp[j] + yp[k]};
        glm::vec2 c = {2.f * a.x, 2.f * a.y};
        glm::vec2 d = {xp[i] - x, yp[i] - y};

        float kk = 1.f / dot(b, b);
        float kx = kk * dot(a, b);
        float ky = kk * (2.f * dot(a, a) + dot(d, b)) / 3.f;
        float kz = kk * dot(d, a);
        float p = ky - kx * kx;
        float p3 = p * p * p;
        float q = kx * (2.f * kx * kx - 3.f * ky) + kz;
        float h = q * q + 4.f * p3;

        if (h >= 0.f)
        {
          h = std::sqrt(h);
          glm::vec2 xx = {0.5f * (h - q), 0.5f * (-h - q)};
          glm::vec2 uv = {
              std::copysign(1.f, xx.x) * std::pow(std::abs(xx.x), 1.f / 3.f),
              std::copysign(1.f, xx.y) * std::pow(std::abs(xx.y), 1.f / 3.f)};
          float t = std::clamp(uv.x + uv.y - kx, 0.f, 1.f);

          glm::vec2 dd = {d.x + (c.x + b.x * t) * t, d.y + (c.y + b.y * t) * t};

          dist = dot(dd, dd);
        }
        else
        {
          float     zz = std::sqrt(-p);
          float     v = std::acos(q / (p * zz * 2.f)) / 3.f;
          float     m = std::cos(v);
          float     n = std::sin(v) * 1.732050808f;
          glm::vec3 tt = {std::clamp((m + m) * zz - kx, 0.f, 1.f),
                          std::clamp((-n - m) * zz - kx, 0.f, 1.f),
                          std::clamp((n - m) * zz - kx, 0.f, 1.f)};

          glm::vec2 dd1 = {d.x + (c.x + b.x * tt.x) * tt.x,
                           d.y + (c.y + b.y * tt.x) * tt.x};
          glm::vec2 dd2 = {d.x + (c.x + b.x * tt.y) * tt.y,
                           d.y + (c.y + b.y * tt.y) * tt.y};

          float dist1 = dot(dd1, dd1);
          float dist2 = dot(dd2, dd2);

          dist = std::min(dist1, dist2);
        }

        dmin = std::min(dmin, dist);
      }

      sdf2(i, j) = std::sqrt(dmin);
    }

  return sdf2;
}

} // namespace hmap
