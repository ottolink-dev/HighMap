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

struct Bezier2D
{
  glm::vec2 p0{0.f};
  glm::vec2 p1{0.f};
  glm::vec2 p2{0.f};

  constexpr AABB2D bounds(float padding = 0.f) const noexcept
  {
    return AABB2D(std::min({p0.x, p1.x, p2.x}) - padding,
                  std::min({p0.y, p1.y, p2.y}) - padding,
                  std::max({p0.x, p1.x, p2.x}) + padding,
                  std::max({p0.y, p1.y, p2.y}) + padding);
  }
};

inline float bezier_distance_sq(const glm::vec2 &pos, const Bezier2D &bez)
{
  glm::vec2 a = bez.p1 - bez.p0;
  glm::vec2 b = bez.p0 - 2.f * bez.p1 + bez.p2;
  glm::vec2 c = 2.f * a;
  glm::vec2 d = bez.p0 - pos;

  float b_len_sq = dot(b, b);
  if (b_len_sq < 1e-10f)
  {
    // Colinear degenerate quadratic bezier: distance to segment [p0, p2]
    glm::vec2 e = bez.p2 - bez.p0;
    glm::vec2 w = pos - bez.p0;
    float     lsq = dot(e, e);
    float     t = (lsq > 0.f) ? std::clamp(dot(w, e) / lsq, 0.f, 1.f) : 0.f;
    glm::vec2 diff = w - e * t;
    return dot(diff, diff);
  }

  float kk = 1.f / b_len_sq;
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
    float     t = std::clamp(uv.x + uv.y - kx, 0.f, 1.f);
    glm::vec2 dd = {d.x + (c.x + b.x * t) * t, d.y + (c.y + b.y * t) * t};
    return dot(dd, dd);
  }
  else
  {
    float zz = std::sqrt(std::max(0.f, -p));
    float denom = p * zz * 2.f;
    float arg = (std::abs(denom) > 1e-12f) ? std::clamp(q / denom, -1.f, 1.f)
                                           : 0.f;
    float v = std::acos(arg) / 3.f;
    float m = std::cos(v);
    float n = std::sin(v) * 1.732050808f;
    glm::vec3 tt = {std::clamp((m + m) * zz - kx, 0.f, 1.f),
                    std::clamp((-n - m) * zz - kx, 0.f, 1.f),
                    std::clamp((n - m) * zz - kx, 0.f, 1.f)};

    glm::vec2 dd1 = {d.x + (c.x + b.x * tt.x) * tt.x,
                     d.y + (c.y + b.y * tt.x) * tt.x};
    glm::vec2 dd2 = {d.x + (c.x + b.x * tt.y) * tt.y,
                     d.y + (c.y + b.y * tt.y) * tt.y};

    float dist1 = dot(dd1, dd1);
    float dist2 = dot(dd2, dd2);
    return std::min(dist1, dist2);
  }
}

// Ray-casting parity test for quadratic bezier curve:
// counts intersections with a horizontal ray cast from p towards +x.
int count_bezier_ray_crossings(const Bezier2D &b, const glm::vec2 &p)
{
  float a = b.p0.y - 2.f * b.p1.y + b.p2.y;
  float b_coeff = 2.f * (b.p1.y - b.p0.y);
  float c = b.p0.y - p.y;

  int crossings = 0;

  auto check_t = [&](float t)
  {
    // Use [0, 1) interval to avoid double-counting at shared curve endpoints
    if (t >= 0.f && t < 1.f)
    {
      float bx = (1.f - t) * (1.f - t) * b.p0.x + 2.f * t * (1.f - t) * b.p1.x +
                 t * t * b.p2.x;
      if (bx > p.x) crossings++;
    }
  };

  if (std::abs(a) < 1e-9f)
  {
    if (std::abs(b_coeff) > 1e-9f) check_t(-c / b_coeff);
  }
  else
  {
    float discr = b_coeff * b_coeff - 4.f * a * c;
    if (discr >= 0.f)
    {
      float s = std::sqrt(discr);
      check_t((-b_coeff - s) / (2.f * a));
      check_t((-b_coeff + s) / (2.f * a));
    }
  }

  return crossings;
}

std::vector<Bezier2D> extract_bezier_curves(const Path &path)
{
  std::vector<float> x_input = path.get_x();
  std::vector<float> y_input = path.get_y();

  if (path.is_closed() && !x_input.empty())
  {
    x_input.push_back(x_input.front());
    y_input.push_back(y_input.front());
  }

  if (x_input.size() % 2 == 0)
  {
    size_t k = x_input.size();
    float  new_x = lerp(x_input[k - 1], x_input[k - 2], 0.1f);
    float  new_y = lerp(y_input[k - 1], y_input[k - 2], 0.1f);

    x_input.insert(x_input.end() - 1, new_x);
    y_input.insert(y_input.end() - 1, new_y);
  }

  std::vector<Bezier2D> curves;
  curves.reserve((x_input.size() - 1) / 2);

  for (size_t k = 0; k < x_input.size() - 2; k += 2)
  {
    Bezier2D bez;
    bez.p0 = {x_input[k], y_input[k]};
    bez.p1 = {x_input[k + 1], y_input[k + 1]};
    bez.p2 = {x_input[k + 2], y_input[k + 2]};
    curves.push_back(bez);
  }

  return curves;
}

class BVH2DBezier
{
public:
  struct Node
  {
    AABB2D   bbox;
    int32_t  left{-1};
    int32_t  right{-1};
    uint32_t curve_start{0};
    uint32_t curve_count{0};

    constexpr bool is_leaf() const noexcept
    {
      return left < 0;
    }
  };

  BVH2DBezier() = default;

  explicit BVH2DBezier(const std::vector<Bezier2D> &input_curves)
  {
    this->build(input_curves);
  }

  void build(const std::vector<Bezier2D> &input_curves)
  {
    this->curves = input_curves;
    this->nodes.clear();
    if (this->curves.empty()) return;

    this->nodes.reserve(this->curves.size() * 2);
    this->build_recursive(0, static_cast<uint32_t>(this->curves.size()));
  }

  const std::vector<Node> &get_nodes() const noexcept
  {
    return nodes;
  }

  const std::vector<Bezier2D> &get_curves() const noexcept
  {
    return curves;
  }

  bool empty() const noexcept
  {
    return nodes.empty();
  }

private:
  std::vector<Node>     nodes;
  std::vector<Bezier2D> curves;

  int32_t build_recursive(uint32_t start, uint32_t count)
  {
    Node node;
    node.curve_start = start;
    node.curve_count = count;

    AABB2D node_bbox;
    AABB2D centroid_bbox;

    for (uint32_t k = 0; k < count; ++k)
    {
      const auto &c = this->curves[start + k];
      node_bbox.expand(c.bounds());
      glm::vec2 centroid = (c.p0 + c.p1 + c.p2) / 3.f;
      centroid_bbox.expand(centroid);
    }

    node.bbox = node_bbox;

    if (count <= 4)
    {
      int32_t node_idx = static_cast<int32_t>(this->nodes.size());
      this->nodes.push_back(node);
      return node_idx;
    }

    glm::vec2 ext = centroid_bbox.extent();
    int       axis = (ext.x >= ext.y) ? 0 : 1;

    auto begin_it = this->curves.begin() + start;
    auto mid_it = begin_it + (count / 2);
    auto end_it = begin_it + count;

    if (axis == 0)
    {
      std::nth_element(begin_it,
                       mid_it,
                       end_it,
                       [](const Bezier2D &c1, const Bezier2D &c2)
                       {
                         return (c1.p0.x + c1.p1.x + c1.p2.x) <
                                (c2.p0.x + c2.p1.x + c2.p2.x);
                       });
    }
    else
    {
      std::nth_element(begin_it,
                       mid_it,
                       end_it,
                       [](const Bezier2D &c1, const Bezier2D &c2)
                       {
                         return (c1.p0.y + c1.p1.y + c1.p2.y) <
                                (c2.p0.y + c2.p1.y + c2.p2.y);
                       });
    }

    uint32_t left_count = count / 2;
    uint32_t right_count = count - left_count;

    int32_t node_idx = static_cast<int32_t>(this->nodes.size());
    this->nodes.push_back(node);

    int32_t left_child = this->build_recursive(start, left_count);
    int32_t right_child = this->build_recursive(start + left_count,
                                                right_count);

    this->nodes[node_idx].left = left_child;
    this->nodes[node_idx].right = right_child;

    return node_idx;
  }
};

float query_nearest_bezier_dist_sq_bvh(const BVH2DBezier &bvh,
                                       const glm::vec2   &p)
{
  const auto &nodes = bvh.get_nodes();
  const auto &curves = bvh.get_curves();
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
      for (uint32_t k = 0; k < node.curve_count; ++k)
      {
        float d2 = bezier_distance_sq(p, curves[node.curve_start + k]);
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

Array sdf_2d_polyline_bezier(const Path  &path,
                             glm::ivec2   shape,
                             glm::vec4    bbox,
                             const Array *p_noise_x,
                             const Array *p_noise_y,
                             bool         use_bvh)
{
  if (use_bvh)
    return sdf_2d_polyline_bezier_bvh(path, shape, bbox, p_noise_x, p_noise_y);

  return sdf_2d_polyline_bezier_brute_force(path,
                                            shape,
                                            bbox,
                                            p_noise_x,
                                            p_noise_y);
}

Array sdf_2d_polyline_bezier_brute_force(const Path  &path,
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

  std::vector<Bezier2D> curves = extract_bezier_curves(path);
  if (curves.empty()) return Array(shape);

  Array sdf2(shape);

  // array cell coordinates
  std::vector<float> xg(shape.x);
  std::vector<float> yg(shape.y);
  grid_xy_vector(xg, yg, shape, bbox, false);

  const bool is_closed = path.is_closed();

  // compute sdf
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float dx = p_noise_x ? (*p_noise_x)(i, j) : 0.f;
      float dy = p_noise_y ? (*p_noise_y)(i, j) : 0.f;

      float x = xg[i] + dx;
      float y = yg[j] + dy;

      glm::vec2 p(x, y);
      float     dist_sq = std::numeric_limits<float>::max();
      int       crossings = 0;

      for (const auto &bez : curves)
      {
        float d2 = bezier_distance_sq(p, bez);
        dist_sq = std::min(dist_sq, d2);

        if (is_closed) crossings += count_bezier_ray_crossings(bez, p);
      }

      float sign = (is_closed && (crossings % 2 == 1)) ? -1.f : 1.f;
      sdf2(i, j) = sign * std::sqrt(dist_sq);
    }

  return sdf2;
}

Array sdf_2d_polyline_bezier_bvh(const Path  &path,
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

  std::vector<Bezier2D> curves = extract_bezier_curves(path);
  if (curves.empty()) return Array(shape);

  BVH2DBezier bvh(curves);

  Array sdf2(shape);

  // array cell coordinates
  std::vector<float> xg(shape.x);
  std::vector<float> yg(shape.y);
  grid_xy_vector(xg, yg, shape, bbox, false);

  const bool is_closed = path.is_closed();

  // compute sdf
#pragma omp parallel for collapse(2) schedule(static)
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float dx = p_noise_x ? (*p_noise_x)(i, j) : 0.f;
      float dy = p_noise_y ? (*p_noise_y)(i, j) : 0.f;

      float x = xg[i] + dx;
      float y = yg[j] + dy;

      glm::vec2 p(x, y);
      float     dist_sq = query_nearest_bezier_dist_sq_bvh(bvh, p);

      float sign = 1.f;
      if (is_closed)
      {
        int crossings = 0;
        for (const auto &bez : curves)
          crossings += count_bezier_ray_crossings(bez, p);

        if (crossings % 2 == 1) sign = -1.f;
      }

      sdf2(i, j) = sign * std::sqrt(dist_sq);
    }

  return sdf2;
}

} // namespace hmap
