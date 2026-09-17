/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

namespace
{

inline Point get_midpoint(const Point &p1, const Point &p2)
{
  return Point((p1.x + p2.x) * 0.5f,
               (p1.y + p2.y) * 0.5f,
               (p1.v + p2.v) * 0.5f);
}

// 64-bit deterministic hash with Murmur3 / SplitMix64 finalizer
uint64_t hash_edge_coords(int64_t       x1,
                          int64_t       y1,
                          int64_t       x2,
                          int64_t       y2,
                          std::uint32_t seed)
{
  uint64_t hash = 14695981039346656037ULL;
  auto     add_v = [&](int64_t v)
  {
    hash ^= static_cast<uint64_t>(v);
    hash *= 1099511628211ULL;
  };
  add_v(x1);
  add_v(y1);
  add_v(x2);
  add_v(y2);
  add_v(static_cast<int64_t>(seed));

  hash ^= hash >> 33;
  hash *= 0xff51afd7ed558ccdULL;
  hash ^= hash >> 33;
  hash *= 0xc4ceb9fe1a85ec53ULL;
  hash ^= hash >> 33;

  return hash;
}

// Selects crossing segment 0 or 1 for the canonical edge [p_min, p_max]
int select_canonical_segment(const Point  &p_min,
                             const Point  &p_max,
                             std::uint32_t seed,
                             const Array  *p_weights,
                             const Array  *p_mask,
                             glm::vec4     bbox)
{
  Point m = get_midpoint(p_min, p_max);
  Point m0 = get_midpoint(p_min, m);
  Point m1 = get_midpoint(m, p_max);

  // Mask constraint: if mask is provided, avoid forbidden regions
  if (p_mask)
  {
    bool valid0 = p_mask->get_value_nearest(m0.x, m0.y, bbox) > 0.5f;
    bool valid1 = p_mask->get_value_nearest(m1.x, m1.y, bbox) > 0.5f;

    if (valid0 && !valid1) return 0;
    if (valid1 && !valid0) return 1;
  }

  // Quantize coordinates to integers for robust canonical hashing
  int64_t qx1 = static_cast<int64_t>(std::round(p_min.x * 1000000.0));
  int64_t qy1 = static_cast<int64_t>(std::round(p_min.y * 1000000.0));
  int64_t qx2 = static_cast<int64_t>(std::round(p_max.x * 1000000.0));
  int64_t qy2 = static_cast<int64_t>(std::round(p_max.y * 1000000.0));

  uint64_t h = hash_edge_coords(qx1, qy1, qx2, qy2, seed);
  float    u = static_cast<float>((h & 0xFFFFFFFFULL)) / 4294967296.0f;

  float p0 = 0.5f;
  if (p_weights)
  {
    float w0 = std::max(0.0001f,
                        p_weights->get_value_nearest(m0.x, m0.y, bbox));
    float w1 = std::max(0.0001f,
                        p_weights->get_value_nearest(m1.x, m1.y, bbox));
    p0 = w0 / (w0 + w1);
  }

  return (u < p0) ? 0 : 1;
}

// Selects crossing segment for the directed edge (p_start -> p_end):
// returns 0 for [p_start, midpoint], 1 for [midpoint, p_end]
int select_edge_segment(const Point  &p_start,
                        const Point  &p_end,
                        std::uint32_t seed,
                        const Array  *p_weights,
                        const Array  *p_mask,
                        glm::vec4     bbox)
{
  bool start_is_min = (p_start.x < p_end.x) ||
                      (p_start.x == p_end.x && p_start.y <= p_end.y);
  const Point &p_min = start_is_min ? p_start : p_end;
  const Point &p_max = start_is_min ? p_end : p_start;

  int canonical_choice = select_canonical_segment(p_min,
                                                  p_max,
                                                  seed,
                                                  p_weights,
                                                  p_mask,
                                                  bbox);
  return start_is_min ? canonical_choice : (1 - canonical_choice);
}

// Subdivides an edge (p1 -> p2) within triangle (apex, p1, p2)
// using the 4 productions of the squiggle curve (Prusinkiewicz et al.)
std::vector<Point> subdivide_edge_squiggle(const Point  &p1,
                                           const Point  &p2,
                                           const Point  &apex,
                                           std::uint32_t seed,
                                           const Array  *p_weights,
                                           const Array  *p_mask,
                                           glm::vec4     bbox)
{
  Point m_in = get_midpoint(p1, apex);
  Point m_out = get_midpoint(p2, apex);
  Point m_neut = get_midpoint(p1, p2);

  int c_in = select_edge_segment(p1, apex, seed, p_weights, p_mask, bbox);
  int c_out = select_edge_segment(p2,
                                  apex,
                                  seed + 1013904223U,
                                  p_weights,
                                  p_mask,
                                  bbox);

  if (c_in == 0 && c_out == 0)
  {
    // Production 1: T_top only
    return {m_in, apex, m_out};
  }
  else if (c_in == 0 && c_out == 1)
  {
    // Production 2: T_top -> T_mid -> T_right
    return {m_in, apex, m_out, m_neut};
  }
  else if (c_in == 1 && c_out == 0)
  {
    // Production 3: T_left -> T_mid -> T_top
    return {m_neut, m_in, apex, m_out};
  }
  else
  {
    // Production 4: T_left -> T_mid -> T_right (s-curve)
    return {m_in, m_neut, m_out};
  }
}

} // namespace

Path squiggle(const Path   &path,
              int           iterations,
              std::uint32_t seed,
              float         height_ratio,
              int           orientation,
              const Array  *p_weights,
              const Array  *p_mask,
              glm::vec4     bbox)
{
  if (p_weights && !validate_non_empty(*p_weights)) return path;
  if (p_mask && !validate_non_empty(*p_mask)) return path;
  if (!validate_min_size(path, 2, "Path points")) return path;

  iterations = std::max(1, iterations);
  bool is_closed = path.is_closed();

  std::vector<Point> current_points(path.begin(), path.end());

  for (int it = 0; it < iterations; ++it)
  {
    std::vector<Point> new_points;
    size_t             n = current_points.size();
    size_t             num_edges = is_closed ? n : (n > 0 ? n - 1 : 0);
    float              current_sign = 1.0f;

    for (size_t k = 0; k < num_edges; ++k)
    {
      const Point &p1 = current_points[k];
      const Point &p2 = current_points[(k + 1) % n];

      float dx = p2.x - p1.x;
      float dy = p2.y - p1.y;
      float dist = std::hypot(dx, dy);
      if (dist < 1e-7f) continue;

      float nx = -dy / dist;
      float ny = dx / dist;

      float sign = 1.0f;
      if (orientation == 0)
      {
        sign = current_sign;
        current_sign *= -1.0f;
      }
      else if (orientation > 0)
      {
        sign = 1.0f;
      }
      else
      {
        sign = -1.0f;
      }

      Point apex(0.5f * (p1.x + p2.x) + nx * dist * height_ratio * sign,
                 0.5f * (p1.y + p2.y) + ny * dist * height_ratio * sign,
                 0.5f * (p1.v + p2.v));

      std::uint32_t edge_seed = seed + static_cast<std::uint32_t>(it * 10007 +
                                                                  k * 7919);
      std::vector<Point> sub_pts = subdivide_edge_squiggle(p1,
                                                           p2,
                                                           apex,
                                                           edge_seed,
                                                           p_weights,
                                                           p_mask,
                                                           bbox);

      if (new_points.empty())
      {
        new_points.push_back(p1);
      }
      for (const auto &sp : sub_pts)
      {
        if (distance(sp, new_points.back()) > 1e-6f)
        {
          new_points.push_back(sp);
        }
      }
      if (distance(p2, new_points.back()) > 1e-6f)
      {
        new_points.push_back(p2);
      }
    }

    if (is_closed && !new_points.empty())
    {
      if (distance(new_points.front(), new_points.back()) > 1e-6f)
      {
        new_points.push_back(new_points.front());
      }
    }

    // Clean duplicate consecutive points
    if (!new_points.empty())
    {
      std::vector<Point> clean_pts;
      clean_pts.reserve(new_points.size());
      clean_pts.push_back(new_points[0]);
      for (size_t i = 1; i < new_points.size(); ++i)
      {
        if (distance(new_points[i], clean_pts.back()) > 1e-6f)
        {
          clean_pts.push_back(new_points[i]);
        }
      }
      current_points = std::move(clean_pts);
    }
  }

  Path out_path(current_points);
  out_path.set_closed(is_closed);
  return out_path;
}

std::vector<Path> squiggle_branches(const Path   &path,
                                    int           iterations,
                                    std::uint32_t seed,
                                    float         branch_probability,
                                    int           max_branches,
                                    float         height_ratio,
                                    int           orientation,
                                    const Array  *p_weights,
                                    const Array  *p_mask,
                                    glm::vec4     bbox)
{
  if (p_weights && !validate_non_empty(*p_weights)) return {path};
  if (p_mask && !validate_non_empty(*p_mask)) return {path};
  if (!validate_min_size(path, 2, "Path points")) return {path};

  iterations = std::max(1, iterations);
  bool is_closed = path.is_closed();

  std::vector<Point> current_points(path.begin(), path.end());
  std::vector<Path>  branches;

  for (int it = 0; it < iterations; ++it)
  {
    std::vector<Point> new_points;
    size_t             n = current_points.size();
    size_t             num_edges = is_closed ? n : (n > 0 ? n - 1 : 0);
    float              current_sign = 1.0f;

    for (size_t k = 0; k < num_edges; ++k)
    {
      const Point &p1 = current_points[k];
      const Point &p2 = current_points[(k + 1) % n];

      float dx = p2.x - p1.x;
      float dy = p2.y - p1.y;
      float dist = std::hypot(dx, dy);
      if (dist < 1e-7f) continue;

      float nx = -dy / dist;
      float ny = dx / dist;

      float sign = 1.0f;
      if (orientation == 0)
      {
        sign = current_sign;
        current_sign *= -1.0f;
      }
      else if (orientation > 0)
      {
        sign = 1.0f;
      }
      else
      {
        sign = -1.0f;
      }

      Point apex(0.5f * (p1.x + p2.x) + nx * dist * height_ratio * sign,
                 0.5f * (p1.y + p2.y) + ny * dist * height_ratio * sign,
                 0.5f * (p1.v + p2.v));

      std::uint32_t edge_seed = seed + static_cast<std::uint32_t>(it * 10007 +
                                                                  k * 7919);

      // Check branch generation if allowed
      if (branch_probability > 0.f &&
          static_cast<int>(branches.size()) < max_branches && it >= 1)
      {
        uint64_t bh = hash_edge_coords(
            static_cast<int64_t>(std::round(p1.x * 1000.0)),
            static_cast<int64_t>(std::round(p1.y * 1000.0)),
            static_cast<int64_t>(it),
            static_cast<int64_t>(k),
            edge_seed + 777U);
        float bu = static_cast<float>(bh & 0xFFFFFFFFULL) / 4294967296.0f;
        if (bu < branch_probability)
        {
          Point branch_start = get_midpoint(p1, p2);
          Point branch_end = apex;
          Path  branch_base({branch_start, branch_end});
          int   branch_iters = std::max(1, iterations - it);
          Path  branch_path = squiggle(branch_base,
                                      branch_iters,
                                      edge_seed ^ 0xA5A5A5A5U,
                                      height_ratio * 0.7f,
                                      -orientation,
                                      p_weights,
                                      p_mask,
                                      bbox);
          if (branch_path.size() >= 2)
          {
            branches.push_back(branch_path);
          }
        }
      }

      std::vector<Point> sub_pts = subdivide_edge_squiggle(p1,
                                                           p2,
                                                           apex,
                                                           edge_seed,
                                                           p_weights,
                                                           p_mask,
                                                           bbox);

      if (new_points.empty())
      {
        new_points.push_back(p1);
      }
      for (const auto &sp : sub_pts)
      {
        if (distance(sp, new_points.back()) > 1e-6f)
        {
          new_points.push_back(sp);
        }
      }
      if (distance(p2, new_points.back()) > 1e-6f)
      {
        new_points.push_back(p2);
      }
    }

    if (is_closed && !new_points.empty())
    {
      if (distance(new_points.front(), new_points.back()) > 1e-6f)
      {
        new_points.push_back(new_points.front());
      }
    }

    if (!new_points.empty())
    {
      std::vector<Point> clean_pts;
      clean_pts.reserve(new_points.size());
      clean_pts.push_back(new_points[0]);
      for (size_t i = 1; i < new_points.size(); ++i)
      {
        if (distance(new_points[i], clean_pts.back()) > 1e-6f)
        {
          clean_pts.push_back(new_points[i]);
        }
      }
      current_points = std::move(clean_pts);
    }
  }

  Path main_path(current_points);
  main_path.set_closed(is_closed);

  std::vector<Path> all_paths;
  all_paths.push_back(main_path);
  for (const auto &b : branches)
  {
    all_paths.push_back(b);
  }

  return all_paths;
}

} // namespace hmap
