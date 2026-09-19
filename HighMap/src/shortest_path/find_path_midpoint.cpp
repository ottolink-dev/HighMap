/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

Point helper_pick_lowest(const Point &pa,
                         const Point &pb,
                         const Array &z,
                         float        max_offset,
                         int          steps)
{
  const glm::ivec2 &shape = z.shape;

  // midpoint
  Point pm(0.5f * (pa.x + pb.x), 0.5f * (pa.y + pb.y));

  // direction
  float dx = pb.x - pa.x;
  float dy = pb.y - pa.y;
  float len = std::sqrt(dx * dx + dy * dy);

  if (len < 1e-5f) return pm;

  // normalized perpendicular
  float px = -dy / len;
  float py = dx / len;

  // helper to clamp + sample
  auto sample = [&](const Point &p)
  {
    int i = int(std::round(p.x));
    int j = int(std::round(p.y));

    i = std::clamp(i, 0, shape.x - 1);
    j = std::clamp(j, 0, shape.y - 1);

    return z(i, j);
  };

  // initialize with midpoint
  Point best_pt = pm;
  float best_val = sample(pm);

  // ensure valid steps
  steps = std::max(steps, 1);

  // sample along perpendicular
  for (int s = 0; s <= steps; ++s)
  {
    float t = -max_offset + (2.f * max_offset) * (float(s) / float(steps));

    Point p(pm.x + px * t, pm.y + py * t);

    float val = sample(p);

    if (val < best_val)
    {
      best_val = val;
      best_pt = p;
    }
  }

  return best_pt;
}

std::vector<glm::ivec2> find_path_midpoint(const Array &z,
                                           glm::ivec2   ij_start,
                                           glm::ivec2   ij_end,
                                           float        offset_ratio,
                                           int          max_it,
                                           int          steps)
{
  if (!validate_non_empty(z)) return {};

  const glm::ivec2 &shape = z.shape;

  // --- adaptive number of iterations if not provided

  float dist = glm::length(glm::vec2(ij_end - ij_start));

  if (max_it == 0) max_it = int(std::ceil(std::log2(std::max(dist, 1.f))));

  int iterations = std::clamp(max_it, 2, 14);

  // --- midpoint

  std::vector<glm::ivec2> idx;
  idx.push_back(ij_start);
  idx.push_back(ij_end);

  for (int it = 0; it < iterations; ++it)
  {
    std::vector<glm::ivec2> new_idx;
    new_idx.reserve(idx.size() * 2);

    bool any_subdivided = false;
    for (size_t k = 0; k < idx.size() - 1; ++k)
    {
      const glm::ivec2 &a = idx[k];
      const glm::ivec2 &b = idx[k + 1];

      new_idx.push_back(a);

      // skip if already adjacent pixels (8-neighborhood)
      if (std::max(std::abs(a.x - b.x), std::abs(a.y - b.y)) <= 1) continue;

      Point pa(float(a.x), float(a.y));
      Point pb(float(b.x), float(b.y));

      // adaptive offset based on segment length
      float seg_len = glm::length(glm::vec2(b - a));
      float max_offset = offset_ratio * seg_len;

      // use helper
      Point best = helper_pick_lowest(pa, pb, z, max_offset, steps);

      // convert to grid
      glm::ivec2 ij(int(std::round(best.x)), int(std::round(best.y)));

      ij.x = std::clamp(ij.x, 0, shape.x - 1);
      ij.y = std::clamp(ij.y, 0, shape.y - 1);

      // avoid duplicates
      if (ij != new_idx.back() && ij != b)
      {
        new_idx.push_back(ij);
        any_subdivided = true;
      }
    }

    new_idx.push_back(idx.back());

    // early exit if stabilized or no segment needed subdivision
    if (!any_subdivided || new_idx.size() == idx.size()) break;

    idx = std::move(new_idx);
  }

  return idx;
}

void find_path_midpoint(const Array      &z,
                        glm::ivec2        ij_start,
                        glm::ivec2        ij_end,
                        std::vector<int> &i_path,
                        std::vector<int> &j_path,
                        float             offset_ratio,
                        int               max_it,
                        int               steps)
{
  std::vector<glm::ivec2> indices = find_path_midpoint(z,
                                                       ij_start,
                                                       ij_end,
                                                       offset_ratio,
                                                       max_it,
                                                       steps);

  i_path.clear();
  j_path.clear();
  i_path.reserve(indices.size());
  j_path.reserve(indices.size());

  for (const auto &p : indices)
  {
    i_path.push_back(p.x);
    j_path.push_back(p.y);
  }
}

Path find_path_midpoint(const Array &z,
                        glm::ivec2   ij_start,
                        glm::ivec2   ij_end,
                        glm::vec4    bbox,
                        float        offset_ratio,
                        int          max_it,
                        int          steps)
{
  if (!validate_non_empty(z)) return Path();

  std::vector<glm::ivec2> indices = find_path_midpoint(z,
                                                       ij_start,
                                                       ij_end,
                                                       offset_ratio,
                                                       max_it,
                                                       steps);

  std::vector<Point> points;
  points.reserve(indices.size());

  const float lx = bbox.y - bbox.x;
  const float ly = bbox.w - bbox.z;
  const float denom_x = (z.shape.x > 1) ? float(z.shape.x - 1) : 1.f;
  const float denom_y = (z.shape.y > 1) ? float(z.shape.y - 1) : 1.f;

  for (const auto &p : indices)
  {
    float px = (float(p.x) / denom_x) * lx + bbox.x;
    float py = (float(p.y) / denom_y) * ly + bbox.z;
    float pv = z(p);

    points.emplace_back(px, py, pv);
  }

  return Path(points);
}

} // namespace hmap
