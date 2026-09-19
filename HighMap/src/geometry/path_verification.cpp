/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

#include "highmap/geometry/kd_tree.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/logger.hpp"

namespace hmap
{

// --- Functions

bool assert_start_end_points(const Path &path1,
                             const Path &path2,
                             float       tol,
                             bool        verbose)
{
  if (path1.size() < 1 || path2.size() < 1)
  {
    hmap::log::error("not enough Path points");
    return false;
  }

  const Point &p1s = path1.front();
  const Point &p1e = path1.back();
  const Point &p2s = path2.front();
  const Point &p2e = path2.back();

  float ds = distance(p1s, p2s);
  float de = distance(p1e, p2e);
  bool  assert = ds < tol && de < tol;

  if (verbose)
    hmap::log::info("ds: {}, de: {}, assert: {}", ds, de, assert ? "T" : "F");

  return assert;
}

float chamfer_distance(const Path &a, const Path &b)
{
  if (a.empty() || b.empty()) return 0.f;

  KDTree tree_b(b);
  float  sum_a = 0.f;
  for (const auto &pa : a)
  {
    auto [idx, dist_sq] = tree_b.nearest_with_distance_squared(pa);
    sum_a += std::sqrt(dist_sq);
  }

  KDTree tree_a(a);
  float  sum_b = 0.f;
  for (const auto &pb : b)
  {
    auto [idx, dist_sq] = tree_a.nearest_with_distance_squared(pb);
    sum_b += std::sqrt(dist_sq);
  }

  return (sum_a / float(a.size())) + (sum_b / float(b.size()));
}

bool has_duplicates(const Path &path, float tol)
{
  if (path.size() < 2) return false;

  KDTree tree(path);
  for (const auto &p : path)
  {
    auto neighbors = tree.radius_search(p, tol);
    // if radius_search finds more than 1 point within tolerance, there is a
    // duplicate
    if (neighbors.size() > 1) return true;
  }

  return false;
}

} // namespace hmap
