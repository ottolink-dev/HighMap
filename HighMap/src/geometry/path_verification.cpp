/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cstddef>
#include <vector>

#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/logger.hpp"

#include <cfloat>

namespace hmap
{

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
  auto avg = [](const Path &p, const Path &q)
  {
    float sum = 0.f;
    for (const auto &pa : p)
    {
      float min_d = FLT_MAX;
      for (const auto &pb : q)
        min_d = std::min(min_d, distance(pa, pb));
      sum += min_d;
    }
    return sum / float(p.size());
  };

  return avg(a, b) + avg(b, a);
}

bool has_duplicates(const Path &path, float tol)
{
  for (size_t i = 0; i < path.size(); ++i)
    for (size_t j = i + 1; j < path.size(); ++j)
      if (distance(path[i], path[j]) < tol) return true;

  return false;
}

} // namespace hmap
