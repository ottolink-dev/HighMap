/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <random>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate1d.hpp"
#include "highmap/interpolate/interpolate_curve.hpp"
#include "highmap/operator.hpp"

namespace hmap
{

// --- HELPERS

Path helper_build_path(const std::vector<Point> &points,
                       InterpolationMethodCurve  method,
                       Path::EdgeDivisionMode    edm,
                       int                       edge_divisions,
                       int                       point_count)
{
  InterpolatorCurve fitp(points, method);

  int npts = edge_divisions;
  if (edm == Path::EdgeDivisionMode::EDM_PER_EDGE) npts *= point_count;

  std::vector<float> t = hmap::linspace(0.f, 1.f, npts);

  return Path(fitp(t));
}

Path bezier(const Path            &path,
            float                  curvature_ratio,
            int                    edge_divisions,
            Path::EdgeDivisionMode edm)
{
  if (!validate_min_size(path, 2, "Path points")) return path;

  // --- generate a new set of points by adding control points
  // --- inbetween path points

  std::vector<Point> new_points = {};
  new_points.reserve(path.size() + 2 * (path.size() - 1));

  size_t npoints = path.size();
  size_t end = path.is_closed() ? npoints : npoints - 1;

  for (size_t k = 0; k < end; k++)
  {
    size_t knext = (k + 1) % npoints;
    size_t knext_after = (k + 2) % npoints;

    Point pc1 = lerp(path[k], path[knext], curvature_ratio);
    Point pc2 = lerp(path[knext], path[knext_after], -curvature_ratio);

    new_points.push_back(path[k]);
    new_points.push_back(pc1);
    new_points.push_back(pc2);
  }

  // preserve endpoint for open paths
  if (!path.is_closed())
  {
    new_points.push_back(path.back());
  }
  else
  {
    // keep loop continuity
    new_points.push_back(path.front());
  }

  // --- interpolate

  Path new_path = helper_build_path(new_points,
                                    InterpolationMethodCurve::BEZIER,
                                    edm,
                                    edge_divisions,
                                    path.size());

  if (path.is_closed()) new_path.pop_back();
  new_path.set_closed(path.is_closed());

  return new_path;
}

Path bezier_round(const Path            &path,
                  float                  curvature_ratio,
                  int                    edge_divisions,
                  Path::EdgeDivisionMode edm)
{
  if (!validate_min_size(path, 2, "Path points")) return path;

  // --- generate a new set of points by adding control points
  // --- inbetween path points

  std::vector<Point> new_points = {};
  new_points.reserve(path.size() + 2 * (path.size() - 1));

  size_t npoints = path.size();
  size_t end = path.is_closed() ? npoints : npoints - 1;

  for (size_t k = 0; k < end; k++)
  {
    size_t kprev = (k - 1) % npoints;
    size_t knext = (k + 1) % npoints;
    size_t knext_after = (k + 2) % npoints;

    Point delta_p1 = path[knext] - path[kprev];
    Point delta_p2 = path[k] - path[knext_after];

    Point pc1 = path[k] + curvature_ratio * delta_p1;
    Point pc2 = path[knext] + curvature_ratio * delta_p2;

    new_points.push_back(path[k]);
    new_points.push_back(pc1);
    new_points.push_back(pc2);
  }

  // preserve endpoint for open paths
  if (!path.is_closed())
  {
    new_points.push_back(path.back());
  }
  else
  {
    // keep loop continuity
    new_points.push_back(path.front());
  }

  // --- interpolate

  Path new_path = helper_build_path(new_points,
                                    InterpolationMethodCurve::BEZIER,
                                    edm,
                                    edge_divisions,
                                    path.size());

  if (path.is_closed()) new_path.pop_back();
  new_path.set_closed(path.is_closed());

  return new_path;
}

Path bspline(const Path &path, int edge_divisions, Path::EdgeDivisionMode edm)
{
  if (!validate_min_size(path, 2, "Path points")) return path;

  std::vector<Point> new_points(path.begin(), path.end());

  if (path.is_closed()) new_points.push_back(path.front());

  Path new_path = helper_build_path(new_points,
                                    InterpolationMethodCurve::BSPLINE,
                                    edm,
                                    edge_divisions,
                                    path.size());

  if (path.is_closed()) new_path.pop_back();
  new_path.set_closed(path.is_closed());

  return new_path;
}

Path catmullrom(const Path            &path,
                int                    edge_divisions,
                Path::EdgeDivisionMode edm)
{
  if (!validate_min_size(path, 2, "Path points")) return path;

  std::vector<Point> new_points(path.begin(), path.end());

  if (path.is_closed()) new_points.push_back(path.front());

  Path new_path = helper_build_path(new_points,
                                    InterpolationMethodCurve::CATMULLROM,
                                    edm,
                                    edge_divisions,
                                    path.size());

  if (path.is_closed()) new_path.pop_back();
  new_path.set_closed(path.is_closed());

  return new_path;
}

Path decasteljau(const Path            &path,
                 int                    edge_divisions,
                 Path::EdgeDivisionMode edm)
{
  if (!validate_min_size(path, 2, "Path points")) return path;

  std::vector<Point> new_points(path.begin(), path.end());

  if (path.is_closed()) new_points.push_back(path.front());

  Path new_path = helper_build_path(new_points,
                                    InterpolationMethodCurve::DECASTELJAU,
                                    edm,
                                    edge_divisions,
                                    path.size());

  if (path.is_closed()) new_path.pop_back();
  new_path.set_closed(path.is_closed());

  return new_path;
}

Path decimate_vw(const Path &path, int n_points_target)
{
  size_t n = path.size();
  if (n < 3 || n <= (size_t)n_points_target) return path;

  Path new_path = path;

  while (new_path.size() > (size_t)n_points_target)
  {
    size_t remove_idx = 1;
    float  min_area = std::numeric_limits<float>::max();

    // find smallest effective triangle
    for (size_t i = 1; i + 1 < new_path.size(); ++i)
    {
      float area = triangle_area(new_path[i - 1], new_path[i], new_path[i + 1]);
      if (area < min_area)
      {
        min_area = area;
        remove_idx = i;
      }
    }

    new_path.erase(new_path.begin() + remove_idx);
  }

  return new_path;
}

Path fractalize(const Path   &path,
                int           iterations,
                std::uint32_t seed,
                float         sigma,
                int           orientation,
                float         persistence,
                Array        *p_ctrl_array,
                glm::vec4     bbox,
                bool          bounded)
{
  if (!validate_min_size(path, 2, "Path points")) return path;
  if (p_ctrl_array && !validate_non_empty(*p_ctrl_array)) return path;

  Path new_path = path;

  std::mt19937                    gen(seed);
  std::normal_distribution<float> dis(0.f, 1.f);

  struct EdgeBBox
  {
    float min_x, max_x, min_y, max_y;
  };
  std::vector<EdgeBBox> initial_bboxes;

  if (bounded)
  {
    size_t n_init = path.size();
    size_t end_init = path.is_closed() ? n_init : (n_init > 0 ? n_init - 1 : 0);
    initial_bboxes.reserve(end_init);

    for (size_t i = 0; i < end_init; ++i)
    {
      size_t      inext = (i + 1) % n_init;
      const auto &p1 = path[i];
      const auto &p2 = path[inext];
      initial_bboxes.push_back({
          std::min(p1.x, p2.x),
          std::max(p1.x, p2.x),
          std::min(p1.y, p2.y),
          std::max(p1.y, p2.y),
      });
    }
  }

  for (int it = 0; it < iterations; it++)
  {
    std::vector<Point> new_points = {};

    // determine the ending index based on whether the list is
    // closed (circular)
    size_t npoints = new_path.size();
    size_t end = new_path.is_closed() ? npoints : npoints - 1;

    for (size_t k = 0; k < end; k++)
    {
      // determine the index of the next point, wrapping around if circular
      size_t knext = (k + 1) % npoints;

      // generate random displacement amplitude (as a ratio to the
      // point distance)
      float amp = sigma * dis(gen);

      // if provided, modulate amplitude based on underlying field
      if (p_ctrl_array)
        amp *= p_ctrl_array->get_value_nearest(new_path[k].x,
                                               new_path[k].y,
                                               bbox);

      // insert midpoint between current edge start and end
      Point pnew = midpoint(new_path[k], new_path[knext], orientation, amp);

      if (bounded && !initial_bboxes.empty())
      {
        size_t edge_idx = k >> it;
        if (edge_idx < initial_bboxes.size())
        {
          const auto &ibox = initial_bboxes[edge_idx];
          pnew.x = std::clamp(pnew.x, ibox.min_x, ibox.max_x);
          pnew.y = std::clamp(pnew.y, ibox.min_y, ibox.max_y);
        }
      }

      new_points.push_back(new_path[k]);
      new_points.push_back(pnew);
    }

    // if the path is not closed, ensure the last original point is added
    if (!new_path.is_closed()) new_points.push_back(new_path.back());

    // replace the original points with the resampled points
    new_path.clear();
    new_path.insert(new_path.end(), new_points.begin(), new_points.end());

    // update sigma by multiplying with persistence factor
    sigma *= persistence;
  }

  return new_path;
}

Path inflate(const Path &path, float radius, bool resample)
{
  if (!validate_min_size(path, 3, "Path points")) return path;

  Path new_path = path;

  std::vector<float>     curvature = new_path.get_curvature(true);
  std::vector<glm::vec2> normals = new_path.get_normals();

  for (size_t k = 0; k < new_path.size(); ++k)
  {
    float amp = -curvature[k] * radius;
    new_path[k].x += amp * normals[k].x;
    new_path[k].y += amp * normals[k].y;
  }

  // preserve the spatial resolution
  if (resample)
  {
    float dmax_prev = path.get_cumulative_distance().back();
    float dmax_new = new_path.get_cumulative_distance().back();
    float ratio = std::max(1.f, dmax_new / dmax_prev);

    int npts = int(path.size() * ratio);
    new_path.resample_interp(npts, InterpolationMethod1D::CUBIC);
  }

  return new_path;
}

Path meanderize(const Path            &path,
                float                  ratio,
                float                  noise_ratio,
                std::uint32_t          seed,
                int                    iterations,
                int                    edge_divisions,
                Path::EdgeDivisionMode edm)
{
  if (!validate_min_size(path, 2, "Path points")) return path;

  Path path_wrk = path;

  std::mt19937                    gen(seed);
  std::normal_distribution<float> dis(-noise_ratio, noise_ratio);

  for (int it = 0; it < iterations; it++)
  {
    Path new_path = Path();

    float cross_product;

    if (path.size() > 1)
      cross_product = (path[2].y - path[0].y) * (path[1].x - path[0].x) -
                      (path[2].x - path[0].x) * (path[1].y - path[0].y);
    else
      cross_product = 1.f;

    float orientation = -std::copysign(1.f, cross_product);

    size_t ks = path.is_closed() ? 0 : 1;
    for (size_t k = 0; k < path.size() - ks; k++)
    {
      size_t kp1 = (k + 1) % path.size();

      new_path.push_back(Point(path[k].x, path[k].y, path[k].v));

      float alpha = angle(path[kp1], path[k]);
      float dist = distance(path[kp1], path[k]);

      Point p = lerp(path[k], path[kp1], 0.5f);

      if (orientation >= 0.f)
        alpha += M_PI_2;
      else
        alpha -= M_PI_2;

      dist *= ratio * (1.f + dis(gen));

      p.x += dist * std::cos(alpha);
      p.y += dist * std::sin(alpha);

      new_path.push_back(p);
      orientation *= -1.f;
    }

    if (path.is_closed())
      new_path.push_back(path[0]);
    else
      new_path.push_back(path.back());

    path_wrk = new_path;
  }

  return bspline(path_wrk, edge_divisions, edm);
}

Array path_sdf_to_array(const Path  &path,
                        glm::ivec2   shape,
                        glm::vec4    bbox_array,
                        const Array *p_noise_x,
                        const Array *p_noise_y)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();

  Array array(shape);
  if (!validate_min_size(path, 2, "Path points")) return array;

  // array base grid
  std::vector<float> xg, yg;
  grid_xy_vector(xg, yg, shape, bbox_array, /* endpoint */ false);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float dx = p_noise_x ? (*p_noise_x)(i, j) : 0.f;
      float dy = p_noise_y ? (*p_noise_y)(i, j) : 0.f;
      float xi = xg[i] + dx;
      float yi = yg[j] + dy;

      // brute force
      float d = std::numeric_limits<float>::max();
      float s = 1.f;

      for (size_t i = 0, j = path.size() - 1; i < path.size(); j = i, i++)
      {
        glm::vec2 e = {path[j].x - path[i].x, path[j].y - path[i].y};
        glm::vec2 w = {xi - path[i].x, yi - path[i].y};
        float     coeff = std::clamp(dot(w, e) / dot(e, e), 0.f, 1.f);
        glm::vec2 b = {w.x - e.x * coeff, w.y - e.y * coeff};
        d = std::min(d, dot(b, b));

        std::array<bool, 3> c = {yi >= path[i].y,
                                 yi<path[j].y, e.x * w.y> e.y * w.x};

        if ((c[0] && c[1] && c[2]) || (not(c[0]) && not(c[1]) && not(c[2])))
          s *= -1.f;
      }

      // disccard in/out info for an open path
      if (!path.is_closed()) s = 1.f;

      array(i, j) = s * std::sqrt(d);
    }

  return array;
}

Path remove_geometric_loops(const Path &path)
{
  if (path.size() < 4) return path; // no loops possible

  Path result = path;
  bool changed;

  do
  {
    changed = false;
    for (size_t i = 0; i + 1 < result.size(); ++i)
    {
      for (size_t j = i + 2; j + 1 < result.size(); ++j)
      {
        auto inter = segment_intersection(result[i],
                                          result[i + 1],
                                          result[j],
                                          result[j + 1]);
        if (inter)
        {
          // remove points between segments i+1 and j
          result.erase(result.begin() + i + 1, result.begin() + j + 1);

          // insert intersection point
          result.insert(result.begin() + i + 1, *inter);
          changed = true;
          break;
        }
      }
      if (changed) break;
    }
  } while (changed);

  return result;
}

Path scale(const Path &path, glm::vec2 scale, glm::vec2 center)
{
  Path result = path;

  for (auto &p : result)
    p = hmap::scale(p, scale, center);

  return result;
}

Path scale(const Path &path, float scale_factor, glm::vec2 center)
{
  return scale(path, glm::vec2(scale_factor, scale_factor), center);
}

Path smooth(const Path &path,
            int         navg,
            float       averaging_intensity,
            float       inertia)
{
  Path new_path = path;

  if (!validate_non_empty(new_path, "Path points")) return new_path;

  const int  n = (int)new_path.size();
  const bool is_closed = new_path.is_closed();

  auto get_index = [&](int i) -> int
  {
    if (is_closed)
    {
      // wrap around
      i = i % n;
      if (i < 0) i += n;
      return i;
    }
    else
    {
      // clamp
      return std::max(0, std::min(i, n - 1));
    }
  };

  // --- Moving average

  std::vector<Point> smooth_points;
  smooth_points.reserve(n);

  for (int i = 0; i < n; i++)
  {
    int is;

    if (is_closed)
    {
      is = navg; // full window always valid
    }
    else
    {
      int i1 = std::min(navg, i);
      int i2 = std::min(navg, n - 1 - i);
      is = std::min(i1, i2);
    }

    Point psum(0.f, 0.f);

    for (int k = -is; k <= is; k++)
    {
      int idx = get_index(i + k);
      psum = psum + new_path[idx];
    }

    Point avg = psum / float(2 * is + 1);

    Point new_point = (1.f - averaging_intensity) * new_path[i] +
                      averaging_intensity * avg;

    smooth_points.push_back(new_point);
  }

  new_path.clear();
  new_path.insert(new_path.end(), smooth_points.begin(), smooth_points.end());

  // --- Inertia

  if (inertia > 0.f)
  {
    if (is_closed)
    {
      for (int i = 0; i < n; i++)
      {
        int prev = get_index(i - 1);
        new_path[i] = (1.f - inertia) * new_path[i] + inertia * new_path[prev];
      }
    }
    else
    {
      for (int i = 1; i < n - 1; i++)
      {
        new_path[i] = (1.f - inertia) * new_path[i] + inertia * new_path[i - 1];
      }
    }
  }

  return new_path;
}

} // namespace hmap
