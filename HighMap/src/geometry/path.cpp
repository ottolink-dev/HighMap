/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <list>
#include <string>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/colormaps.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate1d.hpp"
#include "highmap/math/core.hpp"
#include "highmap/morphology.hpp"
#include "highmap/operator.hpp"
#include "highmap/range.hpp"

namespace hmap
{

void Path::clear()
{
  *this = Path();
}

void Path::divide()
{
  size_t npoints = this->size();
  if (npoints == 0) return;

  size_t end = this->is_closed() ? npoints : npoints - 1;

  std::vector<Point> new_points;
  new_points.reserve(2 * npoints);

  for (size_t k = 0; k < end; ++k)
  {
    size_t knext = (k + 1) % npoints;

    // add the original point
    new_points.push_back(this->points[k]);

    // calculate and add the midpoint
    Point midpoint = lerp(this->points[k], this->points[knext], 0.5f);
    new_points.push_back(midpoint);
  }

  if (!this->is_closed()) new_points.push_back(this->points.back());

  this->points = std::move(new_points);
}

void Path::enforce_monotonic_values(bool decreasing)
{
  if (this->points.size() < 2) return;

  if (decreasing)
  {
    for (size_t k = 0; k < this->size() - 1; k++)
    {
      if (this->points[k + 1].v > this->points[k].v)
        this->points[k + 1].v = this->points[k].v;
    }
  }
  else
  {
    for (size_t k = 0; k < this->size() - 1; k++)
    {
      if (this->points[k + 1].v < this->points[k].v)
        this->points[k + 1].v = this->points[k].v;
    }
  }
}

std::vector<float> Path::get_arc_length() const
{
  if (this->empty()) return {};

  std::vector<float> s = this->get_cumulative_distance();
  float              total_len = s.back();

  if (total_len > 0.f)
  {
    for (auto &v : s)
      v /= total_len;
  }
  else
  {
    std::fill(s.begin(), s.end(), 0.f);
  }

  return s;
}

std::vector<float> Path::get_cumulative_distance() const
{
  if (this->empty()) return {};

  size_t             ke = this->is_closed() ? 1 : 0;
  std::vector<float> dacc(this->size() + ke, 0.f);

  for (size_t k = 1; k < this->size() + ke; k++)
  {
    size_t knext = k % this->size();
    float  dist = distance(this->points[k - 1], this->points[knext]);
    dacc[k] = dacc[k - 1] + dist;
  }

  return dacc;
}

std::vector<float> Path::get_curvature(bool normalized) const
{
  if (this->size() < 3) return std::vector<float>(this->size(), 0.f);

  size_t             ke = this->is_closed() ? 1 : 0;
  std::vector<float> cv(this->size() + ke, 0.f);
  float              cmax = 0.f;

  for (size_t k = 1; k < this->size() - 1 + ke; k++)
  {
    size_t km = (k - 1) % this->size();
    size_t kp = (k + 1) % this->size();
    cv[k] = curvature_signed(this->points[km],
                             this->points[k],
                             this->points[kp]);

    cmax = std::max(cmax, std::abs(cv[k]));
  }

  if (normalized && cmax > 0.f)
  {
    for (auto &v : cv)
      v /= cmax;
  }

  return cv;
}

std::vector<Point> Path::get_edge_centers() const
{
  std::vector<Point> centers;
  const size_t       npts = this->size();

  if (npts < 2) return centers;

  size_t ke = this->is_closed() ? 1 : 0;
  centers.reserve(npts - 1 + ke);

  for (size_t k = 0; k < npts - 1 + ke; ++k)
  {
    size_t kp = (k + 1) % npts;
    Point  pmid = lerp(this->points[k], this->points[kp], 0.5f);
    centers.emplace_back(pmid);
  }

  return centers;
}

std::vector<glm::vec2> Path::get_normals() const
{
  std::vector<glm::vec2> normals = this->get_tangents();

  for (auto &n : normals)
    n = glm::vec2(-n.y, n.x);

  return normals;
}

std::vector<glm::vec2> Path::get_tangents() const
{
  if (this->empty()) return {};
  if (this->size() == 1) return {glm::vec2(1.f, 0.f)};

  size_t                 ke = this->is_closed() ? 1 : 0;
  std::vector<glm::vec2> tangents(this->size() + ke, glm::vec2(0.f));

  for (size_t k = 1; k < this->size() - 1 + ke; k++)
  {
    size_t    km = (k - 1) % this->size();
    size_t    kp = (k + 1) % this->size();
    glm::vec2 delta = {this->points[kp].x - this->points[km].x,
                       this->points[kp].y - this->points[km].y};

    float len = glm::length(delta);
    tangents[k] = len > 0.f ? delta / len : glm::vec2(1.f, 0.f);
  }

  if (!this->is_closed())
  {
    glm::vec2 d_start = {this->points[1].x - this->points[0].x,
                         this->points[1].y - this->points[0].y};
    float     len_start = glm::length(d_start);
    tangents.front() = len_start > 0.f ? d_start / len_start
                                       : glm::vec2(1.f, 0.f);

    size_t    n = this->size();
    glm::vec2 d_end = {this->points[n - 1].x - this->points[n - 2].x,
                       this->points[n - 1].y - this->points[n - 2].y};
    float     len_end = glm::length(d_end);
    tangents.back() = len_end > 0.f ? d_end / len_end : glm::vec2(1.f, 0.f);
  }

  return tangents;
}

std::vector<float> Path::get_values() const
{
  std::vector<float> v(this->size());
  for (size_t i = 0; i < this->size(); i++)
    v[i] = this->points[i].v;

  if (this->is_closed() && this->size() > 0) v.push_back(this->points[0].v);
  return v;
}

std::vector<float> Path::get_x() const
{
  std::vector<float> x(this->size());
  for (size_t i = 0; i < this->size(); i++)
    x[i] = this->points[i].x;

  if (this->is_closed() && this->size() > 0) x.push_back(this->points[0].x);
  return x;
}

std::vector<float> Path::get_xy() const
{
  std::vector<float> xy(2 * this->size());
  for (size_t i = 0; i < this->size(); i++)
  {
    xy[2 * i] = this->points[i].x;
    xy[2 * i + 1] = this->points[i].y;
  }

  if (this->is_closed() && this->size() > 0)
  {
    xy.push_back(this->points[0].x);
    xy.push_back(this->points[0].y);
  }
  return xy;
}

std::vector<float> Path::get_y() const
{
  std::vector<float> y(this->size());
  for (size_t i = 0; i < this->size(); i++)
    y[i] = this->points[i].y;

  if (this->is_closed() && this->size() > 0) y.push_back(this->points[0].y);
  return y;
}

bool Path::is_closed() const
{
  return this->path_closure == PathClosure::PT_CLOSE;
}

void Path::reorder_nns(int start_index)
{
  const size_t n = this->size();
  if (n <= 1) return;

  start_index = std::clamp(start_index, 0, static_cast<int>(n - 1));

  std::vector<Point> reordered;
  reordered.reserve(n);

  std::vector<bool> visited(n, false);
  int               current = start_index;
  visited[current] = true;
  reordered.push_back(this->points[current]);

  for (size_t step = 1; step < n; ++step)
  {
    int   best_idx = -1;
    float best_dist_sq = std::numeric_limits<float>::max();

    const Point &p_curr = this->points[current];

    for (size_t i = 0; i < n; ++i)
    {
      if (visited[i]) continue;

      float dx = this->points[i].x - p_curr.x;
      float dy = this->points[i].y - p_curr.y;
      float dist_sq = dx * dx + dy * dy;

      if (dist_sq < best_dist_sq)
      {
        best_dist_sq = dist_sq;
        best_idx = static_cast<int>(i);
      }
    }

    if (best_idx != -1)
    {
      visited[best_idx] = true;
      reordered.push_back(this->points[best_idx]);
      current = best_idx;
    }
  }

  this->points = std::move(reordered);
}

void Path::resample_by_spacing(float delta, InterpolationMethod1D itp_method)
{
  if (!validate_min_size(this->points, 2, "Path points")) return;

  std::vector<float> cdist = this->get_cumulative_distance();
  int                npoints = std::max(2, int(cdist.back() / delta));

  this->resample_interp(npoints, itp_method);
}

void Path::resample_interp(int npoints, InterpolationMethod1D itp_method)
{
  if (!validate_min_size(this->points, 2, "Path points")) return;

  // work on a copy to manage open/close paths
  Path path_wrk = *this;

  // duplicate 1st/last point for closed path
  if (this->is_closed()) path_wrk.push_back(path_wrk.front());

  // interpolation points along arc
  std::vector<float> t = hmap::linspace(0.f, 1.f, npoints);

  // interpolation functions
  std::vector<float> arc = path_wrk.get_arc_length();
  std::vector<float> x = path_wrk.get_x();
  std::vector<float> y = path_wrk.get_y();
  std::vector<float> v = path_wrk.get_values();

  Interpolator1D itp_x = hmap::Interpolator1D(arc, x, itp_method);
  Interpolator1D itp_y = hmap::Interpolator1D(arc, y, itp_method);
  Interpolator1D itp_v = hmap::Interpolator1D(arc, v, itp_method);

  // interpolate
  std::vector<Point> new_points(npoints);

  for (size_t k = 0; k < size_t(npoints); ++k)
  {
    float ti = t[k];
    Point p(itp_x(ti), itp_y(ti), itp_v(ti));
    new_points[k] = p;
  }

  // remove duplicate 1st/last point
  if (this->is_closed()) new_points.pop_back();

  this->points = std::move(new_points);
}

void Path::resample_uniform(InterpolationMethod1D itp_method)
{
  if (!validate_min_size(this->points, 2, "Path points")) return;

  float dmin = std::numeric_limits<float>::max();
  float dsum = 0.f;

  size_t npoints = this->size();
  size_t end = this->is_closed() ? npoints : npoints - 1;

  for (size_t k = 0; k < end; k++)
  {
    size_t knext = (k + 1) % npoints;
    float  dist = distance(this->points[k], this->points[knext]);
    if (dist < dmin) dmin = dist;

    dsum += dist;
  }

  if (dmin <= 0.f) return;

  this->resample_interp(static_cast<int>(dsum / dmin), itp_method);
}

void Path::reverse()
{
  std::reverse(this->points.begin(), this->points.end());
}

glm::vec3 Path::sample_at(float                     t,
                          const std::vector<float> *p_arc,
                          glm::vec2                *p_tangent) const
{
  if (this->empty())
  {
    if (p_tangent) *p_tangent = glm::vec2(1.f, 0.f);
    return glm::vec3(0.f);
  }

  if (this->size() == 1)
  {
    if (p_tangent) *p_tangent = glm::vec2(1.f, 0.f);
    return glm::vec3(this->points[0].x, this->points[0].y, this->points[0].v);
  }

  const std::vector<float>  arc = p_arc ? *p_arc : this->get_arc_length();
  const std::vector<Point> &pts = this->points;

  t = std::clamp(t, 0.f, 1.f);

  size_t k = 1;
  while (k < arc.size() - 1 && arc[k] < t)
    k++;

  float span = std::max(arc[k] - arc[k - 1], 1e-9f);
  float r = (t - arc[k - 1]) / span;

  const Point &p0 = pts[k - 1];
  const Point &p1 = pts[k];

  glm::vec2 pos = {(1.f - r) * p0.x + r * p1.x, (1.f - r) * p0.y + r * p1.y};
  float     value = lerp(p0.v, p1.v, r);

  if (p_tangent)
  {
    float ddx = p1.x - p0.x;
    float ddy = p1.y - p0.y;
    float dn = std::hypot(ddx, ddy);
    *p_tangent = dn > 0.f ? glm::vec2(ddx / dn, ddy / dn) : glm::vec2(1.f, 0.f);
  }

  return glm::vec3(pos.x, pos.y, value);
}

void Path::set_closed(bool new_value)
{
  this->path_closure = new_value ? PathClosure::PT_CLOSE : PathClosure::PT_OPEN;
}

void Path::subsample(int step)
{
  if (step <= 1 || this->size() <= 2) return;

  std::vector<Point> subsampled;
  subsampled.reserve((this->size() + step - 1) / step + 1);

  for (size_t k = 0; k < this->size(); k += step)
  {
    subsampled.push_back(this->points[k]);
  }

  // always preserve endpoint for open paths
  if (!this->is_closed() && (this->size() - 1) % step != 0)
  {
    subsampled.push_back(this->points.back());
  }

  this->points = std::move(subsampled);
}

void Path::to_array(Array &array, glm::vec4 bbox, bool filled) const
{
  if (!validate_non_empty(array)) return;

  float lx = bbox.y - bbox.x;
  float ly = bbox.w - bbox.z;
  float ppu = std::max(array.shape.x / lx, array.shape.y / ly);

  // create a temporary cloud with the right points density (= 1 ppu)
  Cloud cloud = Cloud(points);

  // project path itself
  size_t ks = this->is_closed() ? 0 : 1;
  for (size_t k = 0; k < this->size() - ks; k++)
  {
    size_t knext = (k + 1) % this->size();
    int    npixels = (int)std::ceil(
                      distance(this->points[k], this->points[knext]) * ppu) +
                  1;

    for (int i = 0; i < npixels; i++)
    {
      float t = (float)i / (float)(npixels - 1);
      Point p = lerp(this->points[k], this->points[knext], t);
      cloud.push_back(p);
    }
  }

  // if filled, set the border to the same value of the filling value
  if (filled) cloud.set_values(1.f);

  cloud.to_array(array, bbox);

  // flood filling
  if (filled)
  {
    Array array_bckp = array;

    // TODO make something more robust
    int i = 0;
    int j = 0;

    flood_fill(array, i, j);

    array = 1.f - array;
    array = maximum(array, array_bckp);
  }
}

Array Path::to_array(glm::ivec2 shape, glm::vec4 bbox, bool filled) const
{
  if (!validate_shape(shape)) return Array();

  Array array(shape);
  this->to_array(array, bbox, filled);
  return array;
}

void Path::to_array_mask(Array &array, glm::vec4 bbox, bool filled) const
{
  if (!validate_non_empty(array)) return;

  Path path_copy = *this;
  path_copy.set_values(1.f);
  path_copy.to_array(array, bbox, filled);
}

void Path::to_png(std::string fname, glm::ivec2 shape)
{
  if (!validate_shape(shape)) return;

  Array array = Array(shape);
  this->to_array(array, this->get_bbox());
  array.to_png(fname, Cmap::INFERNO, false);
}

} // namespace hmap
