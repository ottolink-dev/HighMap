#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <optional>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

// ==========================================================================
//  Point Methods
// ==========================================================================

void Point::print() const
{
  std::cout << "(" << this->x << ", " << this->y << ", " << this->v << ")"
            << std::endl;
}

void Point::set_value_from_array(const Array &array, const glm::vec4 &bbox)
{
  if (!validate_non_empty(array))
  {
    this->v = 0.f;
    return;
  }

  float dx = bbox.y - bbox.x;
  float dy = bbox.w - bbox.z;
  if (std::abs(dx) < 1e-7f || std::abs(dy) < 1e-7f)
  {
    this->v = 0.f;
    return;
  }

  // scale to unit interval
  float xn = (this->x - bbox.x) / dx;
  float yn = (this->y - bbox.z) / dy;

  // scale to array shape
  xn *= static_cast<float>(array.shape.x - 1);
  yn *= static_cast<float>(array.shape.y - 1);

  int i = static_cast<int>(xn);
  int j = static_cast<int>(yn);

  if (i >= 0 && i < array.shape.x && j >= 0 && j < array.shape.y)
  {
    float uu = xn - static_cast<float>(i);
    float vv = yn - static_cast<float>(j);
    this->v = array.get_value_bilinear_at(i, j, uu, vv);
  }
  else
  {
    this->v = 0.f; // if outside array bounding box
  }
}

// ==========================================================================
//  Free Functions
// ==========================================================================

float angle(const Point &p1, const Point &p2)
{
  return std::atan2(p2.y - p1.y, p2.x - p1.x);
}

float angle(const Point &p0, const Point &p1, const Point &p2)
{
  // calculate vectors
  float dx1 = p1.x - p0.x;
  float dy1 = p1.y - p0.y;
  float dx2 = p2.x - p0.x;
  float dy2 = p2.y - p0.y;

  // compute the angle using the atan2 function to get the signed angle
  float angle1 = std::atan2(dy1, dx1);
  float angle2 = std::atan2(dy2, dx2);

  // calculate the angle difference
  float ang = angle2 - angle1;

  // normalize the angle to be in the range [-π, π]
  if (ang > M_PI)
    ang -= 2.f * float(M_PI);
  else if (ang < -M_PI)
    ang += 2.f * float(M_PI);

  return ang;
}

float classify_point(const Point &p_prev,
                     const Point &p,
                     const Point &p_next,
                     const Point &pq)
{
  float c = curvature_signed(p_prev, p, p_next);
  float s = side(p_prev, p, p_next, pq);

  float val = c * s;

  if (val > 0.f) return 1.f;  // convex
  if (val < 0.f) return -1.f; // concave
  return 0.f;
}

float cross_product(const Point &p0, const Point &p1, const Point &p2)
{
  return (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x);
}

float cross_product(const Point &p1, const Point &p2)
{
  return p1.x * p2.y - p1.y * p2.x;
}

float curvature(const Point &p1, const Point &p2, const Point &p3)
{
  float d12_sq = distance_squared(p1, p2);
  float d23_sq = distance_squared(p2, p3);
  float d13_sq = distance_squared(p1, p3);
  float d_prod = d12_sq * d23_sq * d13_sq;

  if (d_prod > 1e-24f)
    return 4.f * triangle_area(p1, p2, p3) / std::sqrt(d_prod);
  else
    return 0.f;
}

float curvature_signed(const Point &p1, const Point &p2, const Point &p3)
{
  float d12_sq = distance_squared(p1, p2);
  float d23_sq = distance_squared(p2, p3);
  float d13_sq = distance_squared(p1, p3);
  float d_prod = d12_sq * d23_sq * d13_sq;

  if (d_prod > 1e-24f)
    return 4.f * triangle_area_signed(p1, p2, p3) / std::sqrt(d_prod);
  else
    return 0.f;
}

float distance(const Point &p1, const Point &p2)
{
  return std::hypot(p1.x - p2.x, p1.y - p2.y);
}

float distance_squared(const Point &p1, const Point &p2)
{
  float dx = p1.x - p2.x;
  float dy = p1.y - p2.y;
  return dx * dx + dy * dy;
}

Point interp_bezier(const Point &p_start,
                    const Point &p_ctrl_start,
                    const Point &p_ctrl_end,
                    const Point &p_end,
                    float        t)
{
  float u = 1.f - t;
  float tt = t * t;
  float uu = u * u;
  float uuu = uu * u;
  float ttt = tt * t;

  float c0 = uuu;
  float c1 = 3.f * uu * t;
  float c2 = 3.f * u * tt;
  float c3 = ttt;

  return Point(
      c0 * p_start.x + c1 * p_ctrl_start.x + c2 * p_ctrl_end.x + c3 * p_end.x,
      c0 * p_start.y + c1 * p_ctrl_start.y + c2 * p_ctrl_end.y + c3 * p_end.y,
      c0 * p_start.v + c1 * p_ctrl_start.v + c2 * p_ctrl_end.v + c3 * p_end.v);
}

Point interp_bspline(const Point &p0,
                     const Point &p1,
                     const Point &p2,
                     const Point &p3,
                     float        t)
{
  float tt = t * t;
  float ttt = tt * t;

  float c0 = (-ttt + 3.f * tt - 3.f * t + 1.f) / 6.f;
  float c1 = (3.f * ttt - 6.f * tt + 4.f) / 6.f;
  float c2 = (-3.f * ttt + 3.f * tt + 3.f * t + 1.f) / 6.f;
  float c3 = ttt / 6.f;

  return Point(c0 * p0.x + c1 * p1.x + c2 * p2.x + c3 * p3.x,
               c0 * p0.y + c1 * p1.y + c2 * p2.y + c3 * p3.y,
               c0 * p0.v + c1 * p1.v + c2 * p2.v + c3 * p3.v);
}

Point interp_catmullrom(const Point &p0,
                        const Point &p1,
                        const Point &p2,
                        const Point &p3,
                        float        t)
{
  float tt = t * t;
  float ttt = tt * t;

  float c0 = (-ttt + 2.f * tt - t) * 0.5f;
  float c1 = (3.f * ttt - 5.f * tt + 2.f) * 0.5f;
  float c2 = (-3.f * ttt + 4.f * tt + t) * 0.5f;
  float c3 = (ttt - tt) * 0.5f;

  return Point(c0 * p0.x + c1 * p1.x + c2 * p2.x + c3 * p3.x,
               c0 * p0.y + c1 * p1.y + c2 * p2.y + c3 * p3.y,
               c0 * p0.v + c1 * p1.v + c2 * p2.v + c3 * p3.v);
}

Point interp_decasteljau(const std::vector<Point> &points, float t)
{
  if (points.empty()) return Point{};
  if (points.size() == 1) return points[0];

  std::vector<Point> temp = points;
  float              u = 1.f - t;
  for (size_t n = temp.size(); n > 1; --n)
  {
    for (size_t i = 0; i < n - 1; ++i)
    {
      temp[i].x = temp[i].x * u + temp[i + 1].x * t;
      temp[i].y = temp[i].y * u + temp[i + 1].y * t;
      temp[i].v = temp[i].v * u + temp[i + 1].v * t;
    }
  }
  return temp[0];
}

glm::vec4 intersect_bounding_boxes(const glm::vec4 &bbox1,
                                   const glm::vec4 &bbox2)
{
  float min_x = std::max(bbox1.x, bbox2.x);
  float max_x = std::min(bbox1.y, bbox2.y);
  float min_y = std::max(bbox1.z, bbox2.z);
  float max_y = std::min(bbox1.w, bbox2.w);

  if (min_x <= max_x && min_y <= max_y)
  {
    return glm::vec4{min_x, max_x, min_y, max_y};
  }

  // else return an "impossible" bounding box with xmin > xmax and ymin > ymax
  return glm::vec4(1.f, -1.f, 1.f, -1.f);
}

bool is_point_within_bounding_box(const Point &p, const glm::vec4 &bbox)
{
  return p.x >= bbox.x && p.x <= bbox.y && p.y >= bbox.z && p.y <= bbox.w;
}

bool is_point_within_bounding_box(float x, float y, const glm::vec4 &bbox)
{
  return x >= bbox.x && x <= bbox.y && y >= bbox.z && y <= bbox.w;
}

Point lerp(const Point &p1, const Point &p2, float t)
{
  return Point(p1.x + t * (p2.x - p1.x),
               p1.y + t * (p2.y - p1.y),
               p1.v + t * (p2.v - p1.v));
}

Point midpoint(const Point &p1,
               const Point &p2,
               int          orientation,
               float        distance_ratio,
               float        t)
{
  Point pmid = lerp(p1, p2, t);

  distance_ratio = orientation == 0 ? distance_ratio
                                    : std::abs(distance_ratio) *
                                          std::copysign(1.0f, orientation);

  float dx = p2.x - p1.x;
  float dy = p2.y - p1.y;

  float perp_x = -dy * distance_ratio;
  float perp_y = dx * distance_ratio;

  return Point(pmid.x + perp_x, pmid.y + perp_y, pmid.v);
}

Point scale(const Point &p, glm::vec2 scale_factors, glm::vec2 center)
{
  return Point(center.x + scale_factors.x * (p.x - center.x),
               center.y + scale_factors.y * (p.y - center.y),
               p.v);
}

Point scale(const Point &p, float scale_factor, glm::vec2 center)
{
  return scale(p, glm::vec2(scale_factor, scale_factor), center);
}

std::optional<Point> segment_intersection(const Point &p1,
                                          const Point &p2,
                                          const Point &q1,
                                          const Point &q2)
{
  float rx = p2.x - p1.x;
  float ry = p2.y - p1.y;
  float sx = q2.x - q1.x;
  float sy = q2.y - q1.y;

  float rxs = rx * sy - ry * sx;
  if (std::abs(rxs) < 1e-10f) return std::nullopt; // parallel or collinear

  float qpx = q1.x - p1.x;
  float qpy = q1.y - p1.y;

  float t = (qpx * sy - qpy * sx) / rxs;
  float u = (qpx * ry - qpy * rx) / rxs;

  if (t >= 0.f && t <= 1.f && u >= 0.f && u <= 1.f)
  {
    return Point(p1.x + t * rx, p1.y + t * ry, p1.v + t * (p2.v - p1.v));
  }

  return std::nullopt; // no intersection
}

float side(const Point &p1,
           const Point &p2,
           const Point &p3,
           const Point &p_query)
{
  // tangent direction at p2 (no need to normalize — only sign matters)
  float tx = p3.x - p1.x;
  float ty = p3.y - p1.y;

  // vector from p2 to query point
  float qx = p_query.x - p2.x;
  float qy = p_query.y - p2.y;

  // z-component of cross product T × Q
  float cross = tx * qy - ty * qx;

  if (cross > 0.f)
    return 1.f; // left of curve (CCW side)
  else if (cross < 0.f)
    return -1.f; // right of curve (CW side)
  else
    return 0.f; // on the tangent line
}

static bool cmp_inf(const Point &a, const Point &b) noexcept
{
  if (a.x != b.x) return a.x < b.x;
  if (a.y != b.y) return a.y < b.y;
  return a.v < b.v;
}

void sort_points(std::vector<Point> &points)
{
  std::sort(points.begin(), points.end(), cmp_inf);
}

float triangle_area(const Point &p1, const Point &p2, const Point &p3)
{
  return std::abs(triangle_area_signed(p1, p2, p3));
}

float triangle_area_signed(const Point &p1, const Point &p2, const Point &p3)
{
  return 0.5f *
         (p1.x * (p2.y - p3.y) + p2.x * (p3.y - p1.y) + p3.x * (p1.y - p2.y));
}

glm::vec4 unit_square_bbox()
{
  return glm::vec4(0.f, 1.f, 0.f, 1.f);
}

} // namespace hmap
