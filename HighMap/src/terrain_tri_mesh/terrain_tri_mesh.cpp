/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <glm/geometric.hpp>

#include "delaunator-cpp.hpp"
#include "hmm/src/heightmap.h"
#include "hmm/src/triangulator.h"

#include "highmap/array.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate2d.hpp"
#include "highmap/terrain_tri_mesh.hpp"

#include <unordered_map>

namespace hmap
{

// ------------------------------------------------------------
// BoundingBox
// ------------------------------------------------------------

bool TerrainTriMesh::BoundingBox::contains(const glm::vec2 &p) const
{
  return (p.x >= this->min.x && p.x <= this->max.x && p.y >= this->min.y &&
          p.y <= this->max.y);
}

glm::vec2 TerrainTriMesh::BoundingBox::clamp(const glm::vec2 &p) const
{
  return {std::clamp(p.x, this->min.x, this->max.x),
          std::clamp(p.y, this->min.y, this->max.y)};
}

// ------------------------------------------------------------
// TerrainTriMesh
// ------------------------------------------------------------

TerrainTriMesh::TerrainTriMesh(const std::vector<glm::vec3> &ref_points)
    : points(ref_points)
{
  this->triangulate_delaunay();
  this->compute_neighbors();
}

TerrainTriMesh::TerrainTriMesh(const std::vector<float> &x,
                               const std::vector<float> &y,
                               const std::vector<float> &z)
{
  this->points.reserve(x.size());

  for (size_t k = 0; k < x.size(); ++k)
    this->points.push_back(glm::vec3(x[k], y[k], z[k]));

  this->triangulate_delaunay();
  this->compute_neighbors();
}

bool TerrainTriMesh::barycentric(const glm::vec2 &p,
                                 size_t           i0,
                                 size_t           i1,
                                 size_t           i2,
                                 float           &w0,
                                 float           &w1,
                                 float           &w2) const
{
  float x0 = this->points[i0].x, y0 = this->points[i0].y;
  float x1 = this->points[i1].x, y1 = this->points[i1].y;
  float x2 = this->points[i2].x, y2 = this->points[i2].y;

  float det = (x1 - x0) * (y2 - y0) - (y1 - y0) * (x2 - x0);
  if (std::abs(det) < 1e-8f) return false;

  w0 = (p.x - x1) * (y2 - y1) - (p.y - y1) * (x2 - x1);
  w1 = (p.x - x2) * (y0 - y2) - (p.y - y2) * (x0 - x2);
  w2 = (p.x - x0) * (y1 - y0) - (p.y - y0) * (x1 - x0);

  const float eps = -1e-6f;
  if (w0 < eps || w1 < eps || w2 < eps) return false;

  float inv = -1.f / det;
  w0 *= inv;
  w1 *= inv;
  w2 *= inv;

  return true;
}

void TerrainTriMesh::compute_gradients()
{
  this->gradients.clear();
  this->gradients.resize(points.size(), {0.f, 0.f});
  std::vector<float> weight_sum(points.size(), 0.f);

  auto angle_at =
      [](const glm::vec2 &a, const glm::vec2 &b, const glm::vec2 &c) -> float
  {
    glm::vec2 u = glm::normalize(b - a);
    glm::vec2 v = glm::normalize(c - a);

    float dot = glm::clamp(glm::dot(u, v), -1.0f, 1.0f);
    return std::acos(dot);
  };

  for (const auto &tri : triangles)
  {
    const auto &p0 = points[tri.a];
    const auto &p1 = points[tri.b];
    const auto &p2 = points[tri.c];

    glm::vec2 a(p0.x, p0.y);
    glm::vec2 b(p1.x, p1.y);
    glm::vec2 c(p2.x, p2.y);

    float dx1 = p1.x - p0.x;
    float dy1 = p1.y - p0.y;
    float dz1 = p1.z - p0.z;

    float dx2 = p2.x - p0.x;
    float dy2 = p2.y - p0.y;
    float dz2 = p2.z - p0.z;

    float det = dx1 * dy2 - dx2 * dy1;
    if (std::abs(det) < 1e-12f) continue;

    // Triangle gradient (constant over triangle)
    float gx = (dz1 * dy2 - dz2 * dy1) / det;
    float gy = (dx1 * dz2 - dx2 * dz1) / det;

    // Compute angles at each vertex
    float angle0 = angle_at(a, b, c);
    float angle1 = angle_at(b, c, a);
    float angle2 = angle_at(c, a, b);

    auto accumulate = [&](size_t i, float w)
    {
      this->gradients[i].x += w * gx;
      this->gradients[i].y += w * gy;
      weight_sum[i] += w;
    };

    accumulate(tri.a, angle0);
    accumulate(tri.b, angle1);
    accumulate(tri.c, angle2);
  }

  // Normalize
  for (size_t i = 0; i < points.size(); ++i)
  {
    if (weight_sum[i] > 1e-12f)
    {
      this->gradients[i].x /= weight_sum[i];
      this->gradients[i].y /= weight_sum[i];
    }
  }
}

void TerrainTriMesh::compute_neighbors()
{
  this->neighbors.clear();
  this->neighbors.resize(this->points.size());

  const auto &pts = this->points;

  auto add_edge = [&](size_t a, size_t b)
  {
    // check if already exists (avoid duplicates)
    for (const auto &nb : this->neighbors.adjacency[a])
    {
      if (nb.index == b) return;
    }

    float vx = pts[a].x - pts[b].x;
    float vy = pts[a].y - pts[b].y;
    float dist = glm::length(glm::vec2(vx, vy));

    this->neighbors.adjacency[a].push_back({b, dist});
    this->neighbors.adjacency[b].push_back({a, dist});
  };

  for (const auto &tri : this->triangles)
  {
    add_edge(tri.a, tri.b);
    add_edge(tri.b, tri.c);
    add_edge(tri.c, tri.a);
  }
}

bool TerrainTriMesh::export_obj(const std::string &filepath) const
{
  std::ofstream file(filepath);
  if (!file.is_open()) return false;

  // vertices
  for (const auto &p : this->points)
  {
    file << "v " << p.x << " " << p.y << " " << p.z << "\n";
  }

  // faces (OBJ is 1-based indexing)
  for (const auto &tri : this->triangles)
  {
    file << "f " << (tri.a + 1) << " " << (tri.b + 1) << " " << (tri.c + 1)
         << "\n";
  }

  return true;
}

int TerrainTriMesh::find_triangle(const glm::vec2 &p,
                                  int              start_tri,
                                  bool             linear_search) const
{
  int       tri = std::max(0, start_tri);
  const int max_iter = int(points.size());

  if (linear_search)
  {
    // fallback for debugging purposes: linear search
    for (size_t k = 0; k < triangles.size(); ++k)
    {
      float       w0, w1, w2;
      const auto &t = triangles[k];
      if (barycentric(p, t.a, t.b, t.c, w0, w1, w2))
      {
        if (w0 >= 0 && w1 >= 0 && w2 >= 0) return int(k);
      }
    }
  }
  else
  {

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"

    for (int iter = 0; iter < max_iter; ++iter)
    {
      const auto &t = triangles[tri];
      float       w0, w1, w2;
      this->barycentric(p, t.a, t.b, t.c, w0, w1, w2);

      if (w0 >= 0 && w1 >= 0 && w2 >= 0) return tri;

      // walk to neighbor across the negative weight
      if (w0 < 0.f)
        tri = this->neighbor_triangle(tri, 1);
      else if (w1 < 0.f)
        tri = this->neighbor_triangle(tri, 2);
      else
        tri = this->neighbor_triangle(tri, 0);

      if (tri < 0) return -1;
    }

#pragma GCC diagnostic pop
  }

  return -1;
}

TerrainTriMesh::BoundingBox TerrainTriMesh::get_bbox() const
{
  glm::vec3 pmin(std::numeric_limits<float>::max());
  glm::vec3 pmax(std::numeric_limits<float>::lowest());

  for (const auto &p : this->points)
  {
    pmin = glm::min(pmin, p);
    pmax = glm::max(pmax, p);
  }

  TerrainTriMesh::BoundingBox bbox;
  bbox.min = pmin;
  bbox.max = pmax;

  return bbox;
}

const std::vector<size_t> &TerrainTriMesh::get_convex_hull() const
{
  return this->convex_hull;
}

const TerrainTriMesh::NeighborData &TerrainTriMesh::get_neighbors() const
{
  if (this->neighbors.adjacency.empty())
    throw std::runtime_error(
        "TerrainTriMesh::NeighborData: neighbors data has not been "
        "initialized, make sure TerrainTriMesh::compute_neighbors has been "
        "called before-hand.");

  return this->neighbors;
}

const std::vector<glm::vec3> &TerrainTriMesh::get_points() const
{
  return this->points;
}

std::vector<glm::vec3> &TerrainTriMesh::get_points()
{
  return this->points;
}

glm::vec2 TerrainTriMesh::get_range_z() const
{
  if (this->points.empty()) return glm::vec2(0.f);

  float zmin = std::numeric_limits<float>::max();
  float zmax = std::numeric_limits<float>::lowest();

  for (const auto &p : this->points)
  {
    if (p.z < zmin) zmin = p.z;
    if (p.z > zmax) zmax = p.z;
  }

  return glm::vec2(zmin, zmax);
}

float TerrainTriMesh::get_reference_area() const
{
  float total = 0.0f;

  for (const auto &tri : this->triangles)
  {
    const glm::vec3 &p0 = this->points[tri.a];
    const glm::vec3 &p1 = this->points[tri.b];
    const glm::vec3 &p2 = this->points[tri.c];

    glm::vec3 e0 = p1 - p0;
    glm::vec3 e1 = p2 - p0;

    float area = 0.5f * glm::length(glm::cross(e0, e1));

    total += area;
  }

  return total;
}

float TerrainTriMesh::get_reference_area_xy() const
{
  float total = 0.0f;

  for (const auto &tri : this->triangles)
  {
    const glm::vec2 p0 = this->to_xy(this->points[tri.a]);
    const glm::vec2 p1 = this->to_xy(this->points[tri.b]);
    const glm::vec2 p2 = this->to_xy(this->points[tri.c]);

    float area = 0.5f * std::abs((p1.x - p0.x) * (p2.y - p0.y) -
                                 (p2.x - p0.x) * (p1.y - p0.y));

    total += area;
  }

  return total;
}

glm::vec3 TerrainTriMesh::get_reference_lengths() const
{
  glm::vec3 total(0.0f);
  int       count = 0;

  for (const auto &tri : this->triangles)
  {
    auto accumulate = [&](const glm::vec3 &a, const glm::vec3 &b)
    {
      glm::vec3 d = glm::abs(b - a);
      total += d;
      count++;
    };

    accumulate(this->points[tri.a], this->points[tri.b]);
    accumulate(this->points[tri.b], this->points[tri.c]);
    accumulate(this->points[tri.c], this->points[tri.a]);
  }

  if (count > 0) total /= static_cast<float>(count);

  return total;
}

const std::vector<TerrainTriMesh::Triangle> &TerrainTriMesh::get_triangles()
    const
{
  return this->triangles;
}

std::vector<float> TerrainTriMesh::get_vertex_areas(bool normalized) const
{
  std::vector<float> areas(this->points.size(), 0.0f);

  glm::vec3 ref(1.f);

  if (normalized) ref = this->get_reference_lengths();

  for (const auto &tri : this->triangles)
  {
    const auto &p0 = this->points[tri.a];
    const auto &p1 = this->points[tri.b];
    const auto &p2 = this->points[tri.c];

    glm::vec3 e0 = p1 - p0;
    glm::vec3 e1 = p2 - p0;

    if (normalized)
    {
      // anisotropic normalization
      e0 /= ref;
      e1 /= ref;
    }

    float triangle_area = 0.5f * glm::length(glm::cross(e0, e1));
    float contribution = triangle_area / 3.f;

    areas[tri.a] += contribution;
    areas[tri.b] += contribution;
    areas[tri.c] += contribution;
  }

  return areas;
}

std::string TerrainTriMesh::info_string() const
{
  std::ostringstream oss;

  oss << "---- TerrainTriMesh Info ----\n";

  // --- counts ---
  const size_t nb_vertices = this->points.size();
  const size_t nb_triangles = this->triangles.size();

  oss << "vertices   : " << nb_vertices << "\n";
  oss << "triangles  : " << nb_triangles << "\n";

  if (this->points.empty())
  {
    oss << "empty mesh\n";
    return oss.str();
  }

  // --- bounding box ---
  glm::vec3 pmin(std::numeric_limits<float>::max());
  glm::vec3 pmax(std::numeric_limits<float>::lowest());

  for (const auto &p : this->points)
  {
    pmin = glm::min(pmin, p);
    pmax = glm::max(pmax, p);
  }

  glm::vec3 extent = pmax - pmin;

  oss << "min        : (" << pmin.x << ", " << pmin.y << ", " << pmin.z
      << ")\n";
  oss << "max        : (" << pmax.x << ", " << pmax.y << ", " << pmax.z
      << ")\n";
  oss << "extent     : (" << extent.x << ", " << extent.y << ", " << extent.z
      << ")\n";

  // --- derived metrics ---
  if (!this->triangles.empty())
  {
    float     area_xy = this->get_reference_area_xy();
    float     area_3d = this->get_reference_area();
    glm::vec3 ref_len = this->get_reference_lengths();

    oss << "area_xy    : " << area_xy << "\n";
    oss << "area_3d    : " << area_3d << "\n";
    oss << "roughness  : " << (area_xy > 0.f ? area_3d / area_xy : 0.f) << "\n";

    oss << "ref_len    : (" << ref_len.x << ", " << ref_len.y << ", "
        << ref_len.z << ")\n";
  }

  oss << "-----------------------------\n";

  return oss.str();
}

float TerrainTriMesh::interpolate_z_linear(const glm::vec2 &p,
                                           int             &last_tri,
                                           float            fill_value) const
{
  int tri = this->find_triangle(p, last_tri);

  if (tri >= 0)
  {
    last_tri = tri;
    const auto &t = this->triangles[tri];

    float w0, w1, w2;
    if (this->barycentric(p, t.a, t.b, t.c, w0, w1, w2))
    {
      // interpolate z using barycentric weights
      return w0 * this->points[t.a].z + w1 * this->points[t.b].z +
             w2 * this->points[t.c].z;
    }
  }

  // outside mesh or degenerate triangle
  return fill_value;
}

float TerrainTriMesh::interpolate_z_linear_gradient(
    const glm::vec2 &p,
    int             &last_tri,
    float            fill_value,
    float            gradient_scaling) const
{
  int tri = this->find_triangle(p, last_tri);

  if (tri >= 0)
  {
    last_tri = tri;
    const auto &t = this->triangles[tri];

    float w0, w1, w2;
    bool  ret = this->barycentric(p, t.a, t.b, t.c, w0, w1, w2);

    if (!ret) return fill_value;

    auto eval = [&](int idx, const glm::vec3 &vertex, float w)
    {
      glm::vec2 grad = gradient_scaling * this->gradients[idx];
      glm::vec2 d = p - glm::vec2(vertex);
      return w * (vertex.z + glm::dot(grad, d));
    };

    const auto &pa = this->points[t.a];
    const auto &pb = this->points[t.b];
    const auto &pc = this->points[t.c];

    float value = 0.0f;

    value += eval(t.a, pa, w0);
    value += eval(t.b, pb, w1);
    value += eval(t.c, pc, w2);

    return value;
  }

  // outside mesh or degenerate triangle
  return fill_value;
}

float TerrainTriMesh::interpolate_z_nearest(const glm::vec2 &p) const
{
  // "true" nearest as opposed to "nearest_approx" method
  float best_d2 = std::numeric_limits<float>::max();
  float best_z = 0.f;

  for (const auto &pt : points)
  {
    glm::vec2 delta = p - to_xy(pt);
    float     d2 = glm::dot(delta, delta);

    if (d2 < best_d2)
    {
      best_d2 = d2;
      best_z = pt.z;
    }
  }

  return best_z;
}

float TerrainTriMesh::interpolate_z_nearest_approx(const glm::vec2 &p,
                                                   int             &last_tri,
                                                   float fill_value) const
{
  // not true global nearest neighbor, pick only the nearest among the
  // 3 vertices of the containing triangle. Usually this is fine (and
  // fast), but: near triangle edges, the actual nearest vertex might
  // belong to a neighboring triangle and fials outside mesh
  int tri = this->find_triangle(p, last_tri);

  if (tri >= 0)
  {
    last_tri = tri;
    const auto &t = this->triangles[tri];

    const glm::vec2 p0 = to_xy(this->points[t.a]);
    const glm::vec2 p1 = to_xy(this->points[t.b]);
    const glm::vec2 p2 = to_xy(this->points[t.c]);

    glm::vec2 diff0 = p - p0;
    glm::vec2 diff1 = p - p1;
    glm::vec2 diff2 = p - p2;

    float d0 = glm::dot(diff0, diff0);
    float d1 = glm::dot(diff1, diff1);
    float d2 = glm::dot(diff2, diff2);

    if (d0 <= d1 && d0 <= d2)
      return this->points[t.a].z;
    else if (d1 <= d2)
      return this->points[t.b].z;
    else
      return this->points[t.c].z;
  }

  return fill_value;
}

int TerrainTriMesh::neighbor_triangle(int tri_index, int edge_index) const
{
  int opp = this->halfedges[3 * tri_index + edge_index];
  return (opp < 0) ? -1 : (opp / 3);
}

void TerrainTriMesh::print_info() const
{
  std::cout << this->info_string();
}

void TerrainTriMesh::remap_z(float vmin, float vmax)
{
  if (this->points.empty()) return;

  // find current z range
  float zmin = std::numeric_limits<float>::max();
  float zmax = std::numeric_limits<float>::lowest();

  for (const auto &p : this->points)
  {
    zmin = std::min(zmin, p.z);
    zmax = std::max(zmax, p.z);
  }

  float dz = zmax - zmin;

  // avoid division by zero (flat terrain)
  if (dz < 1e-30f)
  {
    float zmid = 0.5f * (vmin + vmax);
    for (auto &p : this->points)
      p.z = zmid;
    return;
  }

  float scale = (vmax - vmin) / dz;

  // remap
  for (auto &p : this->points)
    p.z = vmin + (p.z - zmin) * scale;
}

void TerrainTriMesh::relax_xy(float lambda, int iterations, bool preserve_chull)
{
  if (this->neighbors.adjacency.empty()) this->compute_neighbors();

  for (int it = 0; it < iterations; ++it)
  {
    std::vector<glm::vec3> new_points = this->points;

    for (size_t i = 0; i < this->points.size(); ++i)
    {
      const auto &nbrs = this->neighbors.adjacency[i];
      if (nbrs.empty()) continue;

      glm::vec2 avg(0.f);

      for (const auto &nb : nbrs)
        avg += this->to_xy(this->points[nb.index]);

      avg /= static_cast<float>(nbrs.size());

      glm::vec2 current = this->to_xy(this->points[i]);
      glm::vec2 relaxed = glm::mix(current, avg, lambda);

      new_points[i].x = relaxed.x;
      new_points[i].y = relaxed.y;
    }

    // preserve convex hull
    if (preserve_chull)
    {
      for (const auto &i : this->convex_hull)
      {
        new_points[i].x = points[i].x;
        new_points[i].y = points[i].y;
      }
    }

    this->points = std::move(new_points);
  }
}

void TerrainTriMesh::relax_xyz(float lambda,
                               int   iterations,
                               bool  preserve_chull)
{
  if (this->neighbors.adjacency.empty()) this->compute_neighbors();

  for (int it = 0; it < iterations; ++it)
  {
    std::vector<glm::vec3> buffer(this->points.size());

    for (size_t i = 0; i < this->points.size(); ++i)
    {
      const auto &nbrs = this->neighbors.adjacency[i];

      if (nbrs.empty())
      {
        buffer[i] = this->points[i];
        continue;
      }

      glm::vec3 sum(0.f);
      float     sum_w = 0.f;

      const glm::vec3 &pi = this->points[i];

      for (const auto &nb : nbrs)
      {
        size_t           j = nb.index;
        const glm::vec3 &pj = this->points[j];

        float dist = glm::length(pi - pj);

        if (dist > 1e-30f)
        {
          float w = 1.f / dist;
          sum += w * pj;
          sum_w += w;
        }
      }

      if (sum_w > 0.f)
      {
        glm::vec3 avg = sum / sum_w;
        buffer[i] = glm::mix(pi, avg, lambda);
      }
      else
      {
        buffer[i] = pi;
      }
    }

    // preserve convex hull
    if (preserve_chull)
    {
      for (size_t i : this->get_convex_hull())
        buffer[i] = this->points[i];
    }

    this->points = std::move(buffer);
  }
}

void TerrainTriMesh::relax_xyz_taubin(float lambda,
                                      float mu,
                                      int   iterations,
                                      bool  preserve_chull)
{
  for (int it = 0; it < iterations; ++it)
  {
    int sub_iterations = 1;
    this->relax_xyz(lambda, sub_iterations, preserve_chull);
    this->relax_xyz(mu, sub_iterations, preserve_chull);
  }
}

size_t TerrainTriMesh::size() const
{
  return this->points.size();
}

void TerrainTriMesh::slope_limiter(const std::vector<float> &max_slope,
                                   int                       iterations,
                                   float                     sigma)
{
  if (this->neighbors.adjacency.empty()) this->compute_neighbors();

  for (int it = 0; it < iterations; ++it)
  {
    std::vector<glm::vec3> buffer(this->points.size());

    for (size_t i = 0; i < this->points.size(); ++i)
    {
      const auto      &nbrs = this->neighbors.adjacency[i];
      const glm::vec3 &pi = this->points[i];

      buffer[i] = pi;
      float slope_limit = max_slope[i];

      if (nbrs.empty()) continue;

      for (const auto &nb : nbrs)
      {
        size_t           j = nb.index;
        const glm::vec3 &pj = this->points[j];

        float dist = glm::length(glm::vec2(pi.x - pj.x, pi.y - pj.y));
        float dz = pi.z - pj.z;
        float slope = std::abs(dz) / dist;

        if (slope > slope_limit)
        {
          float diff = slope_limit * dist;
          buffer[i].z -= sigma * std::copysignf(1.f, dz) * diff;
        }
      }
    }

    this->points = std::move(buffer);
  }
}

void TerrainTriMesh::slope_limiter(float max_slope, int iterations, float sigma)
{
  std::vector<float> slope_vec(this->size(), max_slope);
  this->slope_limiter(slope_vec, iterations, sigma);
}

void TerrainTriMesh::subdivise()
{
  std::unordered_map<Edge, size_t, EdgeHash> midpoint_map;

  // Reserve a bit to avoid reallocations (optional but nice)
  this->points.reserve(this->points.size() * 2);

  // Lambda to create or reuse midpoint
  auto get_midpoint = [&](size_t i0, size_t i1) -> size_t
  {
    Edge e(i0, i1);
    auto it = midpoint_map.find(e);
    if (it != midpoint_map.end()) return it->second;

    glm::vec3 mid = 0.5f * (points[i0] + points[i1]);
    size_t    idx = points.size();
    points.push_back(mid);

    midpoint_map[e] = idx;
    return idx;
  };

  // Only generate points — no triangles
  for (const auto &tri : triangles)
  {
    get_midpoint(tri.a, tri.b);
    get_midpoint(tri.b, tri.c);
    get_midpoint(tri.c, tri.a);
  }

  // discard old topology
  this->triangles.clear();

  // rebuild triangulation from scratch and recompute adjacency
  this->triangulate_delaunay();
  this->compute_neighbors();
}

Array TerrainTriMesh::to_array(const glm::ivec2         &shape,
                               const std::vector<float> &values,
                               const glm::vec4          &bbox) const
{
  if (!validate_shape(shape)) return Array();

  // prepare data
  std::vector<float> xc, yc, zc;
  const size_t       np = this->points.size();

  if (values.size() != np)
  {
    xc.reserve(np);
    yc.reserve(np);
    zc.reserve(np);

    for (const auto &p : this->points)
    {
      xc.push_back(p.x);
      yc.push_back(p.y);
      zc.push_back(p.z);
    }
  }
  else
  {
    xc.resize(np);
    yc.resize(np);
    zc.resize(np);

    for (size_t k = 0; k < np; ++k)
    {
      xc[k] = this->points[k].x;
      yc[k] = this->points[k].y;
      zc[k] = values[k];
    }
  }

  // interpolate
  Array out = hmap::interpolate2d(shape,
                                  xc,
                                  yc,
                                  zc,
                                  hmap::InterpolationMethod2D::ITP2D_DELAUNAY,
                                  nullptr,
                                  nullptr,
                                  bbox);
  return out;
}

void TerrainTriMesh::to_csv(const std::string &fname) const
{
  std::ofstream f(fname, std::ios::out);
  if (!f.is_open()) return;

  f << "# vertices\n";
  f << "# "
       "vertex_id,x,y,z\n";

  for (size_t i = 0; i < static_cast<size_t>(this->get_points().size()); ++i)
  {
    f << i << "," << this->get_points()[i].x << "," << this->get_points()[i].y
      << "," << this->get_points()[i].z << "\n";
  }

  f.close();
}

glm::vec2 TerrainTriMesh::to_xy(const glm::vec3 &p) const
{
  return {p.x, p.y};
}

void TerrainTriMesh::triangulate_delaunay()
{
  // generate Delaunay triangulation
  std::vector<float> xy;
  xy.reserve(2 * this->points.size());

  for (const auto &p : this->points)
  {
    xy.push_back(p.x);
    xy.push_back(p.y);
  }

  delaunator::Delaunator d(xy);

  // triangles
  {
    const std::vector<size_t> &tri = d.triangles;

    this->triangles.clear();
    this->triangles.reserve(tri.size() / 3);

    this->halfedges.clear();
    this->halfedges.reserve(tri.size());

    for (std::size_t k = 0; k < tri.size(); k += 3)
    {
      this->triangles.push_back({tri[k], tri[k + 1], tri[k + 2]});

      this->halfedges.push_back(d.halfedges[k]);
      this->halfedges.push_back(d.halfedges[k + 1]);
      this->halfedges.push_back(d.halfedges[k + 2]);
    }
  }

  // chull
  {
    this->convex_hull.clear();
    size_t start = d.hull_start;
    size_t i = start;

    do
    {
      this->convex_hull.push_back(i);
      i = d.hull_next[i];
    } while (i != start);
  }
}

// ------------------------------------------------------------
// FUNCTIONS
// ------------------------------------------------------------

std::vector<float> cubic_pulse(const TerrainTriMesh &mesh)
{
  const TerrainTriMesh::BoundingBox &bbox = mesh.get_bbox();

  // center and radius
  float xc = 0.5f * (bbox.max.x + bbox.min.x);
  float yc = 0.5f * (bbox.max.y + bbox.min.y);

  float lx = bbox.max.x - bbox.min.x;
  float ly = bbox.max.y - bbox.min.y;
  float rc = std::min(lx, ly);

  std::vector<float> values;
  values.reserve(mesh.size());

  for (auto &p : mesh.get_points())
  {
    float r = std::hypot(p.x - xc, p.y - yc) / rc;

    if (r < 1.f)
      r = 1.f - r * r * (3.f - 2.f * r);
    else
      r = 0.f;

    values.push_back(r);
  }

  return values;
}

TerrainTriMesh generate_terrain_tri_mesh_from_heightmap(const Array &z,
                                                        float        max_error,
                                                        int max_triangles,
                                                        int max_points)
{
  if (!validate_non_empty(z)) return TerrainTriMesh();

  const glm::ivec2 &shape = z.shape;

  const auto   p_hmap = std::make_shared<Heightmap>(shape.y,
                                                  shape.x,
                                                  z.get_vector());
  Triangulator tri(p_hmap);
  tri.Run(max_error, max_triangles, max_points);

  const auto &points = tri.Points(1.f);

  // x, y normalization coefficients
  const float ax = 1.f / float(shape.y - 1);
  const float ay = 1.f / float(shape.x - 1);

  std::vector<glm::vec3> xyz;
  xyz.reserve(points.size());

  for (const auto &p : points)
    xyz.push_back({ay * p.y, ax * p.x, p.z});

  return TerrainTriMesh(xyz);
}

TerrainTriMesh generate_terrain_tri_mesh_from_heightmap_random(
    const Array  &z,
    int           control_points_count,
    std::uint32_t seed)
{
  if (!validate_non_empty(z)) return TerrainTriMesh();

  const glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f};

  Cloud cloud = random_cloud_jittered(control_points_count,
                                      {0.5f, 0.5f},
                                      {0.f, 0.f},
                                      seed,
                                      bbox);
  cloud.snap_points_to_bounding_box(bbox);
  cloud.set_values_from_array(z, bbox);

  return TerrainTriMesh(cloud.to_vec3());
}

} // namespace hmap
