/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <vector>

#include "delaunator-cpp.hpp"

#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/geometry/point_sampling.hpp"
#include "highmap/math/core.hpp"
#include "highmap/primitives/random.hpp"
#include "highmap/random.hpp"
#include "highmap/range.hpp"

namespace hmap
{

namespace
{

uint32_t hash21u_seed(int ix, int iy, uint32_t seed)
{
  uint32_t h = static_cast<uint32_t>(ix) * 0x8da6b343u ^
               static_cast<uint32_t>(iy) * 0xd8163841u ^ seed * 0x9e3779b9u;
  return hash32(h);
}

float hash12u_f(int ix, int iy, uint32_t seed)
{
  return static_cast<float>(hash21u_seed(ix, iy, seed) >> 8) *
         (1.f / 16777216.f);
}

glm::vec2 hash22u_f(int ix, int iy, uint32_t seed)
{
  uint32_t h0 = hash21u_seed(ix, iy, seed);
  uint32_t h1 = hash32(h0 ^ 0x9e3779b9u);
  return {static_cast<float>(h0 >> 8) * (1.f / 16777216.f),
          static_cast<float>(h1 >> 8) * (1.f / 16777216.f)};
}

glm::vec2 grad22u_f(int ix, int iy, uint32_t seed)
{
  float angle = 6.2831853f * hash12u_f(ix, iy, seed);
  return {std::cos(angle), std::sin(angle)};
}

float smoothstep3_f(float t)
{
  return t * t * (3.f - 2.f * t);
}

float smoothstep5_f(float t)
{
  return t * t * t * (t * (t * 6.f - 15.f) + 10.f);
}

float cubic_interp_f(float p0, float p1, float p2, float p3, float t)
{
  float a = p3 - p2 - p0 + p1;
  float b = p0 - p1 - a;
  float c = p2 - p0;
  float d = p1;
  return a * t * t * t + b * t * t + c * t + d;
}

float base_perlin(float px, float py, uint32_t seed)
{
  float pfx = std::floor(px);
  float pfy = std::floor(py);
  int   ix = static_cast<int>(pfx);
  int   iy = static_cast<int>(pfy);
  float fx = px - pfx;
  float fy = py - pfy;

  glm::vec2 g00 = grad22u_f(ix, iy, seed);
  glm::vec2 g10 = grad22u_f(ix + 1, iy, seed);
  glm::vec2 g01 = grad22u_f(ix, iy + 1, seed);
  glm::vec2 g11 = grad22u_f(ix + 1, iy + 1, seed);

  float n00 = g00.x * fx + g00.y * fy;
  float n10 = g10.x * (fx - 1.f) + g10.y * fy;
  float n01 = g01.x * fx + g01.y * (fy - 1.f);
  float n11 = g11.x * (fx - 1.f) + g11.y * (fy - 1.f);

  float ux = smoothstep5_f(fx);
  float uy = smoothstep5_f(fy);

  float nx0 = n00 + ux * (n10 - n00);
  float nx1 = n01 + ux * (n11 - n01);
  float nxy = nx0 + uy * (nx1 - nx0);

  return 1.42857f * nxy;
}

float base_simplex2(float px, float py, uint32_t seed)
{
  const float K1 = 0.366025403f;
  const float K2 = 0.211324865f;

  float p2x = 0.5f * px;
  float p2y = 0.5f * py;

  float s = (p2x + p2y) * K1;
  float ipx = std::floor(p2x + s);
  float ipy = std::floor(p2y + s);
  int   ix = static_cast<int>(ipx);
  int   iy = static_cast<int>(ipy);

  float t = (ipx + ipy) * K2;
  float ax = p2x - ipx + t;
  float ay = p2y - ipy + t;

  int   ox = (ax >= ay) ? 1 : 0;
  int   oy = 1 - ox;
  float bx = ax - static_cast<float>(ox) + K2;
  float by = ay - static_cast<float>(oy) + K2;
  float cx = ax - 1.f + 2.f * K2;
  float cy = ay - 1.f + 2.f * K2;

  glm::vec2 g0 = grad22u_f(ix, iy, seed);
  glm::vec2 g1 = grad22u_f(ix + ox, iy + oy, seed);
  glm::vec2 g2 = grad22u_f(ix + 1, iy + 1, seed);

  float n0 = std::max(0.5f - (ax * ax + ay * ay), 0.f);
  float n1 = std::max(0.5f - (bx * bx + by * by), 0.f);
  float n2 = std::max(0.5f - (cx * cx + cy * cy), 0.f);

  n0 = n0 * n0 * n0 * n0 * (g0.x * ax + g0.y * ay);
  n1 = n1 * n1 * n1 * n1 * (g1.x * bx + g1.y * by);
  n2 = n2 * n2 * n2 * n2 * (g2.x * cx + g2.y * cy);

  return 100.f * (n0 + n1 + n2);
}

float base_value(float px, float py, uint32_t seed)
{
  float pfx = std::floor(px);
  float pfy = std::floor(py);
  int   ix = static_cast<int>(pfx);
  int   iy = static_cast<int>(pfy);
  float fx = px - pfx;
  float fy = py - pfy;

  float v00 = hash12u_f(ix, iy, seed);
  float v10 = hash12u_f(ix + 1, iy, seed);
  float v01 = hash12u_f(ix, iy + 1, seed);
  float v11 = hash12u_f(ix + 1, iy + 1, seed);

  float ux = smoothstep3_f(fx);
  float uy = smoothstep3_f(fy);

  float nx0 = v00 + ux * (v10 - v00);
  float nx1 = v01 + ux * (v11 - v01);
  float nxy = nx0 + uy * (nx1 - nx0);

  return 2.f * nxy - 1.f;
}

float base_value_cubic(float px, float py, uint32_t seed)
{
  float pfx = std::floor(px);
  float pfy = std::floor(py);
  int   ix = static_cast<int>(pfx);
  int   iy = static_cast<int>(pfy);
  float fx = px - pfx;
  float fy = py - pfy;

  float v[4][4];
  for (int dy = -1; dy <= 2; dy++)
    for (int dx = -1; dx <= 2; dx++)
      v[dy + 1][dx + 1] = hash12u_f(ix + dx, iy + dy, seed);

  float interp_row[4];
  for (int dy = 0; dy < 4; dy++)
    interp_row[dy] = cubic_interp_f(v[dy][0], v[dy][1], v[dy][2], v[dy][3], fx);

  float value = cubic_interp_f(interp_row[0],
                               interp_row[1],
                               interp_row[2],
                               interp_row[3],
                               fy);

  return 1.43f * (value - 0.5f);
}

float base_value_linear(float px, float py, uint32_t seed)
{
  float pfx = std::floor(px);
  float pfy = std::floor(py);
  int   ix = static_cast<int>(pfx);
  int   iy = static_cast<int>(pfy);
  float fx = px - pfx;
  float fy = py - pfy;

  float v00 = hash12u_f(ix, iy, seed);
  float v10 = hash12u_f(ix + 1, iy, seed);
  float v01 = hash12u_f(ix, iy + 1, seed);
  float v11 = hash12u_f(ix + 1, iy + 1, seed);

  float nx0 = v00 + fx * (v10 - v00);
  float nx1 = v01 + fx * (v11 - v01);
  float nxy = nx0 + fy * (nx1 - nx0);

  return 2.f * nxy - 1.f;
}

float base_worley(float px, float py, uint32_t seed)
{
  float pfx = std::floor(px);
  float pfy = std::floor(py);
  int   ix = static_cast<int>(pfx);
  int   iy = static_cast<int>(pfy);
  float fx = px - pfx;
  float fy = py - pfy;

  float min_dist = std::numeric_limits<float>::max();

  for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    {
      glm::vec2 fp = glm::vec2(static_cast<float>(dx), static_cast<float>(dy)) +
                     hash22u_f(ix + dx, iy + dy, seed);
      float diffx = fx - fp.x;
      float diffy = fy - fp.y;
      float dist = diffx * diffx + diffy * diffy;

      min_dist = std::min(min_dist, dist);
    }

  return 1.66f * min_dist - 1.f;
}

float base_worley_value(float px, float py, uint32_t seed)
{
  float pfx = std::floor(px);
  float pfy = std::floor(py);
  int   ix = static_cast<int>(pfx);
  int   iy = static_cast<int>(pfy);
  float fx = px - pfx;
  float fy = py - pfy;

  float min_dist = std::numeric_limits<float>::max();
  int   best_x = ix;
  int   best_y = iy;

  for (int dx = -1; dx <= 1; dx++)
    for (int dy = -1; dy <= 1; dy++)
    {
      glm::vec2 fp = glm::vec2(static_cast<float>(dx), static_cast<float>(dy)) +
                     hash22u_f(ix + dx, iy + dy, seed);
      float diffx = fx - fp.x;
      float diffy = fy - fp.y;
      float dist = diffx * diffx + diffy * diffy;

      if (dist < min_dist)
      {
        min_dist = dist;
        best_x = ix + dx;
        best_y = iy + dy;
      }
    }

  return 2.f * hash12u_f(best_x, best_y, seed ^ 0x517cc1b7u) - 1.f;
}

} // namespace

//----------------------------------------------------------------------
// derived from NoiseFunction class
//----------------------------------------------------------------------

PerlinFunction::PerlinFunction(glm::vec2 kw, std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
}

float PerlinFunction::operator()(float x, float y, float) const
{
  return base_perlin(this->kw.x * x, this->kw.y * y, this->seed);
}

PerlinBillowFunction::PerlinBillowFunction(glm::vec2 kw, std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
}

float PerlinBillowFunction::operator()(float x, float y, float) const
{
  float value = base_perlin(this->kw.x * x, this->kw.y * y, this->seed);
  return 2.f * std::abs(value) - 1.f;
}

PerlinHalfFunction::PerlinHalfFunction(glm::vec2     kw,
                                       std::uint32_t seed,
                                       float         k)
    : NoiseFunction(kw, seed), k(k)
{
}

float PerlinHalfFunction::operator()(float x, float y, float) const
{
  float value = base_perlin(this->kw.x * x, this->kw.y * y, this->seed);
  return clamp_min_smooth(value, 0.f, this->k);
}

PerlinMixFunction::PerlinMixFunction(glm::vec2 kw, std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
}

float PerlinMixFunction::operator()(float x, float y, float) const
{
  float value = base_perlin(this->kw.x * x, this->kw.y * y, this->seed);
  return 0.5f * value + std::abs(value) - 0.5f;
}

Simplex2Function::Simplex2Function(glm::vec2 kw, std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
}

float Simplex2Function::operator()(float x, float y, float) const
{
  return base_simplex2(this->kw.x * x, this->kw.y * y, this->seed);
}

Simplex2SFunction::Simplex2SFunction(glm::vec2 kw, std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
}

float Simplex2SFunction::operator()(float x, float y, float) const
{
  return base_simplex2(this->kw.x * x, this->kw.y * y, this->seed);
}

ValueNoiseFunction::ValueNoiseFunction(glm::vec2 kw, std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
}

float ValueNoiseFunction::operator()(float x, float y, float) const
{
  return base_value(this->kw.x * x, this->kw.y * y, this->seed);
}

ValueCubicNoiseFunction::ValueCubicNoiseFunction(glm::vec2     kw,
                                                 std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
}

float ValueCubicNoiseFunction::operator()(float x, float y, float) const
{
  return base_value_cubic(this->kw.x * x, this->kw.y * y, this->seed);
}

ValueDelaunayNoiseFunction::ValueDelaunayNoiseFunction(glm::vec2     kw,
                                                       std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
  this->set_kw(kw);
  this->set_seed(seed);
  this->update_interpolation_function();
}

void ValueDelaunayNoiseFunction::update_interpolation_function()
{
  // --- generate 'n' random grid points
  int n = (int)(this->kw.x * this->kw.y);

  std::vector<float> x(n);
  std::vector<float> y(n);
  std::vector<float> value(n);

  auto xy = random_points(n,
                          this->seed,
                          PointSamplingMethod::RND_LHS,
                          BBox2D{});
  x = xy[0];
  y = xy[1];

  expand_points_domain(x, y, value, BBox2D{});

  // --- interpolation function

  // triangulate
  std::vector<float> coords(2 * x.size());

  for (size_t k = 0; k < x.size(); k++)
  {
    coords[2 * k] = x[k];
    coords[2 * k + 1] = y[k];
  }

  delaunator::Delaunator<float> d(coords);

  // compute and store triangles area
  std::vector<float> area(d.triangles.size());

  for (size_t k = 0; k < d.triangles.size(); k += 3)
  {
    int p0 = d.triangles[k];
    int p1 = d.triangles[k + 1];
    int p2 = d.triangles[k + 2];

    // true area
    area[k] = 0.5f * (-y[p1] * x[p2] + y[p0] * (-x[p1] + x[p2]) +
                      x[p0] * (y[p1] - y[p2]) + x[p1] * y[p2]);

    // but stored like this to avoid doing it at each evaluation while
    // interpolating
    area[k] = 1.f / (2.f * area[k]);
  }

  auto itp_fct = [x, y, value, d, area](float x_, float y_, float)
  {
    // https://stackoverflow.com/questions/2049582

    // compute barycentric coordinates to find in which triangle the
    // point (x_, y_) is inside
    for (size_t k = 0; k < d.triangles.size(); k += 3)
    {
      int p0 = d.triangles[k];
      int p1 = d.triangles[k + 1];
      int p2 = d.triangles[k + 2];

      float s = area[k] * (y[p0] * x[p2] - x[p0] * y[p2] +
                           (y[p2] - y[p0]) * x_ + (x[p0] - x[p2]) * y_);
      float t = area[k] * (x[p0] * y[p1] - y[p0] * x[p1] +
                           (y[p0] - y[p1]) * x_ + (x[p1] - x[p0]) * y_);

      if (s >= 0.f && t >= 0.f && s + t <= 1.f)
        return value[p0] + s * (value[p1] - value[p0]) +
               t * (value[p2] - value[p0]);
    }

    return 1.f;
  };

  this->itp_fct = itp_fct;
}

float ValueDelaunayNoiseFunction::operator()(float x,
                                             float y,
                                             float ctrl_param) const
{
  return this->itp_fct(x, y, ctrl_param);
}

ValueLinearNoiseFunction::ValueLinearNoiseFunction(glm::vec2     kw,
                                                   std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
  this->set_kw(kw);
  this->set_seed(seed);
  this->update_interpolation_function();
}

void ValueLinearNoiseFunction::update_interpolation_function()
{
  glm::vec4 bbox = {-1.f, 2.f, -1.f, 2.f}; // bounding box

  float lx = bbox.y - bbox.x;
  float ly = bbox.w - bbox.z;

  glm::ivec2 shape_base = glm::ivec2((int)(this->kw.x * lx) + 1,
                                     (int)(this->kw.y * ly) + 1);

  Array value = 2.f * white(shape_base, 0.f, 1.f, seed) - 1.f;

  // corresponding grids
  std::vector<float> xv(shape_base.x);
  std::vector<float> yv(shape_base.y);

  for (int i = 0; i < shape_base.x; i++)
    xv[i] = bbox.x + lx * (float)i / (float)(shape_base.x - 1);

  for (int j = 0; j < shape_base.y; j++)
    yv[j] = bbox.z + ly * (float)j / (float)(shape_base.y - 1);

  auto itp_fct = [xv, yv, value, bbox](float x_, float y_, float)
  {
    float xn = (x_ - bbox.x) / (bbox.y - bbox.x) * (float)(value.shape.x - 1);
    float yn = (y_ - bbox.z) / (bbox.w - bbox.z) * (float)(value.shape.y - 1);

    int in = (int)xn;
    int jn = (int)yn;

    float u = xn - (float)in;
    float v = yn - (float)jn;

    if (in == value.shape.x - 1)
    {
      in = value.shape.x - 2;
      u = 1.f;
    }

    if (jn == value.shape.y - 1)
    {
      jn = value.shape.y - 2;
      v = 1.f;
    }

    return value.get_value_bilinear_at(in, jn, u, v);
  };

  this->itp_fct = itp_fct;
}

float ValueLinearNoiseFunction::operator()(float x,
                                           float y,
                                           float /*ctrl_param*/) const
{
  return base_value_linear(this->kw.x * x, this->kw.y * y, this->seed);
}

WorleyDoubleFunction::WorleyDoubleFunction(glm::vec2     kw,
                                           std::uint32_t seed,
                                           float         ratio,
                                           float         k)
    : NoiseFunction(kw, seed), ratio(ratio), k(k)
{
}

float WorleyDoubleFunction::operator()(float x, float y, float ctrl_param) const
{
  float local_ratio = ctrl_param * this->ratio;

  float w1 = base_worley(this->kw.x * x, this->kw.y * y, this->seed);
  float w2 = base_worley(this->kw.x * x * 2.f,
                         this->kw.y * y * 2.f,
                         this->seed + 1u);
  if (this->k)
    return maximum_smooth(local_ratio * w1, (1.f - local_ratio) * w2, this->k);
  else
    return std::max(local_ratio * w1, (1.f - local_ratio) * w2);
}

WorleyFunction::WorleyFunction(glm::vec2     kw,
                               std::uint32_t seed,
                               bool          return_cell_value)
    : NoiseFunction(kw, seed), return_cell_value(return_cell_value)
{
}

float WorleyFunction::operator()(float x, float y, float) const
{
  if (this->return_cell_value)
    return base_worley_value(this->kw.x * x, this->kw.y * y, this->seed);
  else
    return base_worley(this->kw.x * x, this->kw.y * y, this->seed);
}

// --- helper

} // namespace hmap
