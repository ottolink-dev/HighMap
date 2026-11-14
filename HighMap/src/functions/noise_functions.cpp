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

#include "FastNoiseLite.h"
#include "delaunator-cpp.hpp"

#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/geometry/point_sampling.hpp"
#include "highmap/primitives/random.hpp"
#include "highmap/range.hpp"

namespace hmap
{

//----------------------------------------------------------------------
// derived from NoiseFunction class
//----------------------------------------------------------------------

PerlinFunction::PerlinFunction(glm::vec2 kw, std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
  this->set_seed(seed);
  this->noise.SetFrequency(1.f);
  this->noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
}

float PerlinFunction::operator()(float x, float y, float) const
{
  return this->noise.GetNoise(this->kw.x * x, this->kw.y * y);
}

PerlinBillowFunction::PerlinBillowFunction(glm::vec2 kw, std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
  this->set_seed(seed);
  this->noise.SetFrequency(1.f);
  this->noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
}

float PerlinBillowFunction::operator()(float x, float y, float) const
{
  float value = this->noise.GetNoise(this->kw.x * x, this->kw.y * y);
  return 2.f * std::abs(value) - 1.f;
}

PerlinHalfFunction::PerlinHalfFunction(glm::vec2     kw,
                                       std::uint32_t seed,
                                       float         k)
    : NoiseFunction(kw, seed), k(k)
{
  this->set_seed(seed);
  this->noise.SetFrequency(1.f);
  this->noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
}

float PerlinHalfFunction::operator()(float x, float y, float) const
{
  float value = this->noise.GetNoise(this->kw.x * x, this->kw.y * y);
  return clamp_min_smooth(value, 0.f, this->k);
}

PerlinMixFunction::PerlinMixFunction(glm::vec2 kw, std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
  this->set_seed(seed);
  this->noise.SetFrequency(1.f);
  this->noise.SetNoiseType(FastNoiseLite::NoiseType_Perlin);
}

float PerlinMixFunction::operator()(float x, float y, float) const
{
  float value = this->noise.GetNoise(this->kw.x * x, this->kw.y * y);
  return 0.5f * value + std::abs(value) - 0.5f;
}

Simplex2Function::Simplex2Function(glm::vec2 kw, std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
  this->set_seed(seed);
  this->noise.SetFrequency(0.5f);
  this->noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2);
}

float Simplex2Function::operator()(float x, float y, float) const
{
  return this->noise.GetNoise(this->kw.x * x, this->kw.y * y);
}

Simplex2SFunction::Simplex2SFunction(glm::vec2 kw, std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
  this->set_seed(seed);
  this->noise.SetFrequency(0.5f);
  this->noise.SetNoiseType(FastNoiseLite::NoiseType_OpenSimplex2S);
}

float Simplex2SFunction::operator()(float x, float y, float) const
{
  return this->noise.GetNoise(this->kw.x * x, this->kw.y * y);
}

ValueNoiseFunction::ValueNoiseFunction(glm::vec2 kw, std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
  this->set_seed(seed);
  this->noise.SetFrequency(1.f);
  this->noise.SetNoiseType(FastNoiseLite::NoiseType_Value);
}

float ValueNoiseFunction::operator()(float x, float y, float) const
{
  return this->noise.GetNoise(this->kw.x * x, this->kw.y * y);
}

ValueCubicNoiseFunction::ValueCubicNoiseFunction(glm::vec2     kw,
                                                 std::uint32_t seed)
    : NoiseFunction(kw, seed)
{
  this->set_seed(seed);
  this->noise.SetFrequency(1.f);
  this->noise.SetNoiseType(FastNoiseLite::NoiseType_ValueCubic);
}

float ValueCubicNoiseFunction::operator()(float x, float y, float) const
{
  return 1.43f * this->noise.GetNoise(this->kw.x * x, this->kw.y * y);
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
                          {0.f, 1.f, 0.f, 1.f});
  x = xy[0];
  y = xy[1];

  expand_points_domain(x, y, value, {0.f, 1.f, 0.f, 1.f});

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
                                           float ctrl_param) const
{
  return this->itp_fct(x, y, ctrl_param);
}

WorleyDoubleFunction::WorleyDoubleFunction(glm::vec2     kw,
                                           std::uint32_t seed,
                                           float         ratio,
                                           float         k)
    : NoiseFunction(kw, seed), ratio(ratio), k(k)
{
  this->set_seed(seed);
  this->noise1.SetFrequency(1.0f);
  this->noise2.SetFrequency(1.0f);

  this->noise1.SetNoiseType(FastNoiseLite::NoiseType_Cellular);
  this->noise2.SetNoiseType(FastNoiseLite::NoiseType_Cellular);
}

float WorleyDoubleFunction::operator()(float x, float y, float ctrl_param) const
{
  float local_ratio = ctrl_param * this->ratio;

  float w1 = this->noise1.GetNoise(this->kw.x * x, this->kw.y * y);
  float w2 = this->noise2.GetNoise(this->kw.x * x, this->kw.y * y);
  if (this->k)
    return maximum_smooth(local_ratio * w1, (1.f - local_ratio) * w2, this->k);
  else
    return std::max(local_ratio * w1, (1.f - local_ratio) * w2);
}

WorleyFunction::WorleyFunction(glm::vec2     kw,
                               std::uint32_t seed,
                               bool          return_cell_value)
    : NoiseFunction(kw, seed)
{
  this->set_seed(seed);
  this->noise.SetFrequency(1.f);
  this->noise.SetNoiseType(FastNoiseLite::NoiseType_Cellular);
  this->noise.SetCellularJitter(1.f);

  if (return_cell_value)
    this->noise.SetCellularReturnType(
        FastNoiseLite::CellularReturnType_CellValue);
  else
    this->noise.SetCellularReturnType(
        FastNoiseLite::CellularReturnType_Distance);
}

float WorleyFunction::operator()(float x, float y, float) const
{
  return 1.66f * (0.4f + this->noise.GetNoise(this->kw.x * x, this->kw.y * y));
}

// --- helper

} // namespace hmap
