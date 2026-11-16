/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#include <algorithm>
#include <limits>
#include <utility>

#include "highmap/primitives/dendry_array_control_function.hpp"

#include "NoiseLib/include/math2d.h"
#include "NoiseLib/include/utils.h"

namespace hmap
{

ArrayControlFunction::ArrayControlFunction(hmap::Array array)
    : m_array(std::move(array))
{
}

float ArrayControlFunction::EvaluateImpl(float x, float y) const
{
  x = std::clamp(x, 0.f, 1.f);
  y = std::clamp(y, 0.f, 1.f);

  return sample(x, y);
}

bool ArrayControlFunction::InsideDomainImpl(float x, float y) const
{
  return x >= 0.f && x <= 1.f && y >= 0.f && y <= 1.f;
}

float ArrayControlFunction::DistToDomainImpl(float x, float y) const
{
  if (InsideDomainImpl(x, y))
  {
    return 0.f;
  }

  const Point2D p(x, y);

  const Point2D topLeft(0.f, 0.f);
  const Point2D topRight(1.f, 0.f);
  const Point2D bottomLeft(0.f, 1.f);
  const Point2D bottomRight(1.f, 1.f);

  Point2D c; // Useless point for distToLineSegment

  auto dist = std::numeric_limits<float>::max();

  dist = std::min(dist, distToLineSegment(p, topLeft, topRight, c));
  dist = std::min(dist, distToLineSegment(p, topRight, bottomRight, c));
  dist = std::min(dist, distToLineSegment(p, bottomRight, bottomLeft, c));
  dist = std::min(dist, distToLineSegment(p, bottomLeft, topLeft, c));

  return dist;
}

float ArrayControlFunction::MinimumImpl() const
{
  return 0.f;
}

float ArrayControlFunction::MaximumImpl() const
{
  return 1.f;
}

float ArrayControlFunction::get(int i, int j) const
{
  return (float)m_array(i, j);
}

float ArrayControlFunction::sample(float ri, float rj) const
{
  float x = ri * (m_array.shape.x - 1);
  float y = rj * (m_array.shape.y - 1);

  int i = (int)x;
  int j = (int)y;

  float u;
  float v;

  if (i == m_array.shape.x - 1)
  {
    i = m_array.shape.x - 2;
    u = 1.f;
  }
  else
    u = x - (float)i;

  if (j == m_array.shape.y - 1)
  {
    j = m_array.shape.y - 2;
    v = 1.f;
  }
  else
    v = y - (float)j;

  return m_array.get_value_bilinear_at(i, j, u, v);
}

XyControlFunction::XyControlFunction(
    std::function<float(float, float, float)> noise_function,
    float                                     offset,
    float                                     scaling)
    : noise_function(std::move(noise_function)),
      offset(offset),
      scaling(scaling)
{
}

double XyControlFunction::EvaluateImpl(float x, float y) const
{
  return offset + scaling * noise_function(x, y, 0.f);
}

bool XyControlFunction::InsideDomainImpl(float /* x */, float /* y */) const
{
  return true;
}

double XyControlFunction::DistToDomainImpl(float /* x */, float /* y */) const
{
  return 0.f;
}

double XyControlFunction::MinimumImpl() const
{
  return 0.f;
}

double XyControlFunction::MaximumImpl() const
{
  return 1.f;
}

} // namespace hmap
