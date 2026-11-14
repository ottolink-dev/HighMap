/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file dendry_array_control_function.hpp
 * @author  Otto Link (otto.link.bv@gmail.com)
 * @brief
 * @version 0.1
 * @date 2024-01-01
 *
 * @copyright Copyright (c) 2024
 *
 */
#pragma once

#include <utility>

#include "highmap/array.hpp"
#include "highmap/functions.hpp"

#include "NoiseLib/include/controlfunction.h"
#include "NoiseLib/include/math2d.h"
#include "NoiseLib/include/utils.h"

#include <cassert>

namespace hmap
{

class ArrayControlFunction : public ControlFunction<ArrayControlFunction>
{
  friend class ControlFunction<ArrayControlFunction>;

public:
  explicit ArrayControlFunction(hmap::Array array) : m_array(std::move(array))
  {
  }

protected:
  float EvaluateImpl(float x, float y) const
  {
    x = std::clamp(x, 0.f, 1.f);
    y = std::clamp(y, 0.f, 1.f);

    return sample(x, y);
  }

  bool InsideDomainImpl(float x, float y) const
  {
    return x >= 0.f && x <= 1.f && y >= 0.f && y <= 1.f;
  }

  float DistToDomainImpl(float x, float y) const
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

  float MinimumImpl() const
  {
    return 0.f;
  }

  float MaximumImpl() const
  {
    return 1.f;
  }

private:
  float get(int i, int j) const
  {
    return (float)m_array(i, j);
  }

  float sample(float ri, float rj) const
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

  const hmap::Array m_array;
};

class XyControlFunction : public ControlFunction<XyControlFunction>
{
  friend class ControlFunction<XyControlFunction>;

public:
  XyControlFunction(std::function<float(float, float, float)> noise_function,
                    float                                     offset = 0.f,
                    float                                     scaling = 1.f)
      : noise_function(std::move(noise_function)),
        offset(offset),
        scaling(scaling)
  {
  }

protected:
  double EvaluateImpl(float x, float y) const
  {
    return offset + scaling * noise_function(x, y, 0.f);
  }

  bool InsideDomainImpl(float /* x */, float /* y */) const
  {
    return true;
  }

  double DistToDomainImpl(float /* x */, float /* y */) const
  {
    return 0.f;
  }

  double MinimumImpl() const
  {
    return 0.f;
  }

  double MaximumImpl() const
  {
    return 1.f;
  }

private:
  std::function<float(float, float, float)> noise_function;
  float                                     offset;
  float                                     scaling;
};

} // namespace hmap
