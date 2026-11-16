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

#include <functional>

#include "highmap/array.hpp"

#include "NoiseLib/include/controlfunction.h"

namespace hmap
{

class ArrayControlFunction : public ControlFunction<ArrayControlFunction>
{
  friend class ControlFunction<ArrayControlFunction>;

public:
  explicit ArrayControlFunction(hmap::Array array);

protected:
  float EvaluateImpl(float x, float y) const;
  bool  InsideDomainImpl(float x, float y) const;
  float DistToDomainImpl(float x, float y) const;
  float MinimumImpl() const;
  float MaximumImpl() const;

private:
  float get(int i, int j) const;
  float sample(float ri, float rj) const;

  const hmap::Array m_array;
};

class XyControlFunction : public ControlFunction<XyControlFunction>
{
  friend class ControlFunction<XyControlFunction>;

public:
  XyControlFunction(std::function<float(float, float, float)> noise_function,
                    float                                     offset = 0.f,
                    float                                     scaling = 1.f);

protected:
  double EvaluateImpl(float x, float y) const;
  bool   InsideDomainImpl(float /* x */, float /* y */) const;
  double DistToDomainImpl(float /* x */, float /* y */) const;
  double MinimumImpl() const;
  double MaximumImpl() const;

private:
  std::function<float(float, float, float)> noise_function;
  float                                     offset;
  float                                     scaling;
};

} // namespace hmap
