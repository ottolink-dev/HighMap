/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/operator.hpp"

namespace hmap
{

void set_borders(Array &array, glm::vec4 border_values, glm::ivec4 buffer_sizes)
{
  if (!validate_non_empty(array)) return;

  // west
  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < buffer_sizes.x; i++)
    {
      float r = (float)i / (float)buffer_sizes.x;
      r = smoothstep3(r);
      array(i, j) = (1.f - r) * border_values.x + r * array(i, j);
    }

  // east
  for (int j = 0; j < array.shape.y; j++)
    for (int i = array.shape.x - buffer_sizes.y; i < array.shape.x; i++)
    {
      float r = 1.f - (float)(i - array.shape.x + buffer_sizes.y) /
                          (float)buffer_sizes.y;
      r = smoothstep3(r);
      array(i, j) = (1.f - r) * border_values.y + r * array(i, j);
    }

  // south
  for (int j = 0; j < buffer_sizes.z; j++)
    for (int i = 0; i < array.shape.x; i++)
    {
      float r = (float)j / (float)buffer_sizes.z;
      r = smoothstep3(r);
      array(i, j) = (1.f - r) * border_values.z + r * array(i, j);
    }

  // north
  for (int j = array.shape.y - buffer_sizes.w; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
    {
      float r = 1.f - (float)(j - array.shape.y + buffer_sizes.w) /
                          (float)buffer_sizes.w;
      r = smoothstep3(r);
      array(i, j) = (1.f - r) * border_values.w + r * array(i, j);
    }
}

void set_borders(Array    &array,
                 glm::vec4 border_values,
                 glm::vec4 buffer_sizes,
                 glm::vec4 bbox)
{
  if (!validate_non_empty(array)) return;

  std::vector<float> x = linspace(bbox.x, bbox.y, array.shape.x, true);
  std::vector<float> y = linspace(bbox.z, bbox.w, array.shape.y, true);

  // west
  if (buffer_sizes.x > 0.f)
  {
    for (int i = 0; i < array.shape.x; i++)
    {
      float dist = x[i];
      if (dist < buffer_sizes.x)
      {
        float r = smoothstep3(std::clamp(dist / buffer_sizes.x, 0.f, 1.f));
        for (int j = 0; j < array.shape.y; j++)
          array(i, j) = (1.f - r) * border_values.x + r * array(i, j);
      }
    }
  }

  // east
  if (buffer_sizes.y > 0.f)
  {
    for (int i = 0; i < array.shape.x; i++)
    {
      float dist = 1.f - x[i];
      if (dist < buffer_sizes.y)
      {
        float r = smoothstep3(std::clamp(dist / buffer_sizes.y, 0.f, 1.f));
        for (int j = 0; j < array.shape.y; j++)
          array(i, j) = (1.f - r) * border_values.y + r * array(i, j);
      }
    }
  }

  // south
  if (buffer_sizes.z > 0.f)
  {
    for (int j = 0; j < array.shape.y; j++)
    {
      float dist = y[j];
      if (dist < buffer_sizes.z)
      {
        float r = smoothstep3(std::clamp(dist / buffer_sizes.z, 0.f, 1.f));
        for (int i = 0; i < array.shape.x; i++)
          array(i, j) = (1.f - r) * border_values.z + r * array(i, j);
      }
    }
  }

  // north
  if (buffer_sizes.w > 0.f)
  {
    for (int j = 0; j < array.shape.y; j++)
    {
      float dist = 1.f - y[j];
      if (dist < buffer_sizes.w)
      {
        float r = smoothstep3(std::clamp(dist / buffer_sizes.w, 0.f, 1.f));
        for (int i = 0; i < array.shape.x; i++)
          array(i, j) = (1.f - r) * border_values.w + r * array(i, j);
      }
    }
  }
}

void set_borders(Array &array, float border_values, int buffer_sizes)
{
  glm::vec4  bv = glm::vec4(border_values,
                           border_values,
                           border_values,
                           border_values);
  glm::ivec4 bs = glm::ivec4(buffer_sizes,
                             buffer_sizes,
                             buffer_sizes,
                             buffer_sizes);
  set_borders(array, bv, bs);
}

void set_borders(Array    &array,
                 float     border_values,
                 float     buffer_sizes,
                 glm::vec4 bbox)
{
  glm::vec4 bv = glm::vec4(border_values,
                           border_values,
                           border_values,
                           border_values);
  glm::vec4 bs = glm::vec4(buffer_sizes,
                           buffer_sizes,
                           buffer_sizes,
                           buffer_sizes);
  set_borders(array, bv, bs, bbox);
}

} // namespace hmap
