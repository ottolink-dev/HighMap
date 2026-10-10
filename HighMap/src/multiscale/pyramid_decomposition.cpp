/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/export/image.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/core.hpp"
#include "highmap/multiscale/pyramid.hpp"
#include "highmap/operator.hpp"

namespace hmap
{

PyramidDecomposition::PyramidDecomposition(Array &array, int nlevels_)
    : nlevels(nlevels_), p_array(&array)
{
  if (!validate_non_empty(array)) return;

  // check and/or adjust number of levels
  int np2 = std::min(highest_power_of_2(array.shape.x),
                     highest_power_of_2(array.shape.y));

  if (nlevels_ > np2)
    this->nlevels = np2;
  else if (nlevels_ <= 0)
    this->nlevels = np2 + nlevels_;

  if (this->nlevels != nlevels_)
    hmap::log::trace("number of levels adjusted, effective / maximum: {} / {}, "
                     "requested: {}",
                     this->nlevels,
                     np2,
                     nlevels_);

  // default filter is a Laplace filter
  this->low_pass_filter_function = [](const Array &input)
  {
    Array output = input;
    {
      float sigma = 0.25f;
      int   iterations = 1;
      laplace(output, sigma, iterations);
    }
    return output;
  };
}

void PyramidDecomposition::decompose()
{
  if (!p_array || !validate_non_empty(*p_array)) return;

  // reset any existing decomposition
  this->components.clear();

  // working array
  glm::ivec2 level_shape = this->p_array->shape;
  Array      array_low = *this->p_array;

  for (int n = 0; n < this->nlevels; n++)
  {
    // filtering
    Array array_filtered = this->low_pass_filter_function(array_low);

    // store high-pass component
    this->components.push_back(array_low - array_filtered);

    // keep iterating on the low-pass component
    array_low = array_filtered;

    if (n == this->nlevels - 1)
    {
      this->residual = array_low;
      break;
    }

    // downscale and move on to the next shape (use bilinear
    // interpolation even when downscaling to limit field stretching)
    level_shape /= 2;
    array_low = array_low.resample_to_shape(level_shape);
  }
}

Array PyramidDecomposition::reconstruct()
{
  Array array_out = this->residual;

  glm::ivec2 shape = array_out.shape;

  for (int n = this->nlevels; n-- > 0;)
  {
    array_out += this->components[n];

    if (n > 0)
    {
      shape = glm::ivec2(shape.x * 2, shape.y * 2);
      array_out = array_out.resample_to_shape_bicubic(shape);
    }
  }

  return array_out;
}

void PyramidDecomposition::to_png(std::string fname, int cmap, bool hillshading)
{
  // highpass and lowpass storage
  std::vector<Array> banner_arrays = {};

  glm::ivec2 shape_ref = this->p_array->shape;

  Array      array_out = this->residual;
  glm::ivec2 shape = array_out.shape;

  for (int n = this->nlevels; n-- > 0;)
  {
    Array lphp = vstack(
        this->components[n].resample_to_shape_nearest(shape_ref),
        array_out.resample_to_shape_nearest(shape_ref));

    banner_arrays.insert(banner_arrays.begin(), lphp);

    array_out += this->components[n];

    if (n > 0)
    {
      shape = {shape.x * 2, shape.y * 2};
      array_out = array_out.resample_to_shape_bicubic(shape);
    }
  }
  export_banner_png(fname, banner_arrays, cmap, hillshading);
}

Array PyramidDecomposition::transform(
    std::function<Array(const Array &, const int current_level)> function,
    int                                                          support,
    std::vector<float>                                           level_weights,
    int                                                          finest_level)
{
  // if no weights are provided, just a constant one
  if (!level_weights.size())
  {
    level_weights.resize(this->nlevels);
    std::fill(level_weights.begin(), level_weights.end(), 1);
  }

  Array      array_out = this->residual;
  glm::ivec2 shape = array_out.shape;

  for (int n = this->nlevels; n-- > finest_level;)
  {
    switch (support)
    {
    case pyramid_transform_support::FULL:
    {
      array_out += this->components[n];
      Array component_transformed = function(array_out, n);
      array_out = lerp(array_out, component_transformed, level_weights[n]);
    }
    break;

    case pyramid_transform_support::HIGHPASS_ONLY:
    {
      Array component_transformed = function(this->components[n], n);
      array_out += lerp(this->components[n],
                        component_transformed,
                        level_weights[n]);
    }
    break;

    case pyramid_transform_support::LOWPASS_ONLY:
    {
      Array component_transformed = function(array_out, n);
      array_out = lerp(array_out + this->components[n],
                       component_transformed + this->components[n],
                       level_weights[n]);
    }
    break;

    default:
      hmap::log::error("unknown support");
      throw std::runtime_error("unknown support");
    }

    if (n > 0)
    {
      shape = {shape.x * 2, shape.y * 2};
      array_out = array_out.resample_to_shape_bicubic(shape);
    }
  }

  return array_out;
}

// --- HELPER

int highest_power_of_2(int n)
{
  int p = (int)std::log2(n);
  return p;
}
} // namespace hmap
