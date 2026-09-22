/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/functions.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/math/array.hpp"
#include "highmap/operator.hpp"
#include "highmap/primitives/coherent_noise.hpp"
#include "highmap/range.hpp"
#include "highmap/selector.hpp"
#include "highmap/transform.hpp"

namespace hmap
{

void recast_billow(Array &array, float vref, float k)
{
  if (!validate_non_empty(array)) return;
  array = 2.f * (vref + abs_smooth(array - vref, k)) - 1.f;
}

void recast_canyon(Array &array, const Array &vcut, float gamma)
{
  if (!validate_non_empty(array) || !validate_same_shape(array, vcut)) return;

  auto lambda = [&gamma](float a, float b)
  { return a > b ? a : b * std::pow(a / b, gamma); };

  std::transform(array.vector.begin(),
                 array.vector.end(),
                 vcut.vector.begin(),
                 array.vector.begin(),
                 lambda);
}

void recast_canyon(Array       &array,
                   const Array &vcut,
                   const Array *p_mask,
                   float        gamma)
{
  if (!validate_non_empty(array) || !validate_same_shape(array, vcut)) return;
  if (p_mask && !validate_same_shape(array, *p_mask)) return;

  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { recast_canyon(a, vcut, gamma); });
}

void recast_canyon(Array &array, float vcut, float gamma, const Array *p_noise)
{
  if (!validate_non_empty(array)) return;
  if (p_noise && !validate_same_shape(array, *p_noise)) return;

  if (!p_noise)
  {
    auto lambda = [&vcut, &gamma](float a)
    { return a > vcut ? a : vcut * std::pow(a / vcut, gamma); };

    std::transform(array.vector.begin(),
                   array.vector.end(),
                   array.vector.begin(),
                   lambda);
  }
  else
  {
    auto lambda = [&vcut, &gamma](float a, float b) {
      return a > (vcut + b) ? a : (vcut + b) * std::pow(a / (vcut + b), gamma);
    };

    std::transform(array.vector.begin(),
                   array.vector.end(),
                   (*p_noise).vector.begin(),
                   array.vector.begin(),
                   lambda);
  }
}

void recast_canyon(Array       &array,
                   float        vcut,
                   const Array *p_mask,
                   float        gamma,
                   const Array *p_noise)
{
  if (!validate_non_empty(array)) return;
  if (p_noise && !validate_same_shape(array, *p_noise)) return;
  if (p_mask && !validate_same_shape(array, *p_mask)) return;

  if (!p_mask)
    recast_canyon(array, vcut, gamma, p_noise);
  else
  {
    Array array_f = array;
    recast_canyon(array_f, vcut, gamma, p_noise);
    array = lerp(array, array_f, *(p_mask));
  }
}

void recast_cracks(Array &array,
                   float  cut_min,
                   float  cut_max,
                   float  k_smoothing,
                   float  vmin,
                   float  vmax)
{
  if (!validate_non_empty(array)) return;

  // redefine min/max if sentinels values are detected
  if (vmax < vmin)
  {
    vmin = array.min();
    vmax = array.max();
  }

  remap(array, 0.f, 1.f, vmin, vmax);

  Array z1 = array - cut_max;
  Array z2 = cut_max - array;

  array = maximum_smooth(z1, z2, k_smoothing);
  array = minimum_smooth(array, Array(array.shape, cut_min), k_smoothing);
  array /= cut_min;
}

void recast_escarpment(Array &array,
                       int    ir,
                       float  ratio,
                       float  scale,
                       bool   reverse,
                       bool   transpose_effect,
                       float  global_scaling)
{
  if (!validate_non_empty(array)) return;

  if (transpose_effect) array = transpose(array);

  if (global_scaling == 0.f)
    global_scaling = 20.f * (array.max() - array.min()) / array.shape.x;

  // cumulated displacement /x
  Array cdx(array.shape);

  if (!reverse)
  {
    for (int j = 0; j < array.shape.y; j++)
      for (int i = 1; i < array.shape.x; i++)
      {
        float v = array(i, j) > array(i - 1, j) ? -ratio : 1.f;
        cdx(i, j) = std::min(cdx(i - 1, j) + v, 0.f);
      }
    cdx *= scale * global_scaling;
  }
  else
  {
    for (int j = 0; j < array.shape.y; j++)
      for (int i = array.shape.x - 2; i > -1; i--)
      {
        float v = array(i, j) > array(i + 1, j) ? -ratio : 1.f;
        cdx(i, j) = std::min(cdx(i + 1, j) + v, 0.f);
      }
    cdx.infos();
    cdx *= -scale * global_scaling;
  }

  smooth_flat(cdx, ir);

  // warp
  warp(array, &cdx, nullptr);

  if (transpose_effect) array = transpose(array);
}

void recast_escarpment(Array       &array,
                       const Array *p_mask,
                       int          ir,
                       float        ratio,
                       float        scale,
                       bool         reverse,
                       bool         transpose_effect,
                       float        global_scaling)
{
  if (!validate_non_empty(array)) return;
  if (p_mask && !validate_same_shape(array, *p_mask)) return;

  apply_with_mask(array,
                  p_mask,
                  [&](Array &a)
                  {
                    recast_escarpment(a,
                                      ir,
                                      ratio,
                                      scale,
                                      reverse,
                                      transpose_effect,
                                      global_scaling);
                  });
}

void recast_peak(Array &array, int ir, float gamma, float k)
{
  if (!validate_non_empty(array)) return;

  Array ac = array;
  smooth_cpulse(ac, ir);
  array = maximum_smooth(array, ac, k);
  clamp_min(array, 0.f);
  array = ac * pow(array, gamma);
}

void recast_peak(Array       &array,
                 int          ir,
                 const Array *p_mask,
                 float        gamma,
                 float        k)
{
  if (!validate_non_empty(array)) return;
  if (p_mask && !validate_same_shape(array, *p_mask)) return;

  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { recast_peak(a, ir, gamma, k); });
}

void recast_rocky_slopes(Array        &array,
                         float         talus,
                         int           ir,
                         float         amplitude,
                         std::uint32_t seed,
                         float         kw,
                         float         gamma,
                         const Array  *p_noise,
                         glm::vec4     bbox)
{
  if (!validate_non_empty(array)) return;
  if (p_noise && !validate_same_shape(array, *p_noise)) return;

  // slope-based criteria
  Array c = select_gradient_binary(array, talus);
  smooth_cpulse(c, ir);

  if (!p_noise)
  {
    Array noise = noise_fbm(NoiseType::PERLIN,
                            array.shape,
                            glm::vec2(kw, kw),
                            seed,
                            8,
                            0.f,
                            0.5f,
                            2.f,
                            nullptr,
                            nullptr,
                            nullptr,
                            bbox);
    gamma_correction_local(noise, gamma, ir, 0.1f);
    {
      int ir2 = (int)(ir / 4.f);
      if (ir2 > 1) gamma_correction_local(noise, gamma, ir2, 0.1f);
    }
    array += amplitude * noise * c;
  }
  else
    array += amplitude * (*p_noise) * c;
}

void recast_rocky_slopes(Array        &array,
                         float         talus,
                         int           ir,
                         float         amplitude,
                         std::uint32_t seed,
                         float         kw,
                         const Array  *p_mask,
                         float         gamma,
                         const Array  *p_noise,
                         glm::vec4     bbox)
{
  if (!validate_non_empty(array)) return;
  if (p_noise && !validate_same_shape(array, *p_noise)) return;
  if (p_mask && !validate_same_shape(array, *p_mask)) return;

  apply_with_mask(array,
                  p_mask,
                  [&](Array &a)
                  {
                    recast_rocky_slopes(a,
                                        talus,
                                        ir,
                                        amplitude,
                                        seed,
                                        kw,
                                        gamma,
                                        p_noise,
                                        bbox);
                  });
}

void recast_sag(Array &array, float vref, float k)
{
  if (!validate_non_empty(array)) return;
  array = 0.5f * array + vref - abs_smooth(array - vref, k);
}

void recast_sag(Array &array, float vref, float k, const Array *p_mask)
{
  if (!validate_non_empty(array)) return;
  if (p_mask && !validate_same_shape(array, *p_mask)) return;

  apply_with_mask(array, p_mask, [&](Array &a) { recast_sag(a, vref, k); });
}

} // namespace hmap
