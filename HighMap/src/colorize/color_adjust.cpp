/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/colorize.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/texture.hpp"

namespace hmap
{

inline uint32_t helper_hash_u32(uint32_t x)
{
  x ^= x >> 16;
  x *= 0x7feb352d;
  x ^= x >> 15;
  x *= 0x846ca68b;
  x ^= x >> 16;
  return x;
}

inline float helper_hash_01(uint32_t x)
{
  return float(helper_hash_u32(x)) * (1.0f / float(0xffffffffu));
}

glm::vec3 shade_color(const glm::vec3 &color, float t)
{
  t = std::clamp(t, 0.f, 1.f);
  const float     luma = glm::dot(color, glm::vec3(0.299f, 0.587f, 0.114f));
  const glm::vec3 grey(luma);
  const glm::vec3 lit = color * (1.f + 0.4f * t);
  return glm::mix(grey * (1.f + 0.4f * t), lit, 1.f - 0.3f * t);
}

void color_adjust(Array      &r,
                  Array      &g,
                  Array      &b,
                  ColorAdjust param,
                  glm::ivec2  dither_seed)
{
  if (!validate_non_empty(r) || !validate_same_shape(r, g) ||
      !validate_same_shape(r, b))
    return;

  const glm::ivec2 &shape = r.shape;
  const float       inv_gamma = 1.f / param.gamma;

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      glm::vec3 c{r(i, j), g(i, j), b(i, j)};

      // levels
      c = (c - param.in_min) / (param.in_max - param.in_min);
      c = glm::clamp(c, 0.0f, 1.0f);

      // exposure
      c *= std::exp2(param.exposure);

      // tonemap
      if (param.filmic_tonemap) c = c / (glm::vec3(1.0f) + c); // Reinhard

      if (param.aces_tonemap)
      {
        glm::vec3 num = c * (2.51f * c + 0.03f);
        glm::vec3 den = c * (2.43f * c + 0.59f) + 0.14f;
        c = num / den;
      }

      if (param.agx_tonemap)
      {
        // scene-linear → log2
        glm::vec3 x = glm::max(c, glm::vec3(1e-6f));
        x = glm::log2(x);

        // normalize log range (AGX working range)
        constexpr float log_min = -10.0f;
        constexpr float log_max = 0.0f;

        x = (x - log_min) / (log_max - log_min);
        x = glm::clamp(x, 0.0f, 1.0f);

        // AGX contrast curve (smooth filmic sigmoid)
        x = x * x * (3.0f - 2.0f * x);

        // back to linear
        c = glm::exp2(x * (log_max - log_min) + log_min);

        // Gamut / saturation compression (very important for AGX look)
        float luma = glm::dot(c, glm::vec3(0.2126f, 0.7152f, 0.0722f));
        float sat = glm::mix(1.0f, 0.85f, luma);
        c = glm::mix(glm::vec3(luma), c, sat);
      }

      // contrast
      c = (c - 0.5f) * param.contrast + 0.5f;

      // saturation
      float l = glm::dot(c, glm::vec3(0.2126f, 0.7152f, 0.0722f));
      c = glm::mix(glm::vec3(l), c, param.saturation);

      // temperature
      c.r *= 1.0f + param.temperature;
      c.b *= 1.0f - param.temperature;

      // shade
      if (param.shade != 0.0f) c = shade_color(c, param.shade);

      // gamma
      c = glm::pow(c, glm::vec3(inv_gamma));

      // hash dithering
      int x = dither_seed.x * shape.x + i;
      int y = dither_seed.y * shape.y + j;

      float n = helper_hash_01(uint32_t(x) ^ (uint32_t(y) << 16));
      n = (n - 0.5f) * param.dither_amp;

      c += glm::vec3(n);

      // final clamp
      c = glm::clamp(c, 0.0f, 1.0f);

      r(i, j) = c.r;
      g(i, j) = c.g;
      b(i, j) = c.b;
    }
}

void shade_color(Texture &tex, const Array &factor)
{
  if (!validate_non_empty(tex, 3)) return;
  if (!validate_non_empty(factor)) return;
  if (!validate_same_shape(tex, factor)) return;

  const glm::ivec2 &shape = tex.shape;

#pragma omp parallel for collapse(2)
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      const glm::vec3 c{tex(i, j, 0), tex(i, j, 1), tex(i, j, 2)};
      const glm::vec3 shaded = shade_color(c, factor(i, j));

      tex(i, j, 0) = shaded.r;
      tex(i, j, 1) = shaded.g;
      tex(i, j, 2) = shaded.b;
    }
  }
}

void shade_topographic_residual(Texture     &tex,
                                const Array &z,
                                int          ir,
                                float        strength)
{
  if (!validate_non_empty(tex, 3)) return;
  if (!validate_non_empty(z)) return;
  if (!validate_same_shape(tex, z)) return;

  const Array residual = topographic_residual(z, ir, true);
  Array       factor(z.shape, 0.5f);

  const int64_t n = static_cast<int64_t>(residual.vector.size());

#pragma omp parallel for
  for (int64_t idx = 0; idx < n; ++idx)
  {
    const float val = 0.5f + 0.5f * strength * residual.vector[idx];
    factor.vector[idx] = std::clamp(val, 0.0f, 1.0f);
  }

  shade_color(tex, factor);
}

} // namespace hmap
