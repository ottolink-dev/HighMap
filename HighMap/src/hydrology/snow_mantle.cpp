/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/curvature.hpp"
#include "highmap/erosion/thermal_erosion.hpp"
#include "highmap/gradient.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/core.hpp"
#include "highmap/primitives/coherent_noise.hpp"
#include "highmap/range.hpp"

namespace hmap::gpu
{

Array snow_mantle(const Array &z,
                  float        strength,
                  float        snowline,
                  float        falloff,
                  float        talus_coverage,
                  float        talus_threshold,
                  float        directional,
                  float        wind_direction,
                  float        peaks,
                  float        melt,
                  float        sun_azimuth,
                  float        sun_zenith,
                  int          thermal_iterations,
                  float        thermal_talus,
                  float        coverage,
                  const Array *p_existing_snow,
                  const Array *p_mask,
                  const Array *p_noise)
{
  if (!validate_non_empty(z)) return Array();
  if (p_existing_snow && !validate_same_shape(z, *p_existing_snow))
    return Array();
  if (p_mask && !validate_same_shape(z, *p_mask)) return Array();
  if (p_noise && !validate_same_shape(z, *p_noise)) return Array();

  const glm::ivec2 shape = z.shape;
  const int        n_pixels = shape.x * shape.y;

  // --- Altitude accumulation factor: smoothstep transition around snowline

  const float alt_min = snowline - falloff;
  const float alt_max = snowline + falloff;

  Array altitude_factor(shape);

  if (alt_max <= alt_min + 1e-9f)
  {
    for (int k = 0; k < n_pixels; ++k)
      altitude_factor.vector[k] = z.vector[k] >= snowline ? 1.f : 0.f;
  }
  else
  {
    for (int k = 0; k < n_pixels; ++k)
      altitude_factor.vector[k] = threshold_smooth(z.vector[k],
                                                   alt_min,
                                                   alt_max);
  }

  // --- Talus exclusion factor: steep talus slopes shed snow

  Array talus_factor(shape, 1.f);

  if (talus_coverage > 0.f)
  {
    Array talus = gradient_talus(z);

    const float talus_upper = talus_threshold * 1.5f + 1e-9f;

    for (int k = 0; k < n_pixels; ++k)
    {
      float s = threshold_smooth(talus.vector[k], talus_threshold, talus_upper);
      talus_factor.vector[k] = 1.f - talus_coverage * s;
    }
  }

  // --- Directional wind bias

  Array wind_factor(shape, 1.f);

  if (std::abs(directional) > 1e-5f)
  {
    const float wind_rad = M_PI * wind_direction / 180.f;
    const float cos_w = std::cos(wind_rad);
    const float sin_w = std::sin(wind_rad);

    Array dx = gradient_x(z) * float(shape.x);
    Array dy = gradient_y(z) * float(shape.y);

    for (int k = 0; k < n_pixels; ++k)
    {
      float gx = dx.vector[k];
      float gy = dy.vector[k];
      float norm = std::sqrt(gx * gx + gy * gy + 1.f);

      // Unit normal (pointing out of terrain)
      float nx = -gx / norm;
      float ny = -gy / norm;

      // Leeward accumulation / windward erosion
      float dot_wind = -(nx * cos_w + ny * sin_w);
      float w = 1.f + directional * dot_wind;

      wind_factor.vector[k] = std::max(0.f, w);
    }
  }

  // --- Peaks / Ridges bonus accumulation

  Array peak_factor(shape, 1.f);

  if (peaks > 0.f)
  {
    Array curvature = gpu::curvature_quadric(z, 2, CurvatureType::CT_MEAN);

    for (int k = 0; k < n_pixels; ++k)
    {
      // Positive convex curvature (peaks and sharp ridges)
      float conv = std::max(0.f, -curvature.vector[k]);
      peak_factor.vector[k] = 1.f + peaks * std::min(conv * 10.f, 2.f);
    }
  }

  // --- Base accumulation layer

  Array b(shape);

  for (int k = 0; k < n_pixels; ++k)
  {
    b.vector[k] = altitude_factor.vector[k] * talus_factor.vector[k] *
                  wind_factor.vector[k] * peak_factor.vector[k];
  }

  // --- Solar melting

  if (melt > 0.f)
  {
    Array melt_map = snow_melting_map(z,
                                      snowline - falloff,
                                      snowline + falloff,
                                      1.f,
                                      1.f,
                                      sun_azimuth,
                                      sun_zenith,
                                      1.f,
                                      1.f,
                                      0.5f);

    for (int k = 0; k < n_pixels; ++k)
    {
      float m = std::clamp(melt_map.vector[k], 0.f, 1.f);
      b.vector[k] *= std::max(0.f, 1.f - melt * m);
    }
  }

  // --- Micro-texture / Noise modulation

  if (p_noise)
  {
    for (int k = 0; k < n_pixels; ++k)
      b.vector[k] *= std::max(0.f, p_noise->vector[k]);
  }

  // --- Coverage shaping

  if (coverage > 0.f && coverage < 1.f)
  {
    const float cov_min = 1.f - coverage;

    for (int k = 0; k < n_pixels; ++k)
      b.vector[k] = threshold_smooth(b.vector[k], cov_min, 1.f);
  }

  // --- Apply global strength

  b *= strength;

  // --- Downhill redistribution via mass-preserving GPU thermal transport

  if (thermal_iterations > 0)
  {
    Array z_wrk = z + b;
    Array snow_mask = b / std::max(b.max(), 1e-6f);
    Array talus_map(shape, thermal_talus);

    gpu::thermal_conserve(z_wrk,
                          &snow_mask,
                          talus_map,
                          thermal_iterations,
                          0.5f,
                          &z);

    b = z_wrk - z;
    clamp_min(b, 0.f);
  }

  // --- Combine with existing snow if provided

  if (p_existing_snow)
  {
    for (int k = 0; k < n_pixels; ++k)
    {
      float ex = p_existing_snow->vector[k];
      b.vector[k] = ex + (1.f - std::min(1.f, ex)) * b.vector[k];
    }
  }

  // --- Spatial mask

  if (p_mask) b *= *p_mask;

  clamp_min(b, 0.f);

  return b;
}

} // namespace hmap::gpu
