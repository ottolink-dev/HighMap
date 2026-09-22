/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/operator.hpp"
#include "highmap/range.hpp"

namespace hmap
{

// --- Canyonize Implementation

void canyonize(Array        &array,
               std::uint32_t seed,
               int           nlevels,
               float         convex_ratio,
               float         gamma_convex,
               float         gamma_concave,
               float         clamp_min_val,
               float         k_smooth,
               float         noise_ratio,
               const Array  *p_noise,
               float         vmin,
               float         vmax)
{
  if (!validate_non_empty(array)) return;
  if (p_noise && !validate_same_shape(array, *p_noise)) return;

  // clamp min (smooth) to flatten the canyon bottom
  if (clamp_min_val > -std::numeric_limits<float>::infinity())
  {
    if (k_smooth > 0.f)
      clamp_min_smooth(array, clamp_min_val, k_smooth);
    else
      clamp_min(array, clamp_min_val);
  }

  // redefine min/max if sentinels values are detected
  if (vmax <= vmin)
  {
    vmin = array.min();
    vmax = array.max();
  }

  // avoid division by zero or degenerate ranges
  if (vmax <= vmin || nlevels <= 0) return;

  std::mt19937                          gen(seed);
  std::uniform_real_distribution<float> dis(-noise_ratio, noise_ratio);

  // compute nominal interval widths with convex/concave ratio
  float              c_ratio = std::max(convex_ratio, 1e-3f);
  std::vector<float> widths(nlevels);
  float              total_weight = 0.f;
  for (int i = 0; i < nlevels; ++i)
  {
    float w = (i % 2 == 0) ? c_ratio : 1.f;
    widths[i] = w;
    total_weight += w;
  }

  // build cumulative strata elevation levels
  float              total_range = vmax - vmin;
  std::vector<float> levels(nlevels + 1);
  levels[0] = vmin;
  float current_elev = vmin;
  for (int i = 0; i < nlevels; ++i)
  {
    float step = (widths[i] / total_weight) * total_range;
    current_elev += step;
    levels[i + 1] = current_elev;
  }
  levels.back() = vmax;

  // add noise to fluctuate strata elevations, preserving outer boundaries
  for (size_t k = 1; k < levels.size() - 1; ++k)
  {
    float avg_step = total_range / static_cast<float>(nlevels);
    levels[k] += dis(gen) * avg_step;
  }

  // sort to guarantee strictly non-decreasing monotonic intervals
  std::sort(levels.begin(), levels.end());
  levels.front() = vmin;
  levels.back() = vmax;

  // stratification with alternating power laws starting with convex
  // add noise -> apply filter -> subtract noise (with clamping to prevent range
  // escape)
  auto lambda =
      [&levels, vmin, vmax, gamma_convex, gamma_concave](float x,
                                                         float noise = 0.f)
  {
    // 1. Add noise to input elevation
    float y = x + noise;
    y = std::clamp(y, vmin, vmax);

    // 2. Find level interval [levels[n], levels[n+1]]
    size_t n = 1;
    while (n < levels.size() && y > levels[n])
      n++;
    n--;

    if (n >= levels.size() - 1) n = levels.size() - 2;

    float l_min = levels[n];
    float l_max = levels[n + 1];
    float interval_range = l_max - l_min;

    if (interval_range <= 0.f) return x;

    // normalized elevation within interval [0, 1]
    float u = std::clamp((y - l_min) / interval_range, 0.f, 1.f);

    // 1st interval (n = 0) starts convex (> 1), then alternates to concave (<
    // 1), etc.
    float gamma = (n % 2 == 0) ? gamma_convex : gamma_concave;
    if (gamma <= 0.f) gamma = 1e-4f;

    u = std::pow(u, gamma);

    // 3. Map filtered elevation within interval
    float filtered_y = l_min + u * interval_range;

    // 4. Subtract noise to revert coordinate modulation
    float out = filtered_y - noise;

    // Clamp to valid elevation range
    return std::clamp(out, vmin, vmax);
  };

  if (p_noise)
  {
    std::transform(array.vector.begin(),
                   array.vector.end(),
                   p_noise->vector.begin(),
                   array.vector.begin(),
                   lambda);
  }
  else
  {
    std::transform(array.vector.begin(),
                   array.vector.end(),
                   array.vector.begin(),
                   [lambda](float x) { return lambda(x, 0.f); });
  }
}

void canyonize(Array        &array,
               std::uint32_t seed,
               int           nlevels,
               const Array  *p_mask,
               float         convex_ratio,
               float         gamma_convex,
               float         gamma_concave,
               float         clamp_min_val,
               float         k_smooth,
               float         noise_ratio,
               const Array  *p_noise,
               float         vmin,
               float         vmax)
{
  if (!validate_non_empty(array)) return;
  if (p_noise && !validate_same_shape(array, *p_noise)) return;
  if (p_mask && !validate_same_shape(array, *p_mask)) return;

  apply_with_mask(array,
                  p_mask,
                  [&](Array &a)
                  {
                    canyonize(a,
                              seed,
                              nlevels,
                              convex_ratio,
                              gamma_convex,
                              gamma_concave,
                              clamp_min_val,
                              k_smooth,
                              noise_ratio,
                              p_noise,
                              vmin,
                              vmax);
                  });
}

} // namespace hmap
