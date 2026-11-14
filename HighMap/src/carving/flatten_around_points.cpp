/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>

#include "highmap/array.hpp"
#include "highmap/carving.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/profiles.hpp"

namespace hmap
{

void flatten_around_points(Array                    &z,
                           const std::vector<float> &xp,
                           const std::vector<float> &yp,
                           const std::vector<float> &zp,
                           const std::vector<float> &radius,
                           float                     strength,
                           RadialProfile             radial_profile,
                           float                     profile_param,
                           Array                    *p_mask,
                           glm::vec4                 bbox)
{
  if (!validate_non_empty(z)) return;
  if (!validate_non_empty(xp)) return;
  if (!validate_non_empty(yp)) return;
  if (!validate_non_empty(zp)) return;
  if (!validate_non_empty(radius)) return;

  if (xp.size() != yp.size() || xp.size() != zp.size() ||
      xp.size() != radius.size())
    return;

  if (strength <= 0.0f) return;

  const int nx = z.shape.x;
  const int ny = z.shape.y;

  // cell size in bbox units (guard against single-cell axes)
  const float dx = nx > 1 ? (bbox.y - bbox.x) / float(nx - 1) : 1.f;
  const float dy = ny > 1 ? (bbox.w - bbox.z) / float(ny - 1) : 1.f;

  auto profile_fct = get_radial_profile_function(radial_profile, profile_param);

  // per-cell accumulators
  Array weight_sum(z.shape, 0.f);
  Array z_weighted_sum(z.shape, 0.f);
  Array falloff_sum(z.shape, 0.f);

  for (size_t k = 0; k < xp.size(); ++k)
  {
    const float r = radius[k];

    if (r <= 0.f) continue;

    // pixel bounding box of the point influence
    const int i_min = std::max(0, int(std::floor((xp[k] - r - bbox.x) / dx)));
    const int i_max = std::min(nx - 1,
                               int(std::ceil((xp[k] + r - bbox.x) / dx)));
    const int j_min = std::max(0, int(std::floor((yp[k] - r - bbox.z) / dy)));
    const int j_max = std::min(ny - 1,
                               int(std::ceil((yp[k] + r - bbox.z) / dy)));

    for (int j = j_min; j <= j_max; ++j)
      for (int i = i_min; i <= i_max; ++i)
      {
        const float ddx = bbox.x + float(i) * dx - xp[k];
        const float ddy = bbox.z + float(j) * dy - yp[k];
        const float u = std::sqrt(ddx * ddx + ddy * ddy) / r;

        if (u > 1.f) continue;

        // maximum flattening at the point, tapering to 0 at its radius
        const float falloff = 1.f - profile_fct(u);

        // larger radius, larger influence
        const float weight = falloff * r;

        weight_sum(i, j) += weight;
        z_weighted_sum(i, j) += weight * zp[k];
        falloff_sum(i, j) += falloff;
      }
  }

  // --- Blend toward the weighted reference elevation

  if (p_mask) *p_mask = Array(z.shape, 0.f);

  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
    {
      if (weight_sum(i, j) <= 0.f) continue;

      const float z_ref = z_weighted_sum(i, j) / weight_sum(i, j);
      const float w = std::clamp(strength * falloff_sum(i, j), 0.f, 1.f);

      z(i, j) += w * (z_ref - z(i, j));

      if (p_mask) (*p_mask)(i, j) = w;
    }
}

void flatten_around_points(Array                    &z,
                           const std::vector<float> &xp,
                           const std::vector<float> &yp,
                           const std::vector<float> &zp,
                           float                     radius,
                           float                     strength,
                           RadialProfile             radial_profile,
                           float                     profile_param,
                           Array                    *p_mask,
                           glm::vec4                 bbox)
{
  const std::vector<float> radii(xp.size(), radius);

  flatten_around_points(z,
                        xp,
                        yp,
                        zp,
                        radii,
                        strength,
                        radial_profile,
                        profile_param,
                        p_mask,
                        bbox);
}

} // namespace hmap
