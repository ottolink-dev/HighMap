/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/array.hpp"
#include "highmap/erosion/strata_erosion.hpp"
#include "highmap/erosion/thermal_erosion.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/range.hpp"

namespace hmap::gpu
{

void valley_fill(Array       &z,
                 const Array &talus,
                 int          iterations,
                 float        gamma,
                 float        ratio,
                 float        zmin,
                 float        zmax,
                 float        elevation_max_ratio,
                 bool         preserve_elevation_range,
                 const Array *p_noise,
                 Array       *p_deposition_map)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, talus)) return;
  if (p_noise && !validate_same_shape(z, *p_noise)) return;

  if (zmax <= zmin)
  {
    zmin = z.min();
    zmax = z.max();
  }

  Array z_bckp;
  if (p_deposition_map) z_bckp = z;

  // add noise if any
  const Array *p_talus = &talus;
  Array        talus_with_noise;

  if (p_noise)
  {
    talus_with_noise = talus * (1.f + *p_noise);
    clamp_min(talus_with_noise, 0.f);
    p_talus = &talus_with_noise;
  }

  // scree deposition
  Array ze = z;
  gpu::thermal_scree(ze,
                     *p_talus,
                     Array(z.shape, elevation_max_ratio * zmax),
                     iterations);

  // mixing mask
  Array t = z;
  remap(t, 0.f, 1.f, zmin, zmax);
  gamma_correction(t, gamma);
  remap(t, 1.f - ratio, 1.f);

  // apply
  z = lerp(ze, z, t);

  if (p_deposition_map) *p_deposition_map = z - z_bckp;

  if (preserve_elevation_range) remap(z, zmin, zmax);
}

void valley_fill(Array       &z,
                 const Array *p_mask,
                 const Array &talus,
                 int          iterations,
                 float        gamma,
                 float        ratio,
                 float        zmin,
                 float        zmax,
                 float        elevation_max_ratio,
                 bool         preserve_elevation_range,
                 const Array *p_noise,
                 Array       *p_deposition_map)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, talus)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;
  if (p_noise && !validate_same_shape(z, *p_noise)) return;

  if (!p_mask)
  {
    gpu::valley_fill(z,
                     talus,
                     iterations,
                     gamma,
                     ratio,
                     zmin,
                     zmax,
                     elevation_max_ratio,
                     preserve_elevation_range,
                     p_noise,
                     p_deposition_map);
  }
  else
  {
    Array z_f = z;
    gpu::valley_fill(z_f,
                     talus,
                     iterations,
                     gamma,
                     ratio,
                     zmin,
                     zmax,
                     elevation_max_ratio,
                     preserve_elevation_range,
                     p_noise,
                     p_deposition_map);
    z = lerp(z, z_f, *(p_mask));
  }
}

} // namespace hmap::gpu
