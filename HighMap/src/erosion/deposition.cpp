/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/array.hpp"
#include "highmap/erosion/strata_erosion.hpp"
#include "highmap/erosion/thermal_erosion.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/range.hpp"

#define SPAWN_LOW_LIMIT 0.1f
#define GRADIENT_MIN 0.0001f
#define SEDIMENT_MIN 0.001f
#define VELOCITY_MIN 0.001f

namespace hmap::gpu
{

void sediment_deposition(Array       &z,
                         const Array &talus,
                         Array       *p_deposition_map,
                         float        max_deposition,
                         int          iterations,
                         int          thermal_subiterations)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, talus)) return;

  float deposition_step = max_deposition / (int)iterations;
  Array smap = Array(z.shape); // sediment map

  for (int it = 0; it < iterations; it++)
  {
    smap = smap + deposition_step;
    Array z_tot = z + smap;

    gpu::thermal(z_tot, talus, thermal_subiterations);
    smap = maximum(z_tot - z, 0.f);
  }
  z = z + smap;

  if (p_deposition_map) *p_deposition_map = smap;
}

void sediment_deposition(Array       &z,
                         const Array *p_mask,
                         const Array &talus,
                         Array       *p_deposition_map,
                         float        max_deposition,
                         int          iterations,
                         int          thermal_subiterations)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, talus)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;

  if (!p_mask)
    sediment_deposition(z,
                        talus,
                        p_deposition_map,
                        max_deposition,
                        iterations,
                        thermal_subiterations);
  else
  {
    Array z_f = z;
    sediment_deposition(z_f,
                        talus,
                        p_deposition_map,
                        max_deposition,
                        iterations,
                        thermal_subiterations);
    z = lerp(z, z_f, *(p_mask));
  }
}

void sediment_layer(Array       &z,
                    const Array &talus_layer,
                    const Array &talus_upper_limit,
                    int          iterations,
                    bool         apply_post_filter,
                    Array       *p_deposition_map)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, talus_layer) ||
      !validate_same_shape(z, talus_upper_limit))
    return;

  // backup input
  Array z_bckp = z;

  // prepare talus
  Array g_talus = gradient_talus(z);
  Array talus_ref = talus_layer;
  Array fmask = Array(z.shape, 1.f);

  for (int j = 0; j < z.shape.y; j++)
    for (int i = 0; i < z.shape.x; i++)
      if (g_talus(i, j) > talus_upper_limit(i, j))
      {
        fmask(i, j) = 0.f;
        talus_ref(i, j) = 4.f * (g_talus(i, j) - talus_upper_limit(i, j)) +
                          talus_upper_limit(i, j);
      }

  // apply thermal erosion with prepared inputs
  Array sediment_layer_map(z.shape);

  gpu::thermal(z, talus_ref, iterations, nullptr, &sediment_layer_map);

  // smooth transitions
  if (apply_post_filter)
  {
    // hmap::laplace(sediment_layer_map);
    z = z_bckp + sediment_layer_map;

    // filter also the deposition layer only (and also filter the mask
    // itself to avoid kinky spatial transitions)
    hmap::laplace(fmask);
    hmap::laplace(z, &fmask);
  }

  // output deposition map
  if (p_deposition_map) *p_deposition_map = sediment_layer_map;
}

} // namespace hmap::gpu
