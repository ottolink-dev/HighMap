/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include "highmap/array.hpp"
#include "highmap/erosion/deprecated.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/range.hpp"

namespace hmap
{

void hydraulic_algebric(Array &z,
                        float  talus_ref,
                        int    ir,
                        Array *p_bedrock,
                        Array *p_erosion_map,
                        Array *p_deposition_map,
                        float  c_erosion,
                        float  c_deposition,
                        int    iterations)
{
  if (!validate_non_empty(z)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;

  Array talus = Array(z.shape);
  Array zf = Array(z.shape);

  Array z_bckp = Array();
  if ((p_erosion_map != nullptr) | (p_deposition_map != nullptr)) z_bckp = z;

  for (int it = 0; it < iterations; it++)
  {
    if (ir > 0)
    {
      // pre-filter the highmap to avoid "kinky" artifacts
      zf = z;
      smooth_cpulse(zf, ir);
      gradient_talus(zf, talus);
    }
    else
      gradient_talus(z, talus);

    for (int j = 0; j < z.shape.y; j++)
      for (int i = 0; i < z.shape.x; i++)
      {
        if (talus(i, j) > talus_ref)
          z(i, j) -= c_erosion * (talus(i, j) / talus_ref - 1.f);
        else
          z(i, j) += c_deposition * (1.f - talus(i, j) / talus_ref);
      }

    // make sure bedrock is not eroded
    if (p_bedrock) z = maximum(z, *p_bedrock);
  }

  // splatmaps
  if (p_erosion_map)
  {
    *p_erosion_map = z_bckp - z;
    clamp_min(*p_erosion_map, 0.f);
  }

  if (p_deposition_map)
  {
    *p_deposition_map = z - z_bckp;
    clamp_min(*p_deposition_map, 0.f);
  }
}

void hydraulic_algebric(Array &z,
                        Array *p_mask,
                        float  talus_ref,
                        int    ir,
                        Array *p_bedrock,
                        Array *p_erosion_map,
                        Array *p_deposition_map,
                        float  c_erosion,
                        float  c_deposition,
                        int    iterations)
{
  if (!validate_non_empty(z)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;

  if (!p_mask)
    hydraulic_algebric(z,
                       talus_ref,
                       ir,
                       p_bedrock,
                       p_erosion_map,
                       p_deposition_map,
                       c_erosion,
                       c_deposition,
                       iterations);
  else
  {
    Array z_f = z;
    hydraulic_algebric(z_f,
                       talus_ref,
                       ir,
                       p_bedrock,
                       p_erosion_map,
                       p_deposition_map,
                       c_erosion,
                       c_deposition,
                       iterations);
    z = lerp(z, z_f, *(p_mask));
  }
}

} // namespace hmap
