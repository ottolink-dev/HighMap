/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <algorithm>
#include <cstdint>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/erosion.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/primitives/functions.hpp"
#include "highmap/range.hpp"

namespace hmap
{

//----------------------------------------------------------------------
// Main operator(s)
//----------------------------------------------------------------------

void hydraulic_benes(Array &z,
                     int    iterations,
                     Array *p_bedrock,
                     Array *p_moisture_map,
                     Array *p_erosion_map,
                     Array *p_deposition_map,
                     float  c_capacity,
                     float  c_erosion,
                     float  c_deposition,
                     float  water_level,
                     float  evap_rate,
                     float  rain_rate)
{
  if (!validate_non_empty(z)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;
  if (p_moisture_map && !validate_same_shape(z, *p_moisture_map)) return;

  auto                neighbors = neighborhood::MOORE_8;
  const std::uint32_t nb = neighbors.size();

  // keep a backup of the input if the erosion / deposition maps need
  // to be computed
  Array z_bckp = Array();
  if ((p_erosion_map != nullptr) | (p_deposition_map != nullptr)) z_bckp = z;

  Array w = water_level * constant(z.shape, 1.f);

  if (p_moisture_map) w *= (*p_moisture_map); // water map

  Array w_init = w;                 // backup water initial map
  Array s = constant(z.shape, 0.f); // sediment map (disolved in the water)

  float wmin = w.max() * 0.01f; // minimum water limit

  Array vel = Array(z.shape);

  // main loop
  for (int it = 0; it < iterations; it++)
  {
    // modify neighbor search at each iterations to limit numerical
    // artifacts
    std::rotate(neighbors.begin(), neighbors.begin() + 1, neighbors.end());

    w = (1.f - rain_rate) * w + rain_rate * w_init;

    // --- water flow dynamic and sediment transport

    for (int j = 1; j < z.shape.y - 1; j++)
      for (int i = 1; i < z.shape.x - 1; i++)
      {
        float              dsum = 0.f;
        float              zsavg = 0.f;
        int                navg = 0;
        std::vector<float> dz(nb);

        for (std::uint32_t k = 0; k < nb; k++)
        {
          int ia = i + neighbors[k].offset.x;
          int ja = j + neighbors[k].offset.y;
          dz[k] = z(i, j) + w(i, j) - z(ia, ja) - w(ia, ja);

          if (dz[k] > 0.f)
          {
            dsum += dz[k];
            zsavg += z(ia, ja) + w(ia, ja);
            navg++;
          }
          else
            dz[k] = 0.f;
        }

        if (dsum > 0.f and w(i, j) > wmin)
        {
          zsavg /= (float)navg;

          float dw_tot = std::min(w(i, j), z(i, j) + w(i, j) - zsavg);
          float ds_tot;

          if (w(i, j) > 0.f)
            ds_tot = s(i, j) * dw_tot / w(i, j);
          else
            ds_tot = 0.f;

          w(i, j) -= dw_tot;
          s(i, j) -= ds_tot;

          vel(i, j) = dw_tot;

          for (std::uint32_t k = 0; k < nb; k++)
          {
            if (dz[k] > 0.f)
            {
              int   ia = i + neighbors[k].offset.x;
              int   ja = j + neighbors[k].offset.y;
              float r = dz[k] / dsum;

              w(ia, ja) += dw_tot * r;
              s(ia, ja) += ds_tot * r;
            }
          }
        }
      }

    // --- erosion and deposition
    for (int i = 1; i < z.shape.x - 1; i++)
    {
      for (int j = 1; j < z.shape.y - 1; j++)
      {
        float              dsum = 0.f;
        std::vector<float> dz(nb);

        for (std::uint32_t k = 0; k < nb; k++)
        {
          int ia = i + neighbors[k].offset.x;
          int ja = j + neighbors[k].offset.y;
          dz[k] = z(i, j) + w(i, j) - z(ia, ja) - w(ia, ja); // no water

          if (dz[k] > 0.f)
            dsum += dz[k];
          else
            dz[k] = 0.f;
        }

        if (dsum > 0.f)
        {
          float ds_tot = c_capacity * w(i, j) * vel(i, j) - s(i, j);
          float amount;

          if (ds_tot > 0.f)
            amount = c_erosion * ds_tot;
          else
            amount = c_deposition * ds_tot;

          z(i, j) -= amount;
          s(i, j) += amount;

          for (std::uint32_t k = 0; k < nb; k++)
          {
            if (dz[k] > 0.f)
            {
              int   ia = i + neighbors[k].offset.x;
              int   ja = j + neighbors[k].offset.y;
              float r = dz[k] / dsum;

              z(ia, ja) -= amount * r;
              s(ia, ja) += amount * r;
            }
          }
        }

      } // j
    } // i

    w = (1.f - evap_rate) * w;
    chop(w, wmin);

    extrapolate_borders(z);
    fill_borders(w);
    fill_borders(s);

    laplace(w);
    laplace(s);

    // make sure bedrock is not eroded
    if (p_bedrock)
      for (int i = 0; i < z.shape.x; i++)
        for (int j = 0; j < z.shape.y; j++)
          z(i, j) = std::max(z(i, j), (*p_bedrock)(i, j));

  } // it

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

void hydraulic_benes(Array &z,
                     Array *p_mask,
                     int    iterations,
                     Array *p_bedrock,
                     Array *p_moisture_map,
                     Array *p_erosion_map,
                     Array *p_deposition_map,
                     float  c_capacity,
                     float  c_erosion,
                     float  c_deposition,
                     float  water_level,
                     float  evap_rate,
                     float  rain_rate)
{
  if (!validate_non_empty(z)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;
  if (p_moisture_map && !validate_same_shape(z, *p_moisture_map)) return;

  if (!p_mask)
    hydraulic_benes(z,
                    iterations,
                    p_bedrock,
                    p_moisture_map,
                    p_erosion_map,
                    p_deposition_map,
                    c_capacity,
                    c_erosion,
                    c_deposition,
                    water_level,
                    evap_rate,
                    rain_rate);
  else
  {
    Array z_f = z;
    hydraulic_benes(z_f,
                    iterations,
                    p_bedrock,
                    p_moisture_map,
                    p_erosion_map,
                    p_deposition_map,
                    c_capacity,
                    c_erosion,
                    c_deposition,
                    water_level,
                    evap_rate,
                    rain_rate);
    z = lerp(z, z_f, *(p_mask));
  }
}

} // namespace hmap
