/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <array>
#include <cstdint>

#include "highmap/array.hpp"
#include "highmap/erosion/deprecated.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/primitives/functions.hpp"

namespace hmap
{

//----------------------------------------------------------------------
// Main operator
//----------------------------------------------------------------------

void hydraulic_musgrave(Array &z,
                        Array &moisture_map,
                        int    iterations,
                        float  c_capacity,
                        float  c_erosion,
                        float  c_deposition,
                        float  water_level,
                        float  evap_rate)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, moisture_map)) return;

  Array s = constant(z.shape);          // sediment level
  Array w = water_level * moisture_map; // backup initial moisture map

  auto                neighbors = neighborhood::MOORE_8;
  const std::uint32_t nb = neighbors.size();

  for (int it = 0; it < iterations; it++)
  {
    w = (1 - evap_rate) * w + evap_rate * moisture_map * water_level;

    // modify neighbor search at each iterations to limit numerical
    // artifacts
    std::rotate(neighbors.begin(), neighbors.begin() + 1, neighbors.end());

    for (int j = 1; j < z.shape.y - 1; j++)
      for (int i = 1; i < z.shape.x - 1; i++)
        for (std::uint32_t k = 0; k < nb; k++) // loop over 1st neighbors
        {
          int   p = i + neighbors[k].offset.x;
          int   q = j + neighbors[k].offset.y;
          float dw = std::min(w(i, j),
                              (w(i, j) + z(i, j) - w(p, q) - z(p, q)) *
                                  neighbors[k].inv_distance);

          if (dw <= 0.f)
          {
            // static deposition when water is stagnant
            float ds = c_deposition * s(i, j);
            s(i, j) -= ds;
            z(i, j) += ds;
          }
          else
          {
            // water transfer
            w(i, j) -= dw;
            w(p, q) += dw;

            // sediment carrying capacity
            float cs = c_capacity * dw;

            if (s(i, j) > cs)
            {
              // deposition
              float ds = c_deposition * (s(i, j) - cs);
              s(p, q) += cs;
              z(i, j) += ds;
              s(i, j) = (1.f - c_deposition) * (s(i, j) - cs);
            }
            else
            {
              // erosion
              float ds = c_erosion * (cs - s(i, j));
              s(p, q) += s(i, j) + ds;
              z(i, j) -= ds;
              s(i, j) = 0.f;
            }
          }
        }
  }
}

void hydraulic_musgrave(Array &z,
                        int    iterations,
                        float  c_capacity,
                        float  c_erosion,
                        float  c_deposition,
                        float  water_level,
                        float  evap_rate)
{
  if (!validate_non_empty(z)) return;

  Array moisture_map = constant(z.shape, 1.f);
  hydraulic_musgrave(z,
                     moisture_map,
                     iterations,
                     c_capacity,
                     c_erosion,
                     c_deposition,
                     water_level,
                     evap_rate);
}

} // namespace hmap
