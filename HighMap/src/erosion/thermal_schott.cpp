/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <cstdint>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/erosion.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/primitives/functions.hpp"

namespace hmap
{

void thermal_schott(Array       &z,
                    const Array &talus,
                    int          iterations,
                    float        intensity)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, talus)) return;

  // https://www.shadertoy.com/view/XX2XWD

  std::vector<int>    di = HMAP_DI;
  std::vector<int>    dj = HMAP_DJ;
  std::vector<float>  c = HMAP_CD;
  const std::uint32_t nb = di.size();

  for (int it = 0; it < iterations; it++)
  {
    Array z_new = z;

    for (int j = 1; j < z.shape.y - 1; j++)
      for (int i = 1; i < z.shape.x - 1; i++)
      {
        int up = 0;
        int down = 0;

        // check neighbors
        for (std::uint32_t k = 0; k < nb; k++)
        {
          int p = i + di[k];
          int q = j + dj[k];

          float slope = (z(i, j) - z(p, q)) / c[k];

          if (slope > talus(i, j))
            down++;
          else if (slope < -talus(i, j))
            up++;
        }

        z_new(i, j) += intensity * talus(i, j) * (float)(up - down);
      }

    z = z_new;
    extrapolate_borders(z);
  }
}

void thermal_schott(Array       &z,
                    const Array &talus,
                    Array       *p_mask,
                    int          iterations,
                    float        intensity)
{
  if (!validate_non_empty(z) || !validate_same_shape(z, talus)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;

  if (!p_mask)
    thermal_schott(z, talus, iterations, intensity);
  else
  {
    Array z_f = z;
    thermal_schott(z_f, talus, iterations, intensity);
    z = lerp(z, z_f, *(p_mask));
  }
}

void thermal_schott(Array      &z,
                    const float talus,
                    int         iterations,
                    float       intensity)
{
  if (!validate_non_empty(z)) return;

  Array talus_map = constant(z.shape, talus);
  thermal_schott(z, talus_map, iterations, intensity);
}

void thermal_schott(Array      &z,
                    const float talus,
                    Array      *p_mask,
                    int         iterations,
                    float       intensity)
{
  if (!validate_non_empty(z)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;

  if (!p_mask)
    thermal_schott(z, talus, iterations, intensity);
  else
  {
    Array z_f = z;
    thermal_schott(z_f, talus, iterations, intensity);
    z = lerp(z, z_f, *(p_mask));
  }
}

} // namespace hmap
