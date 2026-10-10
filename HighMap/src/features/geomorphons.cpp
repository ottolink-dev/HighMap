/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/features.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

Array geomorphons(const Array &array, int irmin, int irmax, float epsilon)
{
  if (!validate_non_empty(array)) return Array();

  float epsilon_normed = epsilon / (float)array.shape.x;

  // output array containing the geomorphon labels
  Array gm = Array(array.shape);

  // neighborhood search
  const auto         &neighbors = neighborhood::GEOMORPHON_8;
  const std::uint32_t nb = neighbors.size();

  const int dr_max = irmax - irmin + 1;

  Mat<int> code = Mat<int>(array.shape);

  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
    {
      std::vector<int> signature(nb);

      for (size_t k = 0; k < nb; k++)
      {
        float slope_up = 0.f;
        float slope_dw = 0.f;

        for (int dr = 1; dr < dr_max; dr++)
        {
          int ip = i + neighbors[k].offset.x * dr;
          int jp = j + neighbors[k].offset.y * dr;

          if ((ip > -1) && (ip < array.shape.x) && (jp > -1) &&
              (jp < array.shape.y))
          {
            float slope = (array(ip, jp) - array(i, j)) / (float)dr;

            if (slope > 0.f)
              slope_up = std::max(slope_up, slope);
            else
              slope_dw = std::min(slope_dw, slope);
          }
        }

        if (slope_up > -slope_dw && slope_up > epsilon_normed) signature[k] = 1;
        if (slope_up < -slope_dw && slope_dw < -epsilon_normed)
          signature[k] = -1;
      }

      code(i, j) = HMAP_PACK8(signature[0],
                              signature[1],
                              signature[2],
                              signature[3],
                              signature[4],
                              signature[5],
                              signature[6],
                              signature[7]);
    }

  // determine geomorphon ids (see "Geomorphons - a new approach to
  // classification of landforms")

  // https://geomorphometry.org/wp-content/uploads/2021/07/StepinskiJasiewicz2011geomorphometry.pdf

  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
    {
      switch (code(i, j))
      {
      case HMAP_PACK8(0, 0, 0, 0, 0, 0, 0, 0): // FLAT A
                                               //
      case HMAP_PACK8(1, 0, 0, 0, 0, 0, 0, 0): // FLAT B
      case HMAP_PACK8(0, 1, 0, 0, 0, 0, 0, 0):
      case HMAP_PACK8(0, 0, 1, 0, 0, 0, 0, 0):
      case HMAP_PACK8(0, 0, 0, 1, 0, 0, 0, 0):
      case HMAP_PACK8(0, 0, 0, 0, 1, 0, 0, 0):
      case HMAP_PACK8(0, 0, 0, 0, 0, 1, 0, 0):
      case HMAP_PACK8(0, 0, 0, 0, 0, 0, 1, 0):
      case HMAP_PACK8(0, 0, 0, 0, 0, 0, 0, 1): gm(i, j) = 1.f; break;

      case HMAP_PACK8(-1, -1, -1, -1, -1, -1, -1, -1): // PEAK
        gm(i, j) = 2.f;
        break;

      case HMAP_PACK8(1, 1, 1, 1, 1, 1, 1, 1): // PIT
        gm(i, j) = 10.f;
        break;

      case HMAP_PACK8(0, -1, -1, -1, 0, -1, -1, -1): // RIDGE
      case HMAP_PACK8(-1, 0, -1, -1, -1, 0, -1, -1):
      case HMAP_PACK8(-1, -1, 0, -1, -1, -1, 0, -1):
      case HMAP_PACK8(-1, -1, -1, 0, -1, -1, -1, 0):
        //
        gm(i, j) = 3.f;
        break;

      case HMAP_PACK8(0, 1, 1, 1, 0, 1, 1, 1): // VALLEY
      case HMAP_PACK8(1, 0, 1, 1, 1, 0, 1, 1):
      case HMAP_PACK8(1, 1, 0, 1, 1, 1, 0, 1):
      case HMAP_PACK8(1, 1, 1, 0, 1, 1, 1, 0):
        //
        gm(i, j) = 9.f;
        break;

      case HMAP_PACK8(1, 1, 1, 0, -1, -1, -1, 0): // SLOPE A
      case HMAP_PACK8(0, 1, 1, 1, 0, -1, -1, -1):
      case HMAP_PACK8(-1, 0, 1, 1, 1, 0, -1, -1):
      case HMAP_PACK8(-1, -1, 0, 1, 1, 1, 0, -1):
      case HMAP_PACK8(-1, -1, -1, 0, 1, 1, 1, 0):
      case HMAP_PACK8(0, -1, -1, -1, 0, 1, 1, 1):
      case HMAP_PACK8(1, 0, -1, -1, -1, 0, 1, 1):
      case HMAP_PACK8(1, 1, 0, -1, -1, -1, 0, 1):
        //
      case HMAP_PACK8(1, 1, 1, 1, -1, -1, -1, -1): // SLOPE C
      case HMAP_PACK8(-1, 1, 1, 1, 1, -1, -1, -1):
      case HMAP_PACK8(-1, -1, 1, 1, 1, 1, -1, -1):
      case HMAP_PACK8(-1, -1, -1, 1, 1, 1, 1, -1):
      case HMAP_PACK8(-1, -1, -1, -1, 1, 1, 1, 1):
      case HMAP_PACK8(1, -1, -1, -1, -1, 1, 1, 1):
      case HMAP_PACK8(1, 1, -1, -1, -1, -1, 1, 1):
      case HMAP_PACK8(1, 1, 1, -1, -1, -1, -1, 1): gm(i, j) = 6.f; break;

      case HMAP_PACK8(1, 1, 1, 1, 1, -1, -1, -1): // SPUR
      case HMAP_PACK8(-1, 1, 1, 1, 1, 1, -1, -1):
      case HMAP_PACK8(-1, -1, 1, 1, 1, 1, 1, -1):
      case HMAP_PACK8(-1, -1, -1, 1, 1, 1, 1, 1):
      case HMAP_PACK8(1, -1, -1, -1, 1, 1, 1, 1):
      case HMAP_PACK8(1, 1, -1, -1, -1, 1, 1, 1):
      case HMAP_PACK8(1, 1, 1, -1, -1, -1, 1, 1):
      case HMAP_PACK8(1, 1, 1, 1, -1, -1, -1, 1):
        //
        gm(i, j) = 5.f;
        break;

      case HMAP_PACK8(-1, -1, -1, -1, -1, 1, 1, 1): // HOLLOW
      case HMAP_PACK8(1, -1, -1, -1, -1, -1, 1, 1):
      case HMAP_PACK8(1, 1, -1, -1, -1, -1, -1, 1):
      case HMAP_PACK8(1, 1, 1, -1, -1, -1, -1, -1):
      case HMAP_PACK8(-1, 1, 1, 1, -1, -1, -1, -1):
      case HMAP_PACK8(-1, -1, 1, 1, 1, -1, -1, -1):
      case HMAP_PACK8(-1, -1, -1, 1, 1, 1, -1, -1):
      case HMAP_PACK8(-1, -1, -1, -1, 1, 1, 1, -1):
        //
        gm(i, j) = 7.f;
        break;

      case HMAP_PACK8(1, 1, 1, 0, 0, 0, 0, 0): // FOOTSLOPE A
      case HMAP_PACK8(0, 1, 1, 1, 0, 0, 0, 0):
      case HMAP_PACK8(0, 0, 1, 1, 1, 0, 0, 0):
      case HMAP_PACK8(0, 0, 0, 1, 1, 1, 0, 0):
      case HMAP_PACK8(0, 0, 0, 0, 1, 1, 1, 0):
      case HMAP_PACK8(0, 0, 0, 0, 0, 1, 1, 1):
      case HMAP_PACK8(1, 0, 0, 0, 0, 0, 1, 1):
      case HMAP_PACK8(1, 1, 0, 0, 0, 0, 0, 1): gm(i, j) = 8.f; break;

      case HMAP_PACK8(-1, -1, -1, 0, 0, 0, 0, 0): // SHOULDER
      case HMAP_PACK8(0, -1, -1, -1, 0, 0, 0, 0):
      case HMAP_PACK8(0, 0, -1, -1, -1, 0, 0, 0):
      case HMAP_PACK8(0, 0, 0, -1, -1, -1, 0, 0):
      case HMAP_PACK8(0, 0, 0, 0, -1, -1, -1, 0):
      case HMAP_PACK8(0, 0, 0, 0, 0, -1, -1, -1):
      case HMAP_PACK8(-1, 0, 0, 0, 0, 0, -1, -1):
      case HMAP_PACK8(-1, -1, 0, 0, 0, 0, 0, -1): gm(i, j) = 4.f; break;

      default: // unknown patterns are classified as SLOPEs
        gm(i, j) = 6.f;
      }
    }

  gm.infos();

  return gm;
}

} // namespace hmap
