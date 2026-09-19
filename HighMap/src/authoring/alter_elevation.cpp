/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/kernels.hpp"

namespace hmap
{

void alter_elevation(Array       &array,
                     const Cloud &cloud,
                     int          ir,
                     float        footprint_ratio,
                     glm::vec2    shift,
                     glm::vec2    scale)
{
  if (!validate_non_empty(array)) return;

  Array amp = Array(array.shape);
  for (const auto &p : cloud)
  {
    int ic = (int)((p.x - shift.x) / scale.x * array.shape.x);
    int jc = (int)((p.y - shift.y) / scale.y * array.shape.y);

    // kernel size
    int   nk = (int)((2 * ir + 1) * (1.f + footprint_ratio * std::abs(p.v)));
    Array kernel = cubic_pulse(glm::ivec2(nk, nk));

    // --
    // truncate kernel to make it fit into the heightmap array
    int nk_i0 = (int)(std::floor(0.5f * kernel.shape.x)); // left
    int nk_i1 = kernel.shape.x - nk_i0;                   // right
    int nk_j0 = (int)(std::floor(0.5f * kernel.shape.y));
    int nk_j1 = kernel.shape.y - nk_j0;

    int ik0 = std::max(0, nk_i0 - ic);
    int jk0 = std::max(0, nk_j0 - jc);
    int ik1 = std::min(kernel.shape.x,
                       kernel.shape.x - (ic + nk_i1 - amp.shape.x));
    int jk1 = std::min(kernel.shape.y,
                       kernel.shape.y - (jc + nk_j1 - amp.shape.y));

    // where it goes in the array
    int i0 = std::max(ic - nk_i0, 0);
    int j0 = std::max(jc - nk_j0, 0);
    // int i1 = std::min(ic + nk_i1, amp.shape.x);
    // int j1 = std::min(jc + nk_j1, amp.shape.y);

    float sign = (array(ic, jc) > 0.f) - (array(ic, jc) < 0.f);

    for (int j = jk0; j < jk1; j++)
      for (int i = ik0; i < ik1; i++)
        amp(i - ik0 + i0, j - jk0 + j0) += p.v * kernel(i, j) *
                                           array(i - ik0 + i0, j - jk0 + j0) *
                                           sign;
  }

  array += amp;
}

} // namespace hmap
