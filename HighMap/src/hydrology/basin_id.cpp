/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

Array basin_id(const Array &z, bool fast_approx, bool remove_lakes)
{
  if (!validate_non_empty(z)) return Array();

  Array ids(z.shape);

  // --- stream structure

  FlowGraph fg;
  fg.resize(z.shape.x, z.shape.y);
  fg.flood(z, nullptr, fast_approx);
  fg.receivers(1.f, 0);

  // --- basin Ids

  int id_count = 0;

  for (int o = 0; o < fg.n; ++o)
  {
    const int i = fg.order[o];
    const int r = fg.rcv[i];

    if (r == i)
      ids.vector[i] = static_cast<float>(id_count++);
    else
      ids.vector[i] = ids.vector[r];
  }

  return ids;
}

} // namespace hmap
