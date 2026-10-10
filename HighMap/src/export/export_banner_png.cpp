/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <cstdint>
#include <string>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/operator.hpp"

#include <format>

namespace hmap
{

void export_banner_png(const std::string        &fname,
                       const std::vector<Array> &arrays,
                       int                       cmap,
                       bool                      hillshading,
                       bool                      normalize_arrays)
{
  if (!validate_non_empty(arrays, "Array list")) return;
  if (!validate_cmap(cmap)) return;
  for (const auto &arr : arrays)
  {
    if (!validate_non_empty(arr)) return;
  }

  // build up big array by stacking input arrays
  if (arrays.size() > 1)
  {
    Array banner_array;

    if (normalize_arrays)
    {
      banner_array = hmap::hstack(arrays[0].remapped(), arrays[1].remapped());
      for (std::uint32_t i = 2; i < arrays.size(); i++)
        banner_array = hmap::hstack(banner_array, arrays[i].remapped());
    }
    else
    {
      banner_array = hmap::hstack(arrays[0], arrays[1]);
      for (std::uint32_t i = 2; i < arrays.size(); i++)
        banner_array = hmap::hstack(banner_array, arrays[i]);
    }

    banner_array.to_png(fname, cmap, hillshading);
  }
  else
  {
    Array banner_array = arrays[0];
    banner_array.to_png(fname, cmap, hillshading);
  }
}

} // namespace hmap
