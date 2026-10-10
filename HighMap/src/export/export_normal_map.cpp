/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <string>

#include "highmap/array.hpp"
#include "highmap/export/image.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/texture.hpp"

namespace hmap
{

void export_normal_map_png(const std::string &fname,
                           const Array       &array,
                           int                depth)
{
  if (!validate_non_empty(array)) return;

  Texture nmap = normal_map(array);
  nmap.to_png(fname, depth);
}

} // namespace hmap
