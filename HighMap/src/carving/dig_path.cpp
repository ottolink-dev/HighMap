/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cstddef>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/math/array.hpp"
#include "highmap/morphology.hpp"

namespace hmap
{

void dig_path(Array    &z,
              Path     &path,
              int       width,
              int       decay,
              int       flattening_radius,
              bool      force_downhill,
              glm::vec4 bbox,
              float     depth)
{
  if (!validate_non_empty(z)) return;
  if (!validate_non_empty(path, "Path")) return;

  Array mask = Array(z.shape);
  Path  path_copy = path;
  Array zf;

  if (force_downhill)
  {
    // make sure the path is monotically decreasing
    path_copy.set_values_from_array(z, bbox);

    for (size_t k = 1; k < path_copy.size(); k++)
      if (path_copy[k].v > path_copy[k - 1].v)
        path_copy[k].v = path_copy[k - 1].v;

    path_copy.to_array(mask, bbox);
    zf = local_max(mask, 3 * (width + decay));
    smooth_cpulse(zf, width + decay);

    // regenerate the mask
    path_copy.set_values(1.f);
    path_copy.to_array(mask, bbox);
  }
  else
  {
    // make sure values at the path points are non-zero before creating
    // the mask
    Path path_copy = path;
    path_copy.set_values(1.f);

    path_copy.to_array(mask, bbox);
    zf = local_mean(z, flattening_radius);
  }

  mask = local_max(mask, width);
  mask = distance_transform_approx(mask);
  mask = exp(-mask * mask * 0.5f / ((float)(decay * decay)));
  smooth_cpulse(mask, decay);

  zf += depth;
  z = lerp(z, zf, mask);
}

} // namespace hmap
