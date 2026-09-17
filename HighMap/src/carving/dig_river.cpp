/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/morphology.hpp"

namespace hmap
{

void dig_river(Array                   &z,
               const std::vector<Path> &path_list,
               float                    riverbank_talus,
               int                      river_width,
               int                      merging_width,
               float                    depth,
               float                    riverbed_talus,
               float                    noise_ratio,
               std::uint32_t            seed,
               Array                   *p_mask)
{
  if (!validate_non_empty(z)) return;
  if (!validate_non_empty(path_list, "Path list")) return;

  // generate mask where the river path lies and dig rivers
  Array       mask(z.shape);
  hmap::Array z_carved(z.shape);

  for (auto path : path_list)
  {
    Path path_copy = path;
    path_copy.set_values(1.f);
    glm::vec4 bbox(0.f, 1.f, 0.f, 1.f);
    path_copy.to_array(mask, bbox);

    // expand the path
    path_copy = path;
    path_copy.enforce_monotonic_values();

    for (auto &p : path_copy)
      p.v -= depth;

    // add downstream slope
    if (riverbed_talus > 0.f)
    {
      for (size_t k = 0; k < path_copy.size() - 1; k++)
        path_copy[k + 1].v = std::min(path_copy[k + 1].v,
                                      path_copy[k].v - riverbed_talus);
    }

    path_copy.to_array(z_carved, bbox);
  }

  if (river_width)
  {
    z_carved = dilation(z_carved, river_width);
    mask = dilation(mask, river_width);
  }

  for (int j = 0; j < z.shape.y; ++j)
    for (int i = 0; i < z.shape.x; ++i)
      z_carved(i, j) = mask(i, j) ? z_carved(i, j) : z(i, j);

  expand_talus(z_carved, mask, riverbank_talus, seed, noise_ratio);
  laplace(z_carved);

  // use a distance transform to define a merging mask between the
  // input heightmap "z" and the "z_carved"
  Array dist = distance_transform_approx(mask);
  float sigma2 = std::pow((float)(merging_width + river_width), 2.f);
  dist = exp(-0.5f * dist * dist / sigma2);
  laplace(dist);

  // lerp based on distance
  z = lerp(z, z_carved, dist);

  // return mask if requested
  if (p_mask) *p_mask = dist;
}

void dig_river(Array        &z,
               const Path   &path,
               float         riverbank_talus,
               int           river_width,
               int           merging_width,
               float         depth,
               float         riverbed_talus,
               float         noise_ratio,
               std::uint32_t seed,
               Array        *p_mask)
{
  const std::vector<Path> path_list = {path};

  dig_river(z,
            path_list,
            riverbank_talus,
            river_width,
            merging_width,
            depth,
            riverbed_talus,
            noise_ratio,
            seed,
            p_mask);
}

} // namespace hmap
