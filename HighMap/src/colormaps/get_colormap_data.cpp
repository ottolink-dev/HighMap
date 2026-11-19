/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <vector>

#include "highmap/colormaps.hpp"
#include "highmap/internal/validation.hpp"

namespace hmap
{

std::vector<glm::vec3> helper_tov3(
    const std::vector<std::vector<float>> &colors)
{
  std::vector<glm::vec3> out;
  out.reserve(colors.size());
  for (const auto &v : colors)
    out.push_back({v[0], v[1], v[2]});
  return out;
}

std::vector<glm::vec3> get_colormap_data(int cmap)
{
  if (!validate_cmap(cmap))
  {
    return helper_tov3(CMAP_WHITE_UNIFORM);
  }

  switch (cmap)
  {
  case Cmap::BONE: return helper_tov3(CMAP_BONE);
  case Cmap::GRAY: return helper_tov3(CMAP_GRAY);
  case Cmap::HOT: return helper_tov3(CMAP_HOT);
  case Cmap::INFERNO: return helper_tov3(CMAP_INFERNO);
  case Cmap::JET: return helper_tov3(CMAP_JET);
  case Cmap::MAGMA: return helper_tov3(CMAP_MAGMA);
  case Cmap::NIPY_SPECTRAL: return helper_tov3(CMAP_NIPY_SPECTRAL);
  case Cmap::RAINBOW_DEM: return helper_tov3(CMAP_RAINBOW_DEM);
  case Cmap::SEISMIC: return helper_tov3(CMAP_SEISMIC);
  case Cmap::TERRAIN: return helper_tov3(CMAP_TERRAIN);
  case Cmap::TURBO: return helper_tov3(CMAP_TURBO);
  case Cmap::VIRIDIS: return helper_tov3(CMAP_VIRIDIS);
  case Cmap::WHITE_UNIFORM: return helper_tov3(CMAP_WHITE_UNIFORM);
  default: return helper_tov3(CMAP_WHITE_UNIFORM);
  }
}

} // namespace hmap
