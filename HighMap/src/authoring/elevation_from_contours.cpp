/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap/authoring.hpp"
#include "highmap/geometry/cell_path.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"

namespace hmap
{

namespace
{

// Contour vertices converted to (float) pixel coordinates, using the same
// mapping as Cloud::to_array.
std::vector<glm::vec2> contour_to_pixel_coords(const Path &path,
                                               glm::ivec2  shape,
                                               glm::vec4   bbox)
{
  const float ai = (shape.x - 1) / (bbox.y - bbox.x);
  const float aj = (shape.y - 1) / (bbox.w - bbox.z);

  std::vector<glm::vec2> pts;
  pts.reserve(path.points.size());
  for (const Point &p : path.points)
    pts.emplace_back(ai * (p.x - bbox.x), aj * (p.y - bbox.z));
  return pts;
}

// Pixels of the closed polygon outline (Bresenham, closing segment included).
std::vector<glm::ivec2> rasterize_outline(const std::vector<glm::vec2> &pts,
                                          glm::ivec2                    shape)
{
  std::vector<glm::ivec2> cells;
  const size_t            n = pts.size();

  for (size_t k = 0; k < n; ++k)
  {
    const glm::vec2 &a = pts[k];
    const glm::vec2 &b = pts[(k + 1) % n];
    add_line_bresenham(cells,
                       {(int)std::round(a.x), (int)std::round(a.y)},
                       {(int)std::round(b.x), (int)std::round(b.y)});
  }

  // keep in-grid cells only
  std::vector<glm::ivec2> out;
  out.reserve(cells.size());
  for (const glm::ivec2 &c : cells)
    if (c.x >= 0 && c.x < shape.x && c.y >= 0 && c.y < shape.y)
      out.push_back(c);
  return out;
}

} // namespace

Array elevation_from_contours(glm::ivec2                shape,
                              const std::vector<Path>  &contours,
                              const std::vector<float> &elevations,
                              const Array              *p_probability,
                              float                     randomness,
                              uint                      seed,
                              float                     peak_ratio,
                              float                     outside_ratio,
                              glm::vec4                 bbox)
{
  (void)randomness;
  (void)seed;
  (void)peak_ratio;
  (void)outside_ratio;

  // --- validation
  if (!validate_shape(shape)) return Array();

  if (contours.empty())
  {
    log::error("elevation_from_contours: at least one contour is required");
    return Array();
  }

  if (contours.size() != elevations.size())
  {
    log::error("elevation_from_contours: contours ({}) and elevations ({}) "
               "sizes differ",
               contours.size(),
               elevations.size());
    return Array();
  }

  for (size_t k = 0; k < contours.size(); ++k)
    if (contours[k].points.size() < 3)
    {
      log::error("elevation_from_contours: contour {} has fewer than 3 points",
                 k);
      return Array();
    }

  if (p_probability && !validate_same_shape(shape, *p_probability))
    return Array();

  // --- rasterize contours: outline pixels are fixed to the contour elevation
  Array z(shape, 0.f);

  for (size_t k = 0; k < contours.size(); ++k)
  {
    std::vector<glm::vec2> pts = contour_to_pixel_coords(contours[k],
                                                         shape,
                                                         bbox);
    for (const glm::ivec2 &c : rasterize_outline(pts, shape))
      z(c.x, c.y) = elevations[k];
  }

  return z;
}

} // namespace hmap
