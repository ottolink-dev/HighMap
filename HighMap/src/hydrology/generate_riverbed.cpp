/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/filters.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/range.hpp"
#include "highmap/sdf.hpp"

namespace hmap
{

Array generate_riverbed(const Path &path,
                        glm::ivec2  shape,
                        glm::vec4   bbox,
                        bool        bezier_smoothing,
                        float       depth_start,
                        float       depth_end,
                        float       slope_start,
                        float       slope_end,
                        float       shape_exponent_start,
                        float       shape_exponent_end,
                        float       k_smoothing,
                        int         post_filter_ir,
                        Array      *p_noise_x,
                        Array      *p_noise_y,
                        Array      *p_noise_r)
{
  if (!validate_shape(shape)) return Array();
  if (!validate_min_size(path, 2, "Path points")) return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array();
  if (p_noise_r && !validate_same_shape(shape, *p_noise_r)) return Array();

  Array sdf;
  if (bezier_smoothing)
    sdf = sdf_2d_polyline_bezier(path, shape, bbox, p_noise_x, p_noise_y);
  else
    sdf = sdf_2d_polyline(path, shape, bbox, p_noise_x, p_noise_y);

  Array dz(shape);

  // normalize path nodes coord to unit square
  std::vector<float> xp = path.get_x();
  std::vector<float> yp = path.get_y();

  glm::vec2 e = {xp.back() - xp.front(), yp.back() - yp.front()};
  float     dot_ee = dot(e, e);

  // array cell coordinates
  std::vector<float> xg(shape.x);
  std::vector<float> yg(shape.y);
  grid_xy_vector(xg, yg, shape, bbox, false);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      // relative position projected onto the start/end segment (in [0, 1])
      glm::vec2 w = {xg[i] - xp.front(), yg[j] - yp.front()};
      float     t = std::clamp(dot(w, e) / dot_ee, 0.f, 1.f);

      // tweak shape with distance
      float alpha = t * shape_exponent_start + (1.f - t) * shape_exponent_end;

      // radial function
      float r = sdf(i, j);
      if (p_noise_r) r = std::max(0.f, r + (*p_noise_r)(i, j));

      r *= t * slope_start + (1.f - t) * slope_end;
      float dp = std::pow((std::pow(r, alpha) + 1.f), 1.f / alpha) - 1.f;

      dp = minimum_smooth(1.f, dp, k_smoothing);

      // adjust depth with distance
      float depth = t * depth_start + (1.f - t) * depth_end;
      dz(i, j) = (dp - 1.f) * depth;
    }

  if (post_filter_ir > 0) smooth_cpulse(dz, post_filter_ir);

  return dz;
}

} // namespace hmap
