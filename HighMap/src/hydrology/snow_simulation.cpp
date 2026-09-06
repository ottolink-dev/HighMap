/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/erosion.hpp"
#include "highmap/gradient.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/range.hpp"

namespace hmap
{

Array snow_melting_map(const Array &z,
                       float        melt_start_elevation,
                       float        melt_end_elevation,
                       float        elevation_exp,
                       float        elevation_strength,
                       float        sun_azimuth,
                       float        sun_zenith,
                       float        aspect_strength,
                       float        slope_exp,
                       float        slope_strength)
{
  if (!validate_non_empty(z)) return Array();

  const glm::ivec2 shape = z.shape;
  Array            map(shape); // output

  const float azimuth_rad = M_PI * sun_azimuth / 180.f;
  const float zenith_rad = M_PI * sun_zenith / 180.f;
  const float cz = std::cos(zenith_rad);
  const float sz = std::sin(zenith_rad);

  // gradient features
  Array aspect = gradient_angle(z, true);
  Array dn = 0.5f * gradient_norm(z) * shape.x;

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float h = z(i, j);

      // elevation
      float m_elev = (h - melt_start_elevation) /
                     (melt_end_elevation - melt_start_elevation);
      m_elev = std::clamp(m_elev, 0.f, 1.f);
      m_elev = gain(m_elev, elevation_exp);

      // sun
      float slope = std::atan(dn(i, j));
      float m_aspect = cz * std::cos(slope) +
                       sz * std::sin(slope) *
                           std::cos(azimuth_rad - aspect(i, j));
      m_aspect = std::clamp(m_aspect, 0.f, 1.f);

      // slope
      float m_slope = std::pow(slope, slope_exp);
      m_slope = std::clamp(m_slope, 0.f, 1.f);

      // combine
      float mc = elevation_strength * (1.f - m_elev);
      mc = std::max(mc, m_aspect * aspect_strength);
      mc = std::max(mc, m_slope * slope_strength);

      map(i, j) = mc;
    }

  return map;
}

} // namespace hmap

namespace hmap::gpu
{

Array snow_simulation(const Array &z,
                      float        snow_depth,
                      const Array &fall_map,
                      const Array &melting_map,
                      const Array &talus,
                      int          iterations,
                      float        dt,
                      float        fall_iterations_ratio,
                      float        k_snow,
                      float        k_visc,
                      float        k_melt_factor,
                      float        k_depth_ratio,
                      float        k_depth_slope_ratio,
                      bool         post_filter,
                      float        thermal_talus_ratio)
{
  if (!validate_non_empty(z)) return Array();
  if (!validate_same_shape(z, fall_map)) return Array();
  if (!validate_same_shape(z, melting_map)) return Array();
  if (!validate_same_shape(z, talus)) return Array();

  const glm::ivec2 shape = z.shape;
  Array            s(shape); // output

  if (iterations <= 0) return s;

  // computed parameters
  const float k_melt = k_melt_factor / int(iterations);
  const int   iterations_fall = std::max(1,
                                       int(fall_iterations_ratio * iterations));

  // --- Device state: snow depth is a ping-pong pair (A/B) on GPU

  using clwrapper::Direction;

  auto run = clwrapper::Run("snow_simulation");

  run.bind_imagef("z", z.vector, shape.x, shape.y);
  run.bind_imagef("s_a", s.vector, shape.x, shape.y, Direction::INOUT);
  run.bind_imagef("talus", talus.vector, shape.x, shape.y);
  run.bind_imagef("melting_map", melting_map.vector, shape.x, shape.y);
  run.bind_imagef("s_b", s.vector, shape.x, shape.y, Direction::INOUT);

  run.bind_arguments(shape.x,
                     shape.y,
                     dt,
                     snow_depth,
                     k_snow,
                     k_melt,
                     k_visc,
                     k_depth_ratio,
                     k_depth_slope_ratio);

  auto run_fall = std::make_unique<clwrapper::Run>("hydraulic_vpipes_rain_pass",
                                                   run.get_queue());

  run_fall->bind_image2d("s_a", run.get_image2d("s_a"));
  run_fall->bind_imagef("fall_map", fall_map.vector, shape.x, shape.y);
  run_fall->bind_image2d("s_b", run.get_image2d("s_b"));
  run_fall->bind_arguments(shape.x,
                           shape.y,
                           snow_depth / float(iterations_fall),
                           1);

  // ping-pong handles
  const std::array<cl::Image2D, 2> img_s = {run.get_image2d("s_a").cl_image,
                                            run.get_image2d("s_b").cl_image};

  int sc = 0; // index of the current snow depth image

  for (int it = 0; it < iterations; ++it)
  {
    if (it < iterations_fall)
    {
      run_fall->set_argument(0, img_s[sc]);
      run_fall->set_argument(2, img_s[1 - sc]);
      run_fall->execute_async({shape.x, shape.y});
      sc = 1 - sc;
    }

    run.set_argument(1, img_s[sc]);
    run.set_argument(4, img_s[1 - sc]);
    run.execute_async({shape.x, shape.y});

    sc = 1 - sc;
  }

  run.finish();

  // retrieve results (both snow depth images map to the host array `s`)
  run.read_imagef(sc == 0 ? "s_a" : "s_b");

  extrapolate_borders(s);

  if (post_filter)
  {
    Array mask = s / snow_depth;
    clamp_max_smooth(mask, 1.f);
    Array z_wrk = z + s;
    gpu::thermal(z_wrk,
                 &mask,
                 thermal_talus_ratio * hmap::gradient_norm(z_wrk),
                 iterations);
    s = z_wrk - z;
  }

  clamp_min(s, 0.f);

  return s;
}

} // namespace hmap::gpu
