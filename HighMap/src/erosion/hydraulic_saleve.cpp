/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <future>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/erosion/hydraulic_erosion.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/hydrology/drainage_basin.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate2d.hpp"
#include "highmap/math/array.hpp"
#include "highmap/primitives/functions.hpp"
#include "highmap/range.hpp"
#include "highmap/spatial/uv_tile.hpp"
#include "highmap/terrain_tri_mesh.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

void hydraulic_saleve(TerrainTriMesh           &mesh,
                      const std::vector<float> &erodibility,
                      const std::vector<float> &max_slope,
                      float                     m_exp,
                      float                     uplift_rate,
                      float                     tolerance,
                      int                       max_iterations,
                      float                     noise_strength,
                      std::uint32_t             seed,
                      bool                      enable_post_slope_limiter,
                      float                     post_slope_limit,
                      bool                      enable_post_smoothing)
{
  auto db = DrainageBasin(mesh.get_points());
  db.set_outlets(find_border_sinks(db.get_mesh()));

  for (int it = 0; it < max_iterations; ++it)
  {
    std::vector acc(db.size(), 0.f);
    db.update_stream_tree(seed, noise_strength);

    auto area = db.get_mesh().get_vertex_areas(true);
    db.accumulate_area_by_outlet(area, acc);

    auto response_times = db.compute_response_times(acc, erodibility, m_exp);

    float diff = db.update_elevations(response_times, uplift_rate, max_slope);

    if (diff < tolerance) break;
  }

  // post-treatments
  if (enable_post_slope_limiter)
  {
    glm::vec2 zr = db.get_mesh().get_range_z();
    float     zptp = zr.y - zr.x;
    db.get_mesh().slope_limiter(post_slope_limit * zptp, 100);
  }

  if (enable_post_smoothing) db.get_mesh().relax_xyz(0.2f, 5);

  // override input with eroded field
  mesh = TerrainTriMesh(db.get_mesh());
}

Array hydraulic_saleve(const Array          &z,
                       std::uint32_t         seed,
                       size_t                control_points_count,
                       float                 m_exp,
                       float                 uplift_rate,
                       float                 tolerance,
                       int                   max_iterations,
                       float                 smin,
                       float                 smax,
                       float                 strength,
                       bool                  scale_erodibility_with_z,
                       float                 erodibility_distrib_exp,
                       float                 noise_strength,
                       bool                  enable_post_slope_limiter,
                       float                 post_slope_limit,
                       bool                  enable_post_smoothing,
                       InterpolationMethod2D interpolation_method,
                       const Array          *p_noise_x,
                       const Array          *p_noise_y)
{
  if (!validate_non_empty(z)) return Array();
  if (p_noise_x && !validate_same_shape(z, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(z, *p_noise_y)) return Array();

  const glm::ivec2 shape = z.shape;
  const float      zmin = z.min();
  const float      zmax = z.max();

  // safeguard
  smin = std::min(smin, smax);

  // --- generate triangle mesh

  TerrainTriMesh mesh = generate_terrain_tri_mesh_from_heightmap_random(
      z,
      control_points_count,
      seed);

  // --- spatial parameters

  // erodibility align with elevation
  std::vector<float> erodibility;

  if (scale_erodibility_with_z && (zmin != zmax))
  {
    erodibility.reserve(mesh.size());
    for (const auto &p : mesh.get_points())
    {
      float v = (p.z - zmin) / (zmax - zmin);
      v = std::pow(1.f - v, erodibility_distrib_exp);
      erodibility.push_back(v);
    }
  }
  else
  {
    erodibility = std::vector<float>(control_points_count, 1.f);
  }

  // slope varying with distance to the boundary
  std::vector<float> max_slope = cubic_pulse(mesh);

  for (auto &v : max_slope)
    v = (smax - smin) * v + smin;

  // --- erode the triangle mesh

  hydraulic_saleve(mesh,
                   erodibility,
                   max_slope,
                   m_exp,
                   uplift_rate,
                   tolerance,
                   max_iterations,
                   noise_strength,
                   seed,
                   enable_post_slope_limiter,
                   post_slope_limit,
                   enable_post_smoothing);

  // --- interpolate back to an heightmap

  // make sure the input noise displacement do not modify the convex
  // hull limits to avoid issues with the nautral neighbor
  // interpolation
  Array dx;
  Array dy;

  if (p_noise_x)
  {
    dx = (*p_noise_x) * biquad_pulse_x(shape);
    p_noise_x = &dx;
  }

  if (p_noise_y)
  {
    dy = (*p_noise_y) * biquad_pulse_y(shape);
    p_noise_y = &dy;
  }

  // interpolate
  std::vector<float> xc, yc, zc;
  for (const auto &p : mesh.get_points())
  {
    xc.push_back(p.x);
    yc.push_back(p.y);
    zc.push_back(p.z);
  }

  Array ze = interpolate2d(shape,
                           xc,
                           yc,
                           zc,
                           interpolation_method,
                           p_noise_x,
                           p_noise_y);

  remap(ze, zmin, zmax);

  return lerp(z, ze, strength);
}

Array hydraulic_saleve(const Array          &z,
                       const Array          *p_mask,
                       std::uint32_t         seed,
                       size_t                control_points_count,
                       float                 m_exp,
                       float                 uplift_rate,
                       float                 tolerance,
                       int                   max_iterations,
                       float                 smin,
                       float                 smax,
                       float                 strength,
                       bool                  scale_erodibility_with_z,
                       float                 erodibility_distrib_exp,
                       float                 noise_strength,
                       bool                  enable_post_slope_limiter,
                       float                 post_slope_limit,
                       bool                  enable_post_smoothing,
                       InterpolationMethod2D interpolation_method,
                       const Array          *p_noise_x,
                       const Array          *p_noise_y)
{
  if (!validate_non_empty(z)) return Array();
  if (p_mask && !validate_same_shape(z, *p_mask)) return Array();
  if (p_noise_x && !validate_same_shape(z, *p_noise_x)) return Array();
  if (p_noise_y && !validate_same_shape(z, *p_noise_y)) return Array();

  Array ze = hydraulic_saleve(z,
                              seed,
                              control_points_count,
                              m_exp,
                              uplift_rate,
                              tolerance,
                              max_iterations,
                              smin,
                              smax,
                              strength,
                              scale_erodibility_with_z,
                              erodibility_distrib_exp,
                              noise_strength,
                              enable_post_slope_limiter,
                              post_slope_limit,
                              enable_post_smoothing,
                              interpolation_method,
                              p_noise_x,
                              p_noise_y);

  if (!p_mask)
    return ze;
  else
    return lerp(z, ze, *(p_mask));
}

} // namespace hmap

namespace hmap::va
{

VirtualArray hydraulic_saleve(const ComputeMode    &cm,
                              const VirtualArray   &z,
                              std::uint32_t         seed,
                              size_t                control_points_count,
                              float                 m_exp,
                              float                 uplift_rate,
                              float                 tolerance,
                              int                   max_iterations,
                              float                 smin,
                              float                 smax,
                              float                 strength,
                              bool                  scale_erodibility_with_z,
                              float                 erodibility_distrib_exp,
                              float                 noise_strength,
                              bool                  enable_post_slope_limiter,
                              float                 post_slope_limit,
                              bool                  enable_post_smoothing,
                              InterpolationMethod2D interpolation_method,
                              const VirtualArray   *p_noise_x,
                              const VirtualArray   *p_noise_y,
                              const VirtualArray   *p_mask)
{
  const float zmin = z.min(cm);
  const float zmax = z.max(cm);

  // safeguard
  smin = std::min(smin, smax);

  // --- generate triangle mesh

  TerrainTriMesh mesh;

  hmap::for_each_tile(
      {&z},
      {},
      [&](std::vector<const hmap::Array *> p_arrays_in,
          std::vector<hmap::Array *>,
          const hmap::UvTile &region,
          ICommunicator      &comm)
      {
        auto [pa_z] = unpack<1>(p_arrays_in);

        TerrainTriMesh tile_mesh =
            hmap::generate_terrain_tri_mesh_from_heightmap_random(
                *pa_z,
                control_points_count,
                seed,
                region,
                &comm);

        if (comm.rank() == 0) mesh = std::move(tile_mesh);
      },
      cm);

  // --- spatial parameters

  // erodibility align with elevation
  std::vector<float> erodibility;

  if (scale_erodibility_with_z && (zmin != zmax))
  {
    erodibility.reserve(mesh.size());
    for (const auto &p : mesh.get_points())
    {
      float v = (p.z - zmin) / (zmax - zmin);
      v = std::pow(1.f - v, erodibility_distrib_exp);
      erodibility.push_back(v);
    }
  }
  else
  {
    erodibility = std::vector<float>(control_points_count, 1.f);
  }

  // slope varying with distance to the boundary
  std::vector<float> max_slope = cubic_pulse(mesh);

  for (auto &v : max_slope)
    v = (smax - smin) * v + smin;

  // --- erode the triangle mesh

  hydraulic_saleve(mesh,
                   erodibility,
                   max_slope,
                   m_exp,
                   uplift_rate,
                   tolerance,
                   max_iterations,
                   noise_strength,
                   seed,
                   enable_post_slope_limiter,
                   post_slope_limit,
                   enable_post_smoothing);

  // --- interpolate back to an heightmap

  // make sure the input noise displacement do not modify the convex
  // hull limits to avoid issues with the nautral neighbor
  // interpolation

  std::vector<float> xc, yc, zc;
  for (const auto &p : mesh.get_points())
  {
    xc.push_back(p.x);
    yc.push_back(p.y);
    zc.push_back(p.z);
  }

  VirtualArray ze;
  ze.copy_from(z, cm, /* copy_src_data */ false);

  hmap::for_each_tile(
      {p_noise_x, p_noise_y},
      {&ze},
      [&](std::vector<const hmap::Array *> p_arrays_in,
          std::vector<hmap::Array *>       p_arrays_out,
          const hmap::UvTile              &region)
      {
        auto [pa_noise_x, pa_noise_y] = unpack<2>(p_arrays_in);
        auto [pa_ze] = unpack<1>(p_arrays_out);

        Array dx;
        Array dy;

        if (pa_noise_x)
        {
          dx = (*pa_noise_x) *
               biquad_pulse_x(region.shape, region.uv_rect.to_vec4());
          pa_noise_x = &dx;
        }

        if (pa_noise_y)
        {
          dy = (*pa_noise_y) *
               biquad_pulse_y(region.shape, region.uv_rect.to_vec4());
          pa_noise_y = &dy;
        }

        *pa_ze = interpolate2d(region.shape,
                               xc,
                               yc,
                               zc,
                               interpolation_method,
                               pa_noise_x,
                               pa_noise_y,
                               region.uv_rect.to_vec4());
      },
      cm);

  ze.remap(zmin, zmax, cm);

  hmap::for_each_tile(
      {&z, p_mask},
      {&ze},
      [&](std::vector<const hmap::Array *> p_arrays_in,
          std::vector<hmap::Array *>       p_arrays_out,
          const hmap::UvTile &)
      {
        auto [pa_z, pa_mask] = unpack<2>(p_arrays_in);
        auto [pa_ze] = unpack<1>(p_arrays_out);

        if (pa_mask)
          *pa_ze = lerp(*pa_z, *pa_ze, strength * (*pa_mask));
        else
          *pa_ze = lerp(*pa_z, *pa_ze, strength);
      },
      cm);

  return ze;
}

} // namespace hmap::va
