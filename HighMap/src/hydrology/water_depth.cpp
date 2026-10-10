/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cstdint>
#include <iterator>
#include <limits>
#include <queue>
#include <string>
#include <utility>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/curvature.hpp"
#include "highmap/filters.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate2d.hpp"
#include "highmap/math/array.hpp"
#include "highmap/morphology.hpp"
#include "highmap/range.hpp"

namespace hmap
{

Array merge_water_depths(const Array &depth1,
                         const Array &depth2,
                         float        k_smooth)
{
  if (!validate_non_empty(depth1)) return Array();
  if (!validate_same_shape(depth1, depth2)) return Array();

  Array water_depth(depth1.shape);

  if (k_smooth == 0.f)
    water_depth = maximum(depth1, depth2);
  else
    water_depth = maximum_smooth(depth1, depth2, k_smooth) - k_smooth / 6.f;

  return water_depth;
}

void water_depth_dry_out(Array       &water_depth,
                         float        dry_out_ratio,
                         const Array *p_mask,
                         float        depth_max)
{
  if (!validate_non_empty(water_depth)) return;
  if (p_mask && !validate_same_shape(water_depth, *p_mask)) return;

  if (depth_max == std::numeric_limits<float>::max())
    depth_max = water_depth.max();

  if (p_mask)
  {
    for (int j = 0; j < water_depth.shape.y; j++)
      for (int i = 0; i < water_depth.shape.x; i++)
      {
        water_depth(i, j) -= dry_out_ratio * depth_max * (*p_mask)(i, j);
        water_depth(i, j) = std::max(0.f, water_depth(i, j));
      }
  }
  else
  {
    for (int j = 0; j < water_depth.shape.y; j++)
      for (int i = 0; i < water_depth.shape.x; i++)
      {
        water_depth(i, j) -= dry_out_ratio * depth_max;
        water_depth(i, j) = std::max(0.f, water_depth(i, j));
      }
  }
}

Array water_depth_from_mask(const Array &z,
                            const Array &mask,
                            float        mask_threshold,
                            int          iterations_max,
                            float        tolerance)
{
  if (!validate_non_empty(z)) return Array();
  if (!validate_same_shape(z, mask)) return Array();

  Array water_depth(z.shape);

  // transform to binary 0|1 mask
  Array mask_t = mask;
  make_binary(mask_t, mask_threshold);
  mask_t = 1.f - mask_t; // fixed values

  water_depth = harmonic_interpolation(z, mask_t, iterations_max, tolerance) -
                z;

  return water_depth;
}

Array water_depth_increase(const Array &water_depth,
                           const Array &z,
                           float        additional_depth)
{
  if (!validate_non_empty(water_depth)) return Array();
  if (!validate_same_shape(water_depth, z)) return Array();

  const int ni = water_depth.shape.x;
  const int nj = water_depth.shape.y;
  Array     water_depth_extended(water_depth.shape);

  constexpr int di[8] = {-1, 1, 0, 0, -1, -1, 1, 1};
  constexpr int dj[8] = {0, 0, -1, 1, -1, 1, -1, 1};

  auto in_bounds = [&](int i, int j) noexcept
  { return i >= 0 && i < ni && j >= 0 && j < nj; };

  const int            ncells = ni * nj;
  std::vector<int>     queue;
  std::vector<uint8_t> in_queue(ncells, 0);
  queue.reserve(ncells / 4);

  // --- Seed: copy water cells, enqueue border cells

  for (int j = 0; j < nj; ++j)
    for (int i = 0; i < ni; ++i)
    {
      if (water_depth(i, j) <= 0.f) continue;

      water_depth_extended(i, j) = water_depth(i, j) + additional_depth;

      const int flat = j * ni + i;
      if (in_queue[flat]) continue;

      for (int k = 0; k < 8; ++k)
      {
        int qi = i + di[k], qj = j + dj[k];
        if (in_bounds(qi, qj) && water_depth(qi, qj) <= 0.f)
        {
          queue.push_back(flat);
          in_queue[flat] = 1;
          break;
        }
      }
    }

  // --- Upward flood propagation (toward higher terrain)

  for (int head = 0; head < (int)queue.size(); ++head)
  {
    const int flat_p = queue[head];
    in_queue[flat_p] = 0;

    const int   pi = flat_p % ni;
    const int   pj = flat_p / ni;
    const float base_depth = water_depth_extended(pi, pj);
    const float base_z = z(pi, pj);

    for (int k = 0; k < 8; ++k)
    {
      const int ni_ = pi + di[k], nj_ = pj + dj[k];
      if (!in_bounds(ni_, nj_)) continue;

      const float dz = z(ni_, nj_) - base_z;
      if (dz <= 0.f) continue;

      const float propagated = base_depth - dz;
      if (propagated <= 0.f) continue;

      const int flat_n = nj_ * ni + ni_;
      if (propagated > water_depth_extended(ni_, nj_))
      {
        water_depth_extended(ni_, nj_) = propagated;
        if (!in_queue[flat_n])
        {
          queue.push_back(flat_n);
          in_queue[flat_n] = 1;
        }
      }
    }
  }

  // Hole-filling: downward propagation into depressions

  queue.clear();
  std::fill(in_queue.begin(), in_queue.end(), 0);

  for (int j = 0; j < nj; ++j)
    for (int i = 0; i < ni; ++i)
    {
      if (water_depth_extended(i, j) <= 0.f) continue;

      const int flat = j * ni + i;
      if (in_queue[flat]) continue;

      for (int k = 0; k < 8; ++k)
      {
        int qi = i + di[k], qj = j + dj[k];
        if (in_bounds(qi, qj) && water_depth_extended(qi, qj) <= 0.f)
        {
          queue.push_back(flat);
          in_queue[flat] = 1;
          break;
        }
      }
    }

  for (int head = 0; head < (int)queue.size(); ++head)
  {
    const int flat_p = queue[head];
    in_queue[flat_p] = 0;

    const int   pi = flat_p % ni;
    const int   pj = flat_p / ni;
    const float base_depth = water_depth_extended(pi, pj);
    const float base_z = z(pi, pj);

    for (int k = 0; k < 8; ++k)
    {
      const int ni_ = pi + di[k], nj_ = pj + dj[k];
      if (!in_bounds(ni_, nj_)) continue;

      const float dz = z(ni_, nj_) - base_z;
      if (dz >= 0.f) continue;

      const float propagated = base_depth + dz;
      if (propagated <= 0.f) continue;

      const int flat_n = nj_ * ni + ni_;
      if (propagated > water_depth_extended(ni_, nj_))
      {
        water_depth_extended(ni_, nj_) = propagated;
        if (!in_queue[flat_n])
        {
          queue.push_back(flat_n);
          in_queue[flat_n] = 1;
        }
      }
    }
  }

  return water_depth_extended;
}

Array water_depth_increase_with_flooding(const Array &water_depth,
                                         const Array &z,
                                         float        additional_depth)
{
  if (!validate_non_empty(water_depth)) return Array();
  if (!validate_same_shape(water_depth, z)) return Array();

  // WSE == Water Surface Elevation (z + water_depth)

  const int ni = water_depth.shape.x;
  const int nj = water_depth.shape.y;
  Array     water_depth_extended(water_depth.shape); // zero-initialized

  constexpr int DI[8] = {-1, 1, 0, 0, -1, -1, 1, 1};
  constexpr int DJ[8] = {0, 0, -1, 1, -1, 1, -1, 1};

  auto in_bounds = [&](int i, int j) noexcept
  { return (unsigned)i < (unsigned)ni && (unsigned)j < (unsigned)nj; };

  // max-heap: {water_surface_elevation, flat_index}
  using Entry = std::pair<float, int>;
  std::priority_queue<Entry> pq;
  std::vector<bool>          visited(ni * nj, false);

  // --- Seed: every wet cell defines a WSE

  for (int j = 0; j < nj; ++j)
    for (int i = 0; i < ni; ++i)
    {
      if (water_depth(i, j) <= 0.f) continue;
      const float wse = z(i, j) + water_depth(i, j) + additional_depth;
      pq.push({wse, j * ni + i});
    }

  // --- Single propagation pass (highest WSE settled first)

  while (!pq.empty())
  {
    auto [wse, flat] = pq.top();
    pq.pop();

    if (visited[flat]) continue; // already settled at a higher WSE
    visited[flat] = true;

    const int pi = flat % ni, pj = flat / ni;
    water_depth_extended(pi, pj) = wse - z(pi, pj);

    for (int k = 0; k < 8; ++k)
    {
      const int i2 = pi + DI[k], j2 = pj + DJ[k];
      if (!in_bounds(i2, j2)) continue;
      const int flat_n = j2 * ni + i2;
      if (visited[flat_n]) continue;
      if (z(i2, j2) < wse)      // terrain is below water surface
        pq.push({wse, flat_n}); // neighbor inherits the same WSE
    }
  }

  return water_depth_extended;
}

Array water_frontier_curvature(const Array &water_depth,
                               int          prefilter_ir,
                               bool         extend_values_from_interface)
{
  if (!validate_non_empty(water_depth)) return Array();

  const glm::ivec2 &shape = water_depth.shape;

  Mat<glm::ivec2> closest_in(shape);
  Mat<glm::ivec2> closest_out(shape);

  Array dist_in = distance_transform_with_closest(is_zero(water_depth),
                                                  closest_in);
  Array dist_out = distance_transform_with_closest(water_depth, closest_out);
  Array phi = dist_in - dist_out;
  phi = level_set_curvature(phi, prefilter_ir);

  if (extend_values_from_interface)
  {
    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        if (water_depth(i, j) > 0.f)
          phi(i, j) = phi(closest_in(i, j));
        else
          phi(i, j) = phi(closest_out(i, j));
      }
  }

  return phi;
}

Array water_mask(const Array &water_depth)
{
  if (!validate_non_empty(water_depth)) return Array();

  Array mask = water_depth;
  make_binary(mask);
  return mask;
}

Array water_mask(const Array &water_depth,
                 const Array &z,
                 float        additional_depth)
{
  if (!validate_non_empty(water_depth)) return Array();
  if (!validate_same_shape(water_depth, z)) return Array();

  Array mask(water_depth.shape);
  Array water_depth_extended = water_depth_increase(water_depth,
                                                    z,
                                                    additional_depth);

  mask = water_depth_extended - water_depth;
  mask /= additional_depth;
  mask = smoothstep3(mask);

  return mask;
}

} // namespace hmap

namespace hmap::gpu
{

void water_depth_filter(Array       &depth,
                        const Array &z,
                        int          ir,
                        const Array *p_water_mask,
                        bool         smooth_contour,
                        float        transition_ratio)
{
  if (!validate_non_empty(depth)) return;
  if (!validate_same_shape(depth, z)) return;
  if (p_water_mask && !validate_same_shape(depth, *p_water_mask)) return;

  const glm::ivec2 shape = depth.shape;
  Array            zt = z + depth;

  // if no water mask is provided to describe where there should be
  // water, just use the water depth
  if (!p_water_mask) p_water_mask = &depth;

  Array smooth_mask;
  if (smooth_contour)
  {
    smooth_mask = gpu::contour_smoothing(*p_water_mask, ir, transition_ratio);
    p_water_mask = &smooth_mask;
  }

  auto run = clwrapper::Run("water_depth_filter");

  run.bind_imagef("depth", depth.vector, shape.x, shape.y);
  run.bind_imagef("water_mask", p_water_mask->vector, shape.x, shape.y);
  run.bind_imagef("zt", zt.vector, shape.x, shape.y);
  run.bind_imagef("zt_out", zt.vector, shape.x, shape.y, true);

  run.bind_arguments(shape.x, shape.y, ir);

  run.execute({shape.x, shape.y});
  run.read_imagef("zt_out");

  // retrieve depth
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      if ((*p_water_mask)(i, j) != 0.f)
        depth(i, j) = std::max(0.f, zt(i, j) - z(i, j));
}

Array water_depth_from_mask(const Array &z,
                            const Array &mask,
                            float        mask_threshold,
                            int          iterations_max,
                            float        tolerance)
{
  if (!validate_non_empty(z)) return Array();
  if (!validate_same_shape(z, mask)) return Array();

  Array water_depth(z.shape);

  // transform to binary 0|1 mask
  Array mask_t = mask;
  make_binary(mask_t, mask_threshold);
  mask_t = 1.f - mask_t; // fixed values

  water_depth = gpu::harmonic_interpolation(z,
                                            mask_t,
                                            iterations_max,
                                            tolerance) -
                z;

  return water_depth;
}

Array water_frontier_curvature(const Array &water_depth,
                               int          prefilter_ir,
                               bool         extend_values_from_interface)
{
  if (!validate_non_empty(water_depth)) return Array();

  const glm::ivec2 &shape = water_depth.shape;

  Mat<glm::ivec2> closest_in(shape);
  Mat<glm::ivec2> closest_out(shape);

  Array dist_in = distance_transform_with_closest(is_zero(water_depth),
                                                  closest_in);
  Array dist_out = distance_transform_with_closest(water_depth, closest_out);
  Array phi = dist_in - dist_out;
  phi = gpu::level_set_curvature(phi, prefilter_ir);

  if (extend_values_from_interface)
  {
    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        if (water_depth(i, j) > 0.f)
          phi(i, j) = phi(closest_in(i, j));
        else
          phi(i, j) = phi(closest_out(i, j));
      }
  }

  return phi;
}

} // namespace hmap::gpu
