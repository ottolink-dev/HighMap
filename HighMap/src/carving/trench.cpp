/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/carving.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/math/profiles.hpp"
#include "highmap/vectors.hpp"

namespace hmap
{

namespace
{

struct Segment
{
  float  x0, y0, v0;
  float  x1, y1, v1;
  float  dx, dy;
  float  inv_len_sq;
  float  s0, s1;
  size_t k0, k1;
};

} // namespace

void trench(Array                       &z,
            const Path                  &path,
            float                        width,
            bool                         enable_width_depth_scaling,
            bool                         enable_width_distance_scaling,
            bool                         enable_width_curvature_scaling,
            float                        curvature_radius_min,
            float                        curv_width_ratio_min,
            float                        curv_width_ratio_max,
            RadialProfile                radial_profile,
            float                        radial_profile_parameter,
            ElevationLongitudinalProfile longitudinal_profile,
            float                        elevation_shift,
            float                        shift_ramp_start_ratio,
            float                        shift_ramp_end_ratio,
            float                        min_slope,
            size_t                       k_neighbors,
            const Array                 *p_noise_r,
            Array                       *p_bending_mask,
            glm::vec4                    bbox)
{
  (void)k_neighbors;
  if (!validate_non_empty(z)) return;
  if (p_noise_r && !validate_same_shape(z, *p_noise_r)) return;
  if (!validate_non_empty(path, "Path")) return;

  const glm::ivec2 &shape = z.shape;

  // path working copy
  Path   path_copy = path;
  size_t npts = path_copy.size();

  // for width with distance scaling
  const std::vector<float> arc_length = path_copy.get_arc_length();

  // --- Adjust longitudinal elevation

  // overall shift
  for (size_t k = 0; k < npts; ++k)
  {
    // shape factor
    float t = arc_length[k];

    // start and end ramps
    if (t < shift_ramp_start_ratio)
      t /= shift_ramp_start_ratio;
    else if (t > 1.f - shift_ramp_end_ratio)
      t = (1.f - t) / shift_ramp_end_ratio;
    else
      t = 1.f;

    t = smoothstep3(t);

    // apply
    path_copy[k].v += t * elevation_shift;
  }

  // monotonicity
  for (size_t k = 1; k < npts; ++k)
  {
    // minimum elevation delta (use when monotonicity is enforced)
    float dz_min = min_slope *
                   glm::length(glm::vec2(path_copy[k].x - path_copy[k - 1].x,
                                         path_copy[k].y - path_copy[k - 1].y));

    switch (longitudinal_profile)
    {
    case ElevationLongitudinalProfile::ELP_FLAT:
      path_copy[k].v = path_copy[0].v;
      break;

    case ElevationLongitudinalProfile::ELP_DECREASING:
      path_copy[k].v = std::min(path_copy[k].v, path_copy[k - 1].v - dz_min);
      break;

    case ElevationLongitudinalProfile::ELP_INCREASING:
      path_copy[k].v = std::max(path_copy[k].v, path_copy[k - 1].v + dz_min);
      break;

    case ElevationLongitudinalProfile::ELP_UNCHANGED:
      // nothing here
      break;
    }
  }

  // --- Curvature scaling

  Path               path_curv = path_copy;
  std::vector<float> curvature = path_curv.get_curvature();
  std::vector<float> curv_radius;
  std::vector<float> curv_shape_factor;
  curv_radius.reserve(npts);
  curv_shape_factor.reserve(npts);

  if (!curvature.empty())
  {
    float curvature_max = 1.f / curvature_radius_min;
    float cmin = *std::min_element(curvature.begin(), curvature.end());
    float cmax = *std::max_element(curvature.begin(), curvature.end());
    float cscale = std::min(std::max(std::abs(cmax), std::abs(cmin)),
                            curvature_max);

    for (auto &v : curvature)
    {
      // clamp then normalize and smooth transition
      v = std::clamp(v, -curvature_max, curvature_max);
      v = std::copysign(1.f, v) * smoothstep3(std::abs(v / cscale));

      if (v != 0.f)
      {
        float r = 1.f / v;
        r = std::copysign(1.f, r) * std::min(1.f, std::abs(r));
        curv_radius.push_back(r);
      }
      else
      {
        curv_radius.push_back(0.f);
      }
    }

    remap(curv_radius, -1.f, 1.f);

    // longitudinal scaling => zero the scaling at the sign change to
    // avoid numerical artefacts
    std::vector<size_t> sign_changes = find_sign_changes(curvature);

    if (!sign_changes.empty())
    {
      std::vector<float> arc_gap(npts);
      float              arc_gap_max = 0.f;

      size_t sign_count = 0;
      size_t k0 = 0;
      size_t k1 = sign_changes[0];

      for (size_t k = 0; k < npts; ++k)
      {
        // move to next interval
        while (k >= k1 && sign_count < sign_changes.size())
        {
          k0 = k1;

          sign_count++;

          if (sign_count < sign_changes.size())
            k1 = sign_changes[sign_count];
          else
            k1 = npts; // end (not npts - 1)
        }

        // arc length for this section
        arc_gap[k] = arc_length[k1] - arc_length[k0];
        arc_gap_max = std::max(arc_gap_max, arc_gap[k]);

        // safe interval length
        float t = 0.f;
        if (k1 > k0)
          t = float(k - k0) / float(k1 - k0);
        else
          t = 1.f;

        // triangle shape and smoothing
        t = 1.f - std::abs(2.f * t - 1.f);
        t = smoothstep3(t);

        curv_shape_factor.push_back(t);
      }

      for (size_t k = 0; k < npts; ++k)
        curv_shape_factor[k] *= arc_gap[k] / arc_gap_max;

      curv_shape_factor = moving_average(curv_shape_factor, 1);
    }
  }

  // radial profile
  auto profile_fct = get_radial_profile_function(radial_profile,
                                                 radial_profile_parameter);

  // --- Calculate maximum influence width

  float max_effective_width = width;
  if (p_noise_r)
    max_effective_width = std::max(0.f, width * (1.f + p_noise_r->max()));
  if (enable_width_curvature_scaling)
    max_effective_width *= std::max(1.f, curv_width_ratio_max);

  if (max_effective_width <= 0.f) return;

  // interpolation base grid
  std::vector<float> xg, yg;
  grid_xy_vector(xg, yg, shape, bbox, /* endpoint */ false);

  float dx = (bbox.y - bbox.x) / float(shape.x);
  float dy = (bbox.w - bbox.z) / float(shape.y);

  Array zp = z;
  Array blending_mask(shape);

  // --- Continuous segment capsule rasterization with tile culling

  if (npts == 1)
  {
    // single point influence
    float px = path_copy[0].x;
    float py = path_copy[0].y;
    float pv = path_copy[0].v;

    int imin = std::clamp(
        static_cast<int>(std::floor((px - max_effective_width - bbox.x) / dx)),
        0,
        shape.x - 1);
    int imax = std::clamp(
        static_cast<int>(std::ceil((px + max_effective_width - bbox.x) / dx)),
        0,
        shape.x - 1);
    int jmin = std::clamp(
        static_cast<int>(std::floor((py - max_effective_width - bbox.z) / dy)),
        0,
        shape.y - 1);
    int jmax = std::clamp(
        static_cast<int>(std::ceil((py + max_effective_width - bbox.z) / dy)),
        0,
        shape.y - 1);

#pragma omp parallel for schedule(dynamic, 16)
    for (int j = jmin; j <= jmax; ++j)
    {
      float yi = yg[j];
      for (int i = imin; i <= imax; ++i)
      {
        float xi = xg[i];
        float d_sq = (xi - px) * (xi - px) + (yi - py) * (yi - py);
        float dr = p_noise_r ? (*p_noise_r)(i, j) : 0.f;
        float effective_width = std::max(0.f, width * (1.f + dr));

        if (enable_width_depth_scaling)
        {
          float dz = std::abs((z(i, j) - pv) / elevation_shift);
          effective_width *= std::clamp(dz, 0.f, 1.f);
        }

        if (effective_width <= 0.f || d_sq > effective_width * effective_width)
          continue;

        float r = std::sqrt(d_sq) / effective_width;
        if (r >= 0.f && r <= 1.f)
        {
          float t = profile_fct(r);
          zp(i, j) = lerp(pv, z(i, j), t);
          blending_mask(i, j) = 1.f - t;
        }
      }
    }
  }
  else
  {
    // build segments
    size_t               num_segments = npts - 1;
    std::vector<Segment> segments(num_segments);

    for (size_t k = 0; k < num_segments; ++k)
    {
      const auto &p0 = path_copy[k];
      const auto &p1 = path_copy[k + 1];

      Segment seg;
      seg.x0 = p0.x;
      seg.y0 = p0.y;
      seg.v0 = p0.v;
      seg.x1 = p1.x;
      seg.y1 = p1.y;
      seg.v1 = p1.v;
      seg.dx = p1.x - p0.x;
      seg.dy = p1.y - p0.y;
      float len_sq = seg.dx * seg.dx + seg.dy * seg.dy;
      seg.inv_len_sq = (len_sq > 1e-12f) ? (1.f / len_sq) : 0.f;
      seg.s0 = arc_length[k];
      seg.s1 = arc_length[k + 1];
      seg.k0 = k;
      seg.k1 = k + 1;
      segments[k] = seg;
    }

    // tile spatial partitioning
    constexpr int tile_size = 32;
    int           num_tiles_x = (shape.x + tile_size - 1) / tile_size;
    int           num_tiles_y = (shape.y + tile_size - 1) / tile_size;
    int           total_tiles = num_tiles_x * num_tiles_y;

    std::vector<std::vector<uint32_t>> tile_segments(total_tiles);

    for (size_t k = 0; k < num_segments; ++k)
    {
      const auto &seg = segments[k];

      float seg_xmin = std::min(seg.x0, seg.x1) - max_effective_width;
      float seg_xmax = std::max(seg.x0, seg.x1) + max_effective_width;
      float seg_ymin = std::min(seg.y0, seg.y1) - max_effective_width;
      float seg_ymax = std::max(seg.y0, seg.y1) + max_effective_width;

      int imin = std::clamp(
          static_cast<int>(std::floor((seg_xmin - bbox.x) / dx)),
          0,
          shape.x - 1);
      int imax = std::clamp(
          static_cast<int>(std::ceil((seg_xmax - bbox.x) / dx)),
          0,
          shape.x - 1);
      int jmin = std::clamp(
          static_cast<int>(std::floor((seg_ymin - bbox.z) / dy)),
          0,
          shape.y - 1);
      int jmax = std::clamp(
          static_cast<int>(std::ceil((seg_ymax - bbox.z) / dy)),
          0,
          shape.y - 1);

      int t_imin = imin / tile_size;
      int t_imax = imax / tile_size;
      int t_jmin = jmin / tile_size;
      int t_jmax = jmax / tile_size;

      for (int ty = t_jmin; ty <= t_jmax; ++ty)
      {
        for (int tx = t_imin; tx <= t_imax; ++tx)
        {
          tile_segments[ty * num_tiles_x + tx].push_back(
              static_cast<uint32_t>(k));
        }
      }
    }

#pragma omp parallel for schedule(dynamic, 8)
    for (int tile_idx = 0; tile_idx < total_tiles; ++tile_idx)
    {
      const auto &segs = tile_segments[tile_idx];
      if (segs.empty()) continue;

      int tx = tile_idx % num_tiles_x;
      int ty = tile_idx / num_tiles_x;

      int i_start = tx * tile_size;
      int i_end = std::min(i_start + tile_size, shape.x);
      int j_start = ty * tile_size;
      int j_end = std::min(j_start + tile_size, shape.y);

      for (int j = j_start; j < j_end; ++j)
      {
        float yi = yg[j];

        for (int i = i_start; i < i_end; ++i)
        {
          float xi = xg[i];

          // find closest segment
          float    min_dist_sq = std::numeric_limits<float>::max();
          float    best_t = 0.f;
          uint32_t best_seg_idx = segs[0];

          for (uint32_t seg_idx : segs)
          {
            const auto &seg = segments[seg_idx];
            float       qx = xi - seg.x0;
            float       qy = yi - seg.y0;
            float t = std::clamp((qx * seg.dx + qy * seg.dy) * seg.inv_len_sq,
                                 0.f,
                                 1.f);
            float proj_x = seg.x0 + t * seg.dx;
            float proj_y = seg.y0 + t * seg.dy;
            float d_sq = (xi - proj_x) * (xi - proj_x) +
                         (yi - proj_y) * (yi - proj_y);

            if (d_sq < min_dist_sq)
            {
              min_dist_sq = d_sq;
              best_t = t;
              best_seg_idx = seg_idx;
            }
          }

          const auto &seg = segments[best_seg_idx];
          size_t      k0 = (best_t <= 0.5f) ? seg.k0 : seg.k1;
          float       zref = (1.f - best_t) * seg.v0 + best_t * seg.v1;
          float       arc = (1.f - best_t) * seg.s0 + best_t * seg.s1;

          float dr = p_noise_r ? (*p_noise_r)(i, j) : 0.f;
          float effective_width = std::max(0.f, width * (1.f + dr));

          if (enable_width_distance_scaling)
            effective_width *= smoothstep3(arc);

          if (enable_width_depth_scaling)
          {
            float dz = std::abs((z(i, j) - zref) / elevation_shift);
            effective_width *= std::clamp(dz, 0.f, 1.f);
          }

          if (enable_width_curvature_scaling && !curvature.empty())
          {
            size_t km = (k0 > 0) ? (k0 - 1) : 0;
            size_t kp = std::min(k0 + 1, npts - 1);

            float s = -classify_point(path_copy[km],
                                      path_copy[k0],
                                      path_copy[kp],
                                      Point(xi, yi));

            float t_long = arc;
            t_long = t_long * (1.f - t_long) * 4.f;

            float camp = std::abs(curv_radius[k0]) * t_long;

            effective_width *= lerp(
                1.f,
                std::max(curv_width_ratio_min,
                         1.f + (curv_width_ratio_max - 1.f) * s * camp),
                curv_shape_factor[k0]);
          }

          if (effective_width <= 0.f ||
              min_dist_sq > effective_width * effective_width)
          {
            continue;
          }

          float r = std::sqrt(min_dist_sq) / effective_width;
          if (r >= 0.f && r <= 1.f)
          {
            float t = profile_fct(r);
            zp(i, j) = lerp(zref, z(i, j), t);
            blending_mask(i, j) = 1.f - t;
          }
        }
      }
    }
  }

  // --- Outputs

  if (p_bending_mask) *p_bending_mask = std::move(blending_mask);

  z = std::move(zp);
}

} // namespace hmap
