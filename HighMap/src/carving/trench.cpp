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
  float    x0, y0, v0;
  float    x1, y1, v1;
  float    dx, dy;
  float    inv_len_sq;
  float    s0, s1;
  size_t   k0, k1;
  uint32_t path_index;
};

} // namespace

void trench(Array                       &z,
            const std::vector<Path>     &paths,
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
  if (paths.empty()) return;

  const glm::ivec2 &shape = z.shape;

  // --- Preprocess paths

  struct PreprocessedPath
  {
    Path               path;
    std::vector<float> arc_length;
    std::vector<float> curvature;
    std::vector<float> curv_radius;
    std::vector<float> curv_shape_factor;
  };

  std::vector<PreprocessedPath> prep_paths;
  prep_paths.reserve(paths.size());

  for (const auto &path : paths)
  {
    if (path.empty()) continue;

    PreprocessedPath pp;
    pp.path = path;
    size_t npts = pp.path.size();
    pp.arc_length = pp.path.get_arc_length();

    // overall shift
    for (size_t k = 0; k < npts; ++k)
    {
      float t = pp.arc_length[k];

      if (t < shift_ramp_start_ratio)
        t /= shift_ramp_start_ratio;
      else if (t > 1.f - shift_ramp_end_ratio)
        t = (1.f - t) / shift_ramp_end_ratio;
      else
        t = 1.f;

      t = smoothstep3(t);
      pp.path[k].v += t * elevation_shift;
    }

    // monotonicity
    for (size_t k = 1; k < npts; ++k)
    {
      float dz_min = min_slope *
                     glm::length(glm::vec2(pp.path[k].x - pp.path[k - 1].x,
                                           pp.path[k].y - pp.path[k - 1].y));

      switch (longitudinal_profile)
      {
      case ElevationLongitudinalProfile::ELP_FLAT:
        pp.path[k].v = pp.path[0].v;
        break;

      case ElevationLongitudinalProfile::ELP_DECREASING:
        pp.path[k].v = std::min(pp.path[k].v, pp.path[k - 1].v - dz_min);
        break;

      case ElevationLongitudinalProfile::ELP_INCREASING:
        pp.path[k].v = std::max(pp.path[k].v, pp.path[k - 1].v + dz_min);
        break;

      case ElevationLongitudinalProfile::ELP_UNCHANGED: break;
      }
    }

    // curvature scaling
    if (enable_width_curvature_scaling)
    {
      Path path_curv = pp.path;
      pp.curvature = path_curv.get_curvature();
      pp.curv_radius.reserve(npts);
      pp.curv_shape_factor.reserve(npts);

      if (!pp.curvature.empty())
      {
        float curvature_max = 1.f / curvature_radius_min;
        float cmin = *std::min_element(pp.curvature.begin(),
                                       pp.curvature.end());
        float cmax = *std::max_element(pp.curvature.begin(),
                                       pp.curvature.end());
        float cscale = std::min(std::max(std::abs(cmax), std::abs(cmin)),
                                curvature_max);

        for (auto &v : pp.curvature)
        {
          v = std::clamp(v, -curvature_max, curvature_max);
          v = std::copysign(1.f, v) * smoothstep3(std::abs(v / cscale));

          if (v != 0.f)
          {
            float r = 1.f / v;
            r = std::copysign(1.f, r) * std::min(1.f, std::abs(r));
            pp.curv_radius.push_back(r);
          }
          else
          {
            pp.curv_radius.push_back(0.f);
          }
        }

        remap(pp.curv_radius, -1.f, 1.f);

        std::vector<size_t> sign_changes = find_sign_changes(pp.curvature);

        if (!sign_changes.empty())
        {
          std::vector<float> arc_gap(npts);
          float              arc_gap_max = 0.f;

          size_t sign_count = 0;
          size_t k0 = 0;
          size_t k1 = sign_changes[0];

          for (size_t k = 0; k < npts; ++k)
          {
            while (k >= k1 && sign_count < sign_changes.size())
            {
              k0 = k1;
              sign_count++;

              if (sign_count < sign_changes.size())
                k1 = sign_changes[sign_count];
              else
                k1 = npts;
            }

            arc_gap[k] = pp.arc_length[k1] - pp.arc_length[k0];
            arc_gap_max = std::max(arc_gap_max, arc_gap[k]);

            float t = 0.f;
            if (k1 > k0)
              t = float(k - k0) / float(k1 - k0);
            else
              t = 1.f;

            t = 1.f - std::abs(2.f * t - 1.f);
            t = smoothstep3(t);

            pp.curv_shape_factor.push_back(t);
          }

          for (size_t k = 0; k < npts; ++k)
            pp.curv_shape_factor[k] *= arc_gap[k] / arc_gap_max;

          pp.curv_shape_factor = moving_average(pp.curv_shape_factor, 1);
        }
      }
    }

    prep_paths.push_back(std::move(pp));
  }

  if (prep_paths.empty()) return;

  // --- Calculate maximum influence width

  float max_effective_width = width;
  if (p_noise_r)
    max_effective_width = std::max(0.f, width * (1.f + p_noise_r->max()));
  if (enable_width_curvature_scaling)
    max_effective_width *= std::max(1.f, curv_width_ratio_max);

  if (max_effective_width <= 0.f) return;

  // radial profile
  auto profile_fct = get_radial_profile_function(radial_profile,
                                                 radial_profile_parameter);

  // interpolation base grid
  std::vector<float> xg, yg;
  grid_xy_vector(xg, yg, shape, bbox, /* endpoint */ false);

  float dx = (bbox.y - bbox.x) / float(shape.x);
  float dy = (bbox.w - bbox.z) / float(shape.y);

  Array blending_mask;
  if (p_bending_mask) blending_mask = Array(shape);

  // build segments across all paths
  std::vector<Segment> segments;
  for (size_t p_idx = 0; p_idx < prep_paths.size(); ++p_idx)
  {
    const auto &pp = prep_paths[p_idx];
    size_t      npts = pp.path.size();

    if (npts == 1)
    {
      const auto &p0 = pp.path[0];
      Segment     seg;
      seg.x0 = p0.x;
      seg.y0 = p0.y;
      seg.v0 = p0.v;
      seg.x1 = p0.x;
      seg.y1 = p0.y;
      seg.v1 = p0.v;
      seg.dx = 0.f;
      seg.dy = 0.f;
      seg.inv_len_sq = 0.f;
      seg.s0 = 0.f;
      seg.s1 = 0.f;
      seg.k0 = 0;
      seg.k1 = 0;
      seg.path_index = static_cast<uint32_t>(p_idx);
      segments.push_back(seg);
    }
    else
    {
      for (size_t k = 0; k + 1 < npts; ++k)
      {
        const auto &p0 = pp.path[k];
        const auto &p1 = pp.path[k + 1];

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
        seg.s0 = pp.arc_length[k];
        seg.s1 = pp.arc_length[k + 1];
        seg.k0 = k;
        seg.k1 = k + 1;
        seg.path_index = static_cast<uint32_t>(p_idx);
        segments.push_back(seg);
      }
    }
  }

  if (segments.empty()) return;

  // tile spatial partitioning
  constexpr int tile_size = 32;
  int           num_tiles_x = (shape.x + tile_size - 1) / tile_size;
  int           num_tiles_y = (shape.y + tile_size - 1) / tile_size;
  int           total_tiles = num_tiles_x * num_tiles_y;

  std::vector<std::vector<uint32_t>> tile_segments(total_tiles);

  for (size_t k = 0; k < segments.size(); ++k)
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
    int imax = std::clamp(static_cast<int>(std::ceil((seg_xmax - bbox.x) / dx)),
                          0,
                          shape.x - 1);
    int jmin = std::clamp(
        static_cast<int>(std::floor((seg_ymin - bbox.z) / dy)),
        0,
        shape.y - 1);
    int jmax = std::clamp(static_cast<int>(std::ceil((seg_ymax - bbox.z) / dy)),
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
        const auto &pp = prep_paths[seg.path_index];
        size_t      k0 = (best_t <= 0.5f) ? seg.k0 : seg.k1;
        float       zref = (1.f - best_t) * seg.v0 + best_t * seg.v1;
        float       arc = (1.f - best_t) * seg.s0 + best_t * seg.s1;

        float dr = p_noise_r ? (*p_noise_r)(i, j) : 0.f;
        float effective_width = std::max(0.f, width * (1.f + dr));

        if (enable_width_distance_scaling) effective_width *= smoothstep3(arc);

        if (enable_width_depth_scaling)
        {
          float dz = std::abs((z(i, j) - zref) / elevation_shift);
          effective_width *= std::clamp(dz, 0.f, 1.f);
        }

        if (enable_width_curvature_scaling && !pp.curvature.empty())
        {
          size_t km = (k0 > 0) ? (k0 - 1) : 0;
          size_t kp = std::min(k0 + 1, pp.path.size() - 1);

          float s = -classify_point(pp.path[km],
                                    pp.path[k0],
                                    pp.path[kp],
                                    Point(xi, yi));

          float t_long = arc;
          t_long = t_long * (1.f - t_long) * 4.f;

          float camp = std::abs(pp.curv_radius[k0]) * t_long;

          effective_width *= lerp(
              1.f,
              std::max(curv_width_ratio_min,
                       1.f + (curv_width_ratio_max - 1.f) * s * camp),
              pp.curv_shape_factor[k0]);
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
          z(i, j) = lerp(zref, z(i, j), t);
          if (p_bending_mask) blending_mask(i, j) = 1.f - t;
        }
      }
    }
  }

  // --- Outputs

  if (p_bending_mask) *p_bending_mask = std::move(blending_mask);
}

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
  if (!validate_non_empty(path, "Path")) return;

  std::vector<Path> paths = {path};
  trench(z,
         paths,
         width,
         enable_width_depth_scaling,
         enable_width_distance_scaling,
         enable_width_curvature_scaling,
         curvature_radius_min,
         curv_width_ratio_min,
         curv_width_ratio_max,
         radial_profile,
         radial_profile_parameter,
         longitudinal_profile,
         elevation_shift,
         shift_ramp_start_ratio,
         shift_ramp_end_ratio,
         min_slope,
         k_neighbors,
         p_noise_r,
         p_bending_mask,
         bbox);
}

} // namespace hmap
