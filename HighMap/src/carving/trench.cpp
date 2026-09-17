/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/carving.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/geometry/kd_tree.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/math/profiles.hpp"
#include "highmap/vectors.hpp"

namespace hmap
{

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
  if (!validate_non_empty(z)) return;
  if (p_noise_r && !validate_same_shape(z, *p_noise_r)) return;
  if (!validate_non_empty(path, "Path")) return;

  const glm::ivec2 &shape = z.shape;

  // path working copy
  Path   path_copy = path;
  size_t npts = path_copy.size();

  // --- Force path resolution

  // {
  //   float lx = bbox.y - bbox.x;
  //   float ly = bbox.w - bbox.z;
  //   float dmin = std::min(lx / shape.x, ly / shape.y);
  //   path_copy.resample(dmin);
  // }

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
    {
      path_copy[k].v = std::min(path_copy[k].v, path_copy[k - 1].v - dz_min);
    }
    break;

    case ElevationLongitudinalProfile::ELP_INCREASING:
      path_copy[k].v = std::max(path_copy[k].v, path_copy[k - 1].v + dz_min);
      break;

    case ElevationLongitudinalProfile::ELP_UNCHANGED:
      // nothing here
      break;
    }
  }

  // --- SDF-based transform

  Array zp = z;
  Array blending_mask(shape);

  KDTree tree(path_copy);

  // interpolation base grid
  std::vector<float> xg, yg;
  grid_xy_vector(xg, yg, shape, bbox, /* endpoint */ false);

  // for curvature scaling
  Path path_curv = path_copy;
  // path_curv.decimate_vw(40);
  // path_curv.bspline(50);

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

  // calculate bounding box of influence
  float max_effective_width = width;
  if (p_noise_r)
    max_effective_width = std::max(0.f, width * (1.f + p_noise_r->max()));
  if (enable_width_curvature_scaling)
    max_effective_width *= std::max(1.f, curv_width_ratio_max);

  float xmin_path = path_copy[0].x;
  float xmax_path = path_copy[0].x;
  float ymin_path = path_copy[0].y;
  float ymax_path = path_copy[0].y;

  for (size_t k = 1; k < npts; ++k)
  {
    xmin_path = std::min(xmin_path, path_copy[k].x);
    xmax_path = std::max(xmax_path, path_copy[k].x);
    ymin_path = std::min(ymin_path, path_copy[k].y);
    ymax_path = std::max(ymax_path, path_copy[k].y);
  }

  float xmin_roi = xmin_path - max_effective_width;
  float xmax_roi = xmax_path + max_effective_width;
  float ymin_roi = ymin_path - max_effective_width;
  float ymax_roi = ymax_path + max_effective_width;

  float dx = (bbox.y - bbox.x) / float(shape.x);
  float dy = (bbox.w - bbox.z) / float(shape.y);

  int imin = std::clamp(static_cast<int>(std::floor((xmin_roi - bbox.x) / dx)),
                        0,
                        shape.x - 1);
  int imax = std::clamp(static_cast<int>(std::ceil((xmax_roi - bbox.x) / dx)),
                        0,
                        shape.x - 1);
  int jmin = std::clamp(static_cast<int>(std::floor((ymin_roi - bbox.z) / dy)),
                        0,
                        shape.y - 1);
  int jmax = std::clamp(static_cast<int>(std::ceil((ymax_roi - bbox.z) / dy)),
                        0,
                        shape.y - 1);

#pragma omp parallel for schedule(dynamic, 16)
  for (int j = jmin; j <= jmax; ++j)
  {
    std::vector<size_t> indices;
    std::vector<float>  distances;
    indices.reserve(k_neighbors);
    distances.reserve(k_neighbors);

    for (int i = imin; i <= imax; ++i)
    {
      float xi = xg[i];
      float yi = yg[j];

      tree.neighbor_search(xi, yi, k_neighbors, indices, distances);

      float value = 0.f;
      float value_mask = 0.f;
      float dr = p_noise_r ? (*p_noise_r)(i, j) : 0.f;
      float effective_width = std::max(0.f, width * (1.f + dr));

      // use only 1st neighbor for various width scalings
      size_t k0 = indices[0];

      if (enable_width_distance_scaling)
        effective_width *= smoothstep3(arc_length[k0]);

      if (enable_width_depth_scaling)
      {
        float dz = std::abs((z(i, j) - path_copy[k0].v) / elevation_shift);
        effective_width *= std::clamp(dz, 0.f, 1.f);
      }

      if (enable_width_curvature_scaling)
      {
        // determine on which side of the path the cell is
        size_t km = std::max(k0 - 1, size_t(0));
        size_t kp = std::min(k0 + 1, npts);

        float s = -classify_point(path_copy[km],
                                  path_copy[k0],
                                  path_copy[kp],
                                  Point(xi, yi));

        // longitudinal scaling
        float t = arc_length[k0];
        t = t * (1.f - t) * 4.f;

        // float camp = std::abs(curvature[k0]) * t;
        float camp = std::abs(curv_radius[k0]) * t;

        effective_width *= lerp(
            1.f,
            std::max(curv_width_ratio_min,
                     1.f + (curv_width_ratio_max - 1.f) * s * camp),
            curv_shape_factor[k0]);
      }

      if (effective_width <= 0.f ||
          distances[0] > effective_width * effective_width)
      {
        continue;
      }

      for (size_t k = 0; k < k_neighbors; ++k)
      {
        float zref = path_copy[indices[k]].v;
        float r = std::sqrt(distances[k]) / effective_width;

        if (r >= 0.f && r <= 1.f)
        {
          float t = profile_fct(r);
          value += lerp(zref, z(i, j), t);
          value_mask += 1.f - t;
        }
        else
        {
          value += z(i, j);
        }
      }

      zp(i, j) = value / float(k_neighbors);
      blending_mask(i, j) = value_mask / float(k_neighbors);
    }
  }

  // --- outputs

  if (p_bending_mask) *p_bending_mask = std::move(blending_mask);

  z = std::move(zp);
}

} // namespace hmap
