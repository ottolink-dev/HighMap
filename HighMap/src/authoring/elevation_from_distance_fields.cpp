/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <limits>

#include "highmap/authoring.hpp"
#include "highmap/morphology.hpp"

namespace hmap
{

namespace
{

Array helper_compute_dt(const Array &mask, DistanceTransformType dt_type)
{
  switch (dt_type)
  {
  case DistanceTransformType::DT_APPROX: return distance_transform_approx(mask);
  case DistanceTransformType::DT_MANHATTAN:
    return distance_transform_manhattan(mask);
  case DistanceTransformType::DT_EXACT:
  default: return distance_transform(mask);
  }
}

bool helper_has_positive_values(const Array &arr)
{
  for (int j = 0; j < arr.shape.y; ++j)
    for (int i = 0; i < arr.shape.x; ++i)
      if (arr(i, j) > 0.f) return true;
  return false;
}

} // namespace

Array elevation_from_distance_fields(const Array          &mountains,
                                     const Array          *p_boundary,
                                     const Array          *p_coastline,
                                     float                 exponent,
                                     float                 smoothing_radius,
                                     float                 z_mountains,
                                     float                 z_boundary,
                                     float                 z_coastline,
                                     const Array          *p_noise,
                                     float                 noise_scale,
                                     DistanceTransformType dt_type)
{
  glm::ivec2 shape = mountains.shape;
  Array      elev(shape, 0.f);

  // --- Mountain distance field (A)

  bool  has_mountains = helper_has_positive_values(mountains);
  Array dA;

  if (has_mountains) dA = helper_compute_dt(mountains, dt_type);

  // --- Boundary distance field (B)

  bool  has_boundary = false;
  Array dB;

  if (p_boundary != nullptr && helper_has_positive_values(*p_boundary))
  {
    dB = helper_compute_dt(*p_boundary, dt_type);
    has_boundary = true;
  }
  else if (p_boundary == nullptr)
  {
    // Explicit radial distance from center (radius normalized by half diagonal
    // length)
    dB = Array(shape);
    float cx = 0.5f * static_cast<float>(shape.x - 1);
    float cy = 0.5f * static_cast<float>(shape.y - 1);
    float r_max = std::sqrt(cx * cx + cy * cy);
    if (r_max <= 0.f) r_max = 1.f;

    for (int j = 0; j < shape.y; ++j)
    {
      for (int i = 0; i < shape.x; ++i)
      {
        float dx = static_cast<float>(i) - cx;
        float dy = static_cast<float>(j) - cy;
        float r = std::sqrt(dx * dx + dy * dy);
        float d = r_max - r;
        dB(i, j) = (d <= 1e-5f) ? 0.f : d;
      }
    }
    has_boundary = true;
  }

  // --- Coastline distance field (C)

  bool  has_coastline = false;
  Array dC;

  if (p_coastline != nullptr && helper_has_positive_values(*p_coastline))
  {
    dC = helper_compute_dt(*p_coastline, dt_type);
    has_coastline = true;
  }

  // --- Combine distance fields using inverse-distance harmonic weighting

  float eps = std::max(0.f, smoothing_radius);
  float eps2 = eps * eps;

#pragma omp parallel for collapse(2)
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float val = 0.f;

      if (eps == 0.f)
      {
        // Exact piecewise check for zero-distance singularities
        bool on_m = has_mountains && (dA(i, j) == 0.f);
        bool on_c = has_coastline && (dC(i, j) == 0.f);
        bool on_b = has_boundary && (dB(i, j) == 0.f);

        if (on_m || on_c || on_b)
        {
          float sum_z = 0.f;
          int   count = 0;
          if (on_m)
          {
            sum_z += z_mountains;
            count++;
          }
          if (on_c)
          {
            sum_z += z_coastline;
            count++;
          }
          if (on_b)
          {
            sum_z += z_boundary;
            count++;
          }
          val = sum_z / static_cast<float>(count);
        }
        else
        {
          float wA = 0.f;
          float wB = 0.f;
          float wC = 0.f;
          float n_val = (p_noise != nullptr) ? (*p_noise)(i, j) * noise_scale
                                             : 0.f;

          if (has_mountains)
          {
            float distA = dA(i, j);
            if (n_val != 0.f)
            {
              distA = std::max(1e-5f, distA * (1.f + n_val));
            }
            wA = 1.f / std::pow(distA, exponent);
          }

          if (has_boundary)
          {
            float distB = dB(i, j);
            if (n_val != 0.f)
            {
              distB = std::max(1e-5f, distB * (1.f - n_val));
            }
            wB = 1.f / std::pow(distB, exponent);
          }

          if (has_coastline)
          {
            float distC = dC(i, j);
            if (n_val != 0.f)
            {
              distC = std::max(1e-5f, distC * (1.f + n_val));
            }
            wC = 1.f / std::pow(distC, exponent);
          }

          float sum_w = wA + wB + wC;
          if (sum_w > 0.f)
          {
            val = (z_mountains * wA + z_boundary * wB + z_coastline * wC) /
                  sum_w;
          }
        }
      }
      else
      {
        // Softened harmonic mean: smooth C^inf transitions without frontier
        // halos
        float wA = 0.f;
        float wB = 0.f;
        float wC = 0.f;
        float n_val = (p_noise != nullptr) ? (*p_noise)(i, j) * noise_scale
                                           : 0.f;

        if (has_mountains)
        {
          float distA = dA(i, j);
          if (n_val != 0.f)
          {
            distA = std::max(0.f, distA * (1.f + n_val));
          }
          float d2 = distA * distA + eps2;
          wA = 1.f / std::pow(d2, 0.5f * exponent);
        }

        if (has_boundary)
        {
          float distB = dB(i, j);
          if (n_val != 0.f)
          {
            distB = std::max(0.f, distB * (1.f - n_val));
          }
          float d2 = distB * distB + eps2;
          wB = 1.f / std::pow(d2, 0.5f * exponent);
        }

        if (has_coastline)
        {
          float distC = dC(i, j);
          if (n_val != 0.f)
          {
            distC = std::max(0.f, distC * (1.f + n_val));
          }
          float d2 = distC * distC + eps2;
          wC = 1.f / std::pow(d2, 0.5f * exponent);
        }

        float sum_w = wA + wB + wC;
        if (sum_w > 0.f)
        {
          val = (z_mountains * wA + z_boundary * wB + z_coastline * wC) / sum_w;
        }
      }

      elev(i, j) = val;
    }

  return elev;
}

} // namespace hmap
