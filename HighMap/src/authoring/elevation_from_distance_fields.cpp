/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <limits>

#include "highmap/authoring.hpp"
#include "highmap/math/core.hpp"
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
                                     float                 alpha,
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
    // Explicit radial distance from center using domain width and height
    // (semi-axes rx and ry) instead of the diagonal
    dB = Array(shape);
    float cx = 0.5f * static_cast<float>(shape.x - 1);
    float cy = 0.5f * static_cast<float>(shape.y - 1);
    float rx = (cx <= 0.f) ? 1.f : cx;
    float ry = (cy <= 0.f) ? 1.f : cy;

    for (int j = 0; j < shape.y; ++j)
    {
      for (int i = 0; i < shape.x; ++i)
      {
        float dx = static_cast<float>(i) - cx;
        float dy = static_cast<float>(j) - cy;
        float rho2 = (dx * dx) / (rx * rx) + (dy * dy) / (ry * ry);

        if (rho2 >= 1.f)
        {
          dB(i, j) = 0.f;
        }
        else
        {
          float r_euclid = std::sqrt(dx * dx + dy * dy);
          if (r_euclid <= 1e-5f)
          {
            dB(i, j) = std::min(rx, ry);
          }
          else
          {
            float rho = std::sqrt(rho2);
            float R_theta = r_euclid / rho;
            dB(i, j) = R_theta - r_euclid;
          }
        }
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

  // --- Combine distance fields using Rvachev R-functions

#pragma omp parallel for collapse(2)
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float n_val = (p_noise != nullptr) ? (*p_noise)(i, j) * noise_scale : 0.f;

      float distA = has_mountains ? std::max(0.f, dA(i, j) * (1.f + n_val))
                                  : 0.f;
      float distB = has_boundary ? std::max(0.f, dB(i, j) * (1.f - n_val))
                                 : 0.f;
      float distC = has_coastline ? std::max(0.f, dC(i, j) * (1.f + n_val))
                                  : 0.f;

      if (exponent != 1.f)
      {
        if (distA > 0.f) distA = std::pow(distA, exponent);
        if (distB > 0.f) distB = std::pow(distB, exponent);
        if (distC > 0.f) distC = std::pow(distC, exponent);
      }

      float wA = 0.f;
      float wB = 0.f;
      float wC = 0.f;

      if (has_mountains && has_boundary && has_coastline)
      {
        // Omega_{\neg k} using Rvachev Rmin intersection of the other distance
        // fields
        wA = r_min(distB, distC, alpha);
        wB = r_min(distA, distC, alpha);
        wC = r_min(distA, distB, alpha);
      }
      else if (has_mountains && has_boundary)
      {
        wA = distB;
        wB = distA;
      }
      else if (has_mountains && has_coastline)
      {
        wA = distC;
        wC = distA;
      }
      else if (has_coastline && has_boundary)
      {
        wC = distB;
        wB = distC;
      }
      else if (has_mountains)
      {
        wA = 1.f;
      }

      float sum_w = wA + wB + wC;
      if (sum_w > 0.f)
      {
        elev(i, j) = (z_mountains * wA + z_boundary * wB + z_coastline * wC) /
                     sum_w;
      }
      else
      {
        // If constraints intersect at (i, j) where all weights are zero,
        // average them
        float sum_z = 0.f;
        int   count = 0;
        if (has_mountains && distA == 0.f)
        {
          sum_z += z_mountains;
          count++;
        }
        if (has_boundary && distB == 0.f)
        {
          sum_z += z_boundary;
          count++;
        }
        if (has_coastline && distC == 0.f)
        {
          sum_z += z_coastline;
          count++;
        }
        elev(i, j) = (count > 0) ? (sum_z / static_cast<float>(count)) : 0.f;
      }
    }

  return elev;
}

} // namespace hmap
