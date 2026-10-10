/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/array.hpp"
#include "highmap/authoring.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate2d.hpp"

namespace hmap
{

Array elevation_from_sparse_constraints(const Array &mountains,
                                        const Array *p_coastline,
                                        const Array *p_boundary,
                                        float        mountain_elevation,
                                        float        coastline_elevation,
                                        float        boundary_elevation,
                                        int          iterations_max,
                                        float        tolerance,
                                        const Array *p_noise,
                                        float        noise_amplitude)
{
  if (!validate_non_empty(mountains)) return Array();

  if (p_coastline && !validate_same_shape(mountains, *p_coastline))
    return Array();
  if (p_boundary && !validate_same_shape(mountains, *p_boundary))
    return Array();
  if (p_noise && !validate_same_shape(mountains, *p_noise)) return Array();

  const glm::ivec2 shape = mountains.shape;
  Array            z_init(shape, 0.0f);
  Array            mask_fixed(shape, 0.0f);

  // 1. Boundary constraint (outer perimeter default if not provided)
  if (p_boundary)
  {
    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        if ((*p_boundary)(i, j) > 0.0f)
        {
          z_init(i, j) = boundary_elevation;
          mask_fixed(i, j) = 1.0f;
        }
      }
  }
  else
  {
    // Default: fixed border along domain edges
    for (int i = 0; i < shape.x; ++i)
    {
      z_init(i, 0) = boundary_elevation;
      mask_fixed(i, 0) = 1.0f;
      z_init(i, shape.y - 1) = boundary_elevation;
      mask_fixed(i, shape.y - 1) = 1.0f;
    }
    for (int j = 0; j < shape.y; ++j)
    {
      z_init(0, j) = boundary_elevation;
      mask_fixed(0, j) = 1.0f;
      z_init(shape.x - 1, j) = boundary_elevation;
      mask_fixed(shape.x - 1, j) = 1.0f;
    }
  }

  // 2. Coastline constraint (z = coastline_elevation)
  if (p_coastline)
  {
    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        if ((*p_coastline)(i, j) > 0.0f)
        {
          z_init(i, j) = coastline_elevation;
          mask_fixed(i, j) = 1.0f;
        }
      }
  }

  // 3. Mountains constraint (z = mountain_elevation * amplitude)
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      float m = mountains(i, j);
      if (m > 0.0f)
      {
        z_init(i, j) = mountain_elevation * m;
        mask_fixed(i, j) = 1.0f;
      }
    }

  // 4. Solve the Laplace equation Δz = 0 on GPU
  Array z = gpu::harmonic_interpolation(z_init,
                                        mask_fixed,
                                        iterations_max,
                                        tolerance);

  // 5. Optional procedural detail addition modulated away from exact
  // constraints
  if (p_noise && noise_amplitude > 0.0f)
  {
    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
      {
        z(i, j) += noise_amplitude * (*p_noise)(i, j);
      }
  }

  return z;
}

} // namespace hmap
