/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/array.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/math/array.hpp"
#include "highmap/synthesis.hpp"
#include "highmap/transform.hpp"

namespace hmap
{

namespace
{

// Block-mean downsampling to produce the low-resolution guide.
Array downsample_block_mean(const Array &array, int factor)
{
  if (factor <= 1) return array;

  Array out(glm::ivec2(array.shape.x / factor, array.shape.y / factor));

  for (int j = 0; j < out.shape.y; j++)
    for (int i = 0; i < out.shape.x; i++)
    {
      float sum = 0.f;
      for (int q = 0; q < factor; q++)
        for (int p = 0; p < factor; p++)
          sum += array(i * factor + p, j * factor + q);
      out(i, j) = sum / (float)(factor * factor);
    }

  return out;
}

} // namespace

Array loose_symmetry(const Array &array,
                     SymmetryType symmetry_type,
                     float        strength,
                     int          factor,
                     int          patch_size,
                     int          analysis_stride,
                     int          synthesis_stride,
                     int          sparsity,
                     bool         flatten_center,
                     float        flatten_radius)
{
  if (!validate_non_empty(array)) return Array();

  if (factor < 1)
  {
    hmap::log::warn("loose_symmetry: factor must be >= 1 (got {})", factor);
    return Array();
  }

  // --- Macro base downsampling

  Array base = downsample_block_mean(array, factor);

  // --- Symmetrical guide construction with tunable strength

  Array symmetrical = symmetrize(base,
                                 symmetry_type,
                                 flatten_center,
                                 flatten_radius);
  Array guide = (strength >= 1.f)   ? symmetrical
                : (strength <= 0.f) ? base
                                    : lerp(base, symmetrical, strength);

  // --- Exemplar-based feature transfer

  return terrain_super_resolution(guide,
                                  array,
                                  factor,
                                  patch_size,
                                  analysis_stride,
                                  synthesis_stride,
                                  sparsity);
}

} // namespace hmap
