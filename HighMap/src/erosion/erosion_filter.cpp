/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <cstdint>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/erosion.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/opencl/gpu_opencl.hpp"

namespace hmap::gpu
{

// --- Runevision Advanced Terrain Erosion Filter

void erosion_filter(Array        &z,
                    float         scale,
                    float         strength,
                    float         gully_weight,
                    float         detail,
                    glm::vec4     rounding,
                    glm::vec4     onset,
                    glm::vec2     assumed_slope,
                    float         cell_scale,
                    int           octaves,
                    float         gain,
                    float         lacunarity,
                    float         normalization,
                    std::uint32_t seed,
                    const Array  *p_fade_target,
                    Array        *p_ridge_map,
                    glm::vec4     bbox)
{
  if (!validate_non_empty(z)) return;
  if (p_fade_target && !validate_same_shape(z, *p_fade_target)) return;
  if (p_ridge_map && !validate_same_shape(z, *p_ridge_map))
  {
    *p_ridge_map = Array(z.shape);
  }

  const glm::ivec2 shape = z.shape;

  auto run = clwrapper::Run("erosion_filter");

  // input heightmap as image2d_t
  run.bind_imagef("height_in", z.vector, shape.x, shape.y);

  // output heightmap buffer
  run.bind_buffer<float>("output_height", z.vector);

  // output ridge map buffer (or dummy if not requested)
  std::vector<float> dummy_ridge(1, 0.f);
  if (p_ridge_map)
  {
    run.bind_buffer<float>("output_ridge_map", p_ridge_map->vector);
  }
  else
  {
    run.bind_buffer<float>("output_ridge_map", dummy_ridge);
  }

  // optional fade target buffer
  helper_bind_optional_buffer(run, "fade_target_in", p_fade_target);

  run.bind_arguments(shape.x,
                     shape.y,
                     seed,
                     scale,
                     strength,
                     gully_weight,
                     detail,
                     rounding,
                     onset,
                     assumed_slope,
                     cell_scale,
                     octaves,
                     gain,
                     lacunarity,
                     normalization,
                     p_fade_target ? 1 : 0,
                     p_ridge_map ? 1 : 0,
                     bbox);

  run.execute({shape.x, shape.y});

  run.read_buffer("output_height");
  if (p_ridge_map)
  {
    run.read_buffer("output_ridge_map");
  }
}

} // namespace hmap::gpu
