/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/primitives/coherent_noise.hpp"
#include "highmap/primitives/functions.hpp"
#include "highmap/primitives/geo.hpp"
#include "highmap/range.hpp"

namespace hmap::gpu
{

// --- helpers

void helper_bind_optional_buffers(clwrapper::Run &run,
                                  const Array    *p_noise_x,
                                  const Array    *p_noise_y)
{
  std::vector<float> dummy_vector(1);

  if (p_noise_x)
  {
    run.bind_buffer<float>("noise_x", p_noise_x->vector);
    run.write_buffer("noise_x");
  }
  else
    run.bind_buffer<float>("noise_x", dummy_vector);

  if (p_noise_y)
  {
    run.bind_buffer<float>("noise_y", p_noise_y->vector);
    run.write_buffer("noise_y");
  }
  else
    run.bind_buffer<float>("noise_y", dummy_vector);
}

// --- functions

Array gabor_wave(glm::ivec2    shape,
                 glm::vec2     kw,
                 std::uint32_t seed,
                 const Array  &angle,
                 float         angle_spread_ratio,
                 glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (!validate_same_shape(shape, angle)) return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("gabor_wave");

  run.bind_buffer<float>("array", array.vector);
  run.bind_buffer<float>("angle", angle.vector);
  run.write_buffer("angle");
  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     angle_spread_ratio,
                     bbox);

  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array gabor_wave(glm::ivec2    shape,
                 glm::vec2     kw,
                 std::uint32_t seed,
                 float         angle,
                 float         angle_spread_ratio,
                 glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();

  Array array(shape);
  Array array_angle(shape, angle);

  array = gabor_wave(shape, kw, seed, array_angle, angle_spread_ratio, bbox);

  return array;
}

Array gabor_wave_fbm(glm::ivec2    shape,
                     glm::vec2     kw,
                     std::uint32_t seed,
                     const Array  &angle,
                     float         angle_spread_ratio,
                     int           octaves,
                     float         weight,
                     float         persistence,
                     float         lacunarity,
                     const Array  *p_ctrl_param,
                     const Array  *p_noise_x,
                     const Array  *p_noise_y,
                     glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (!validate_same_shape(shape, angle)) return Array(shape);
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("gabor_wave_fbm");

  run.bind_buffer<float>("array", array.vector);
  run.bind_buffer<float>("angle", angle.vector);

  helper_bind_optional_buffer(run, "ctrl_param", p_ctrl_param);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     angle_spread_ratio,
                     octaves,
                     weight,
                     persistence,
                     lacunarity,
                     p_ctrl_param ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.write_buffer("angle");
  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array gabor_wave_fbm(glm::ivec2    shape,
                     glm::vec2     kw,
                     std::uint32_t seed,
                     float         angle,
                     float         angle_spread_ratio,
                     int           octaves,
                     float         weight,
                     float         persistence,
                     float         lacunarity,
                     const Array  *p_ctrl_param,
                     const Array  *p_noise_x,
                     const Array  *p_noise_y,
                     glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array(shape);
  Array array_angle(shape, angle);

  array = gabor_wave_fbm(shape,
                         kw,
                         seed,
                         array_angle,
                         angle_spread_ratio,
                         octaves,
                         weight,
                         persistence,
                         lacunarity,
                         p_ctrl_param,
                         p_noise_x,
                         p_noise_y,
                         bbox);

  return array;
}

Array gavoronoise(glm::ivec2    shape,
                  glm::vec2     kw,
                  std::uint32_t seed,
                  const Array  &angle,
                  float         amplitude,
                  float         angle_spread_ratio,
                  glm::vec2     kw_multiplier,
                  float         slope_strength,
                  float         branch_strength,
                  float         z_cut_min,
                  float         z_cut_max,
                  int           octaves,
                  float         persistence,
                  float         lacunarity,
                  const Array  *p_ctrl_param,
                  const Array  *p_noise_x,
                  const Array  *p_noise_y,
                  glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (!validate_same_shape(shape, angle)) return Array(shape);
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("gavoronoise");

  run.bind_buffer<float>("array", array.vector);
  run.bind_buffer<float>("angle",
                         const_cast<std::vector<float> &>(angle.vector));

  helper_bind_optional_buffer(run, "ctrl_param", p_ctrl_param);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     amplitude,
                     angle_spread_ratio,
                     kw_multiplier,
                     slope_strength,
                     branch_strength,
                     z_cut_min,
                     z_cut_max,
                     octaves,
                     persistence,
                     lacunarity,
                     p_ctrl_param ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.write_buffer("angle");
  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array gavoronoise(glm::ivec2    shape,
                  glm::vec2     kw,
                  std::uint32_t seed,
                  float         angle,
                  float         amplitude,
                  float         angle_spread_ratio,
                  glm::vec2     kw_multiplier,
                  float         slope_strength,
                  float         branch_strength,
                  float         z_cut_min,
                  float         z_cut_max,
                  int           octaves,
                  float         persistence,
                  float         lacunarity,
                  const Array  *p_ctrl_param,
                  const Array  *p_noise_x,
                  const Array  *p_noise_y,
                  glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array(shape);
  Array array_angle(shape, angle);

  array = gavoronoise(shape,
                      kw,
                      seed,
                      array_angle,
                      amplitude,
                      angle_spread_ratio,
                      kw_multiplier,
                      slope_strength,
                      branch_strength,
                      z_cut_min,
                      z_cut_max,
                      octaves,
                      persistence,
                      lacunarity,
                      p_ctrl_param,
                      p_noise_x,
                      p_noise_y,
                      bbox);

  return array;
}

Array gavoronoise(const Array  &base,
                  glm::vec2     kw,
                  std::uint32_t seed,
                  float         amplitude,
                  glm::vec2     kw_multiplier,
                  float         slope_strength,
                  float         branch_strength,
                  float         z_cut_min,
                  float         z_cut_max,
                  int           octaves,
                  float         persistence,
                  float         lacunarity,
                  const Array  *p_ctrl_param,
                  const Array  *p_noise_x,
                  const Array  *p_noise_y,
                  glm::vec4     bbox)
{
  if (!validate_non_empty(base)) return Array();
  if (p_ctrl_param && !validate_same_shape(base, *p_ctrl_param))
    return Array(base.shape);
  if (p_noise_x && !validate_same_shape(base, *p_noise_x))
    return Array(base.shape);
  if (p_noise_y && !validate_same_shape(base, *p_noise_y))
    return Array(base.shape);

  Array array(base.shape);

  auto run = clwrapper::Run("gavoronoise_with_base");

  run.bind_imagef("base",
                  const_cast<std::vector<float> &>(base.vector),
                  base.shape.x,
                  base.shape.y);
  run.bind_buffer<float>("array", array.vector);

  helper_bind_optional_buffer(run, "ctrl_param", p_ctrl_param);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     amplitude,
                     kw_multiplier,
                     slope_strength,
                     branch_strength,
                     z_cut_min,
                     z_cut_max,
                     octaves,
                     persistence,
                     lacunarity,
                     p_ctrl_param ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array hemisphere_field(glm::ivec2    shape,
                       glm::vec2     kw,
                       std::uint32_t seed,
                       float         rmin,
                       float         rmax,
                       float         amplitude_random_ratio,
                       float         density,
                       glm::vec2     jitter,
                       float         shift,
                       const Array  *p_noise_x,
                       const Array  *p_noise_y,
                       const Array  *p_noise_distance,
                       const Array  *p_density_multiplier,
                       const Array  *p_size_multiplier,
                       glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);
  if (p_noise_distance && !validate_same_shape(shape, *p_noise_distance))
    return Array(shape);
  if (p_density_multiplier &&
      !validate_same_shape(shape, *p_density_multiplier))
    return Array(shape);
  if (p_size_multiplier && !validate_same_shape(shape, *p_size_multiplier))
    return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("hemisphere_field");

  run.bind_buffer<float>("array", array.vector);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);
  helper_bind_optional_buffer(run, "p_noise_distance", p_noise_distance);
  helper_bind_optional_buffer(run,
                              "p_density_multiplier",
                              p_density_multiplier);
  helper_bind_optional_buffer(run, "p_size_multiplier", p_size_multiplier);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     rmin,
                     rmax,
                     amplitude_random_ratio,
                     density,
                     jitter,
                     shift,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     p_noise_distance ? 1 : 0,
                     p_density_multiplier ? 1 : 0,
                     p_size_multiplier ? 1 : 0,
                     bbox);

  run.write_buffer("array");
  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array hemisphere_field_fbm(glm::ivec2    shape,
                           glm::vec2     kw,
                           std::uint32_t seed,
                           float         rmin,
                           float         rmax,
                           float         amplitude_random_ratio,
                           float         density,
                           glm::vec2     jitter,
                           float         shift,
                           int           octaves,
                           float         persistence,
                           float         lacunarity,
                           const Array  *p_noise_x,
                           const Array  *p_noise_y,
                           const Array  *p_noise_distance,
                           const Array  *p_density_multiplier,
                           const Array  *p_size_multiplier,
                           glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);
  if (p_noise_distance && !validate_same_shape(shape, *p_noise_distance))
    return Array(shape);
  if (p_density_multiplier &&
      !validate_same_shape(shape, *p_density_multiplier))
    return Array(shape);
  if (p_size_multiplier && !validate_same_shape(shape, *p_size_multiplier))
    return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("hemisphere_field_fbm");

  run.bind_buffer<float>("array", array.vector);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);
  helper_bind_optional_buffer(run, "p_noise_distance", p_noise_distance);
  helper_bind_optional_buffer(run,
                              "p_density_multiplier",
                              p_density_multiplier);
  helper_bind_optional_buffer(run, "p_size_multiplier", p_size_multiplier);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     rmin,
                     rmax,
                     amplitude_random_ratio,
                     density,
                     jitter,
                     shift,
                     octaves,
                     persistence,
                     lacunarity,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     p_noise_distance ? 1 : 0,
                     p_density_multiplier ? 1 : 0,
                     p_size_multiplier ? 1 : 0,
                     bbox);

  run.write_buffer("array");
  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array mountain_range_radial(glm::ivec2    shape,
                            glm::vec2     kw,
                            std::uint32_t seed,
                            float         half_width,
                            float         angle_spread_ratio,
                            float         core_size_ratio,
                            glm::vec2     center,
                            int           octaves,
                            float         weight,
                            float         persistence,
                            float         lacunarity,
                            const Array  *p_ctrl_param,
                            const Array  *p_noise_x,
                            const Array  *p_noise_y,
                            const Array  *p_angle,
                            glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);
  if (p_angle && !validate_same_shape(shape, *p_angle)) return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("mountain_range_radial");

  run.bind_buffer<float>("array", array.vector);

  helper_bind_optional_buffer(run, "ctrl_param", p_ctrl_param);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);
  helper_bind_optional_buffer(run, "angle", p_angle);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     half_width,
                     angle_spread_ratio,
                     core_size_ratio,
                     center,
                     octaves,
                     weight,
                     persistence,
                     lacunarity,
                     p_ctrl_param ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     p_angle ? 1 : 0,
                     bbox);

  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  if (p_angle) run.read_buffer("angle");

  return array;
}

Array noise(NoiseType     noise_type,
            glm::ivec2    shape,
            glm::vec2     kw,
            std::uint32_t seed,
            const Array  *p_noise_x,
            const Array  *p_noise_y,
            glm::vec4     bbox,
            glm::ivec2    period)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array(shape);

  int noise_id = static_cast<int>(noise_type);

  auto run = clwrapper::Run("noise");

  run.bind_buffer<float>("array", array.vector);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     noise_id,
                     kw.x,
                     kw.y,
                     seed,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     period.x,
                     period.y,
                     bbox);

  run.write_buffer("array");
  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array noise_fbm(NoiseType     noise_type,
                glm::ivec2    shape,
                glm::vec2     kw,
                std::uint32_t seed,
                int           octaves,
                float         weight,
                float         persistence,
                float         lacunarity,
                const Array  *p_ctrl_param,
                const Array  *p_noise_x,
                const Array  *p_noise_y,
                glm::vec4     bbox,
                glm::ivec2    period)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array(shape);

  int noise_id = static_cast<int>(noise_type);

  auto run = clwrapper::Run("noise_fbm");

  run.bind_buffer<float>("array", array.vector);
  helper_bind_optional_buffer(run, "ctrl_param", p_ctrl_param);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     noise_id,
                     kw.x,
                     kw.y,
                     seed,
                     octaves,
                     weight,
                     persistence,
                     lacunarity,
                     p_ctrl_param ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     period.x,
                     period.y,
                     bbox);

  run.write_buffer("array");
  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array polygon_field(glm::ivec2    shape,
                    glm::vec2     kw,
                    std::uint32_t seed,
                    float         rmin,
                    float         rmax,
                    float         clamping_dist,
                    float         clamping_k,
                    int           n_vertices_min,
                    int           n_vertices_max,
                    float         density,
                    glm::vec2     jitter,
                    float         shift,
                    const Array  *p_noise_x,
                    const Array  *p_noise_y,
                    const Array  *p_noise_distance,
                    const Array  *p_density_multiplier,
                    const Array  *p_size_multiplier,
                    glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);
  if (p_noise_distance && !validate_same_shape(shape, *p_noise_distance))
    return Array(shape);
  if (p_density_multiplier &&
      !validate_same_shape(shape, *p_density_multiplier))
    return Array(shape);
  if (p_size_multiplier && !validate_same_shape(shape, *p_size_multiplier))
    return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("polygon_field");

  run.bind_buffer<float>("array", array.vector);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);
  helper_bind_optional_buffer(run, "p_noise_distance", p_noise_distance);
  helper_bind_optional_buffer(run,
                              "p_density_multiplier",
                              p_density_multiplier);
  helper_bind_optional_buffer(run, "p_size_multiplier", p_size_multiplier);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     rmin,
                     rmax,
                     clamping_dist,
                     clamping_k,
                     n_vertices_min,
                     n_vertices_max,
                     density,
                     jitter,
                     shift,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     p_noise_distance ? 1 : 0,
                     p_density_multiplier ? 1 : 0,
                     p_size_multiplier ? 1 : 0,
                     bbox);

  run.write_buffer("array");
  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array polygon_field_fbm(glm::ivec2    shape,
                        glm::vec2     kw,
                        std::uint32_t seed,
                        float         rmin,
                        float         rmax,
                        float         clamping_dist,
                        float         clamping_k,
                        int           n_vertices_min,
                        int           n_vertices_max,
                        float         density,
                        glm::vec2     jitter,
                        float         shift,
                        int           octaves,
                        float         persistence,
                        float         lacunarity,
                        const Array  *p_noise_x,
                        const Array  *p_noise_y,
                        const Array  *p_noise_distance,
                        const Array  *p_density_multiplier,
                        const Array  *p_size_multiplier,
                        glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);
  if (p_noise_distance && !validate_same_shape(shape, *p_noise_distance))
    return Array(shape);
  if (p_density_multiplier &&
      !validate_same_shape(shape, *p_density_multiplier))
    return Array(shape);
  if (p_size_multiplier && !validate_same_shape(shape, *p_size_multiplier))
    return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("polygon_field_fbm");

  run.bind_buffer<float>("array", array.vector);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);
  helper_bind_optional_buffer(run, "p_noise_distance", p_noise_distance);
  helper_bind_optional_buffer(run,
                              "p_density_multiplier",
                              p_density_multiplier);
  helper_bind_optional_buffer(run, "p_size_multiplier", p_size_multiplier);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     rmin,
                     rmax,
                     clamping_dist,
                     clamping_k,
                     n_vertices_min,
                     n_vertices_max,
                     density,
                     jitter,
                     shift,
                     octaves,
                     persistence,
                     lacunarity,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     p_noise_distance ? 1 : 0,
                     p_density_multiplier ? 1 : 0,
                     p_size_multiplier ? 1 : 0,
                     bbox);

  run.write_buffer("array");
  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array vorolines(glm::ivec2        shape,
                float             density,
                std::uint32_t     seed,
                float             k_smoothing,
                float             exp_sigma,
                float             alpha,
                float             alpha_span,
                VoronoiReturnType return_type,
                const Array      *p_noise_x,
                const Array      *p_noise_y,
                glm::vec4         bbox,
                glm::vec4         bbox_points)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  // --- generate random set of points

  // density is the number of pts per unit surface
  int npoints = static_cast<int>(density * (bbox_points.y - bbox_points.x) *
                                 (bbox_points.w - bbox_points.z));
  npoints = std::max(1, npoints);
  Cloud cloud = Cloud(npoints, seed, bbox_points);

  std::vector<float> xp = cloud.get_x();
  std::vector<float> yp = cloud.get_y();
  std::vector<float> v = cloud.get_values();

  // generate secondary set of points to define lines
  std::vector<float> xs, ys;
  xs.reserve(xp.size());
  ys.reserve(xp.size());

  for (size_t k = 0; k < v.size(); ++k)
  {
    // use random values at points generated by the cloud to
    // determine line angles
    float theta = alpha + (2.f * v[k] - 1.f) * alpha_span;
    float x = xp[k] + std::cos(theta);
    float y = yp[k] + std::sin(theta);

    xs.push_back(x);
    ys.push_back(y);
  }

  // --- generate

  Array array(shape);

  auto run = clwrapper::Run("vorolines");

  run.bind_buffer<float>("array", array.vector);

  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_buffer<float>("xp", xp);
  run.bind_buffer<float>("yp", yp);
  run.bind_buffer<float>("xs", xs);
  run.bind_buffer<float>("ys", ys);

  run.write_buffer("xp");
  run.write_buffer("yp");
  run.write_buffer("xs");
  run.write_buffer("ys");

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     static_cast<int>(xp.size()),
                     k_smoothing,
                     exp_sigma,
                     (int)return_type,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array vorolines_fbm(glm::ivec2        shape,
                    float             density,
                    std::uint32_t     seed,
                    float             k_smoothing,
                    float             exp_sigma,
                    float             alpha,
                    float             alpha_span,
                    VoronoiReturnType return_type,
                    int               octaves,
                    float             weight,
                    float             persistence,
                    float             lacunarity,
                    const Array      *p_noise_x,
                    const Array      *p_noise_y,
                    glm::vec4         bbox,
                    glm::vec4         bbox_points)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array n = Array(shape);
  Array na = Array(shape, 0.6f);
  float nf = 1.f;

  for (int i = 0; i < octaves; i++)
  {
    Array v = vorolines(shape,
                        nf * density,
                        seed++,
                        k_smoothing,
                        exp_sigma,
                        alpha,
                        alpha_span,
                        return_type,
                        p_noise_x,
                        p_noise_y,
                        bbox,
                        bbox_points);

    n += v * na;
    na *= (1.f - weight) + weight * minimum(v + 1.f, 2.f) * 0.5f;
    na *= persistence;
    nf *= lacunarity;
  }

  return n;
}

Array voronoi(glm::ivec2        shape,
              glm::vec2         kw,
              std::uint32_t     seed,
              glm::vec2         jitter,
              float             bias,
              float             k_smoothing,
              float             exp_sigma,
              VoronoiReturnType return_type,
              const Array      *p_ctrl_param,
              const Array      *p_noise_x,
              const Array      *p_noise_y,
              glm::vec4         bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("voronoi");

  run.bind_buffer<float>("array", array.vector);

  helper_bind_optional_buffer(run, "ctrl_param", p_ctrl_param);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     jitter,
                     bias,
                     k_smoothing,
                     exp_sigma,
                     (int)return_type,
                     p_ctrl_param ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array voronoi_fbm(glm::ivec2        shape,
                  glm::vec2         kw,
                  std::uint32_t     seed,
                  glm::vec2         jitter,
                  float             bias,
                  float             k_smoothing,
                  float             exp_sigma,
                  VoronoiReturnType return_type,
                  int               octaves,
                  float             weight,
                  float             persistence,
                  float             lacunarity,
                  const Array      *p_ctrl_param,
                  const Array      *p_noise_x,
                  const Array      *p_noise_y,
                  glm::vec4         bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("voronoi_fbm");

  run.bind_buffer<float>("array", array.vector);

  helper_bind_optional_buffer(run, "ctrl_param", p_ctrl_param);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     jitter,
                     bias,
                     k_smoothing,
                     exp_sigma,
                     (int)return_type,
                     octaves,
                     weight,
                     persistence,
                     lacunarity,
                     p_ctrl_param ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array voronoise(glm::ivec2    shape,
                glm::vec2     kw,
                float         u_param,
                float         v_param,
                std::uint32_t seed,
                const Array  *p_noise_x,
                const Array  *p_noise_y,
                glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("voronoise");

  run.bind_buffer<float>("array", array.vector);

  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     u_param,
                     v_param,
                     seed,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array voronoise_fbm(glm::ivec2    shape,
                    glm::vec2     kw,
                    float         u_param,
                    float         v_param,
                    std::uint32_t seed,
                    int           octaves,
                    float         weight,
                    float         persistence,
                    float         lacunarity,
                    const Array  *p_ctrl_param,
                    const Array  *p_noise_x,
                    const Array  *p_noise_y,
                    glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("voronoise_fbm");

  run.bind_buffer<float>("array", array.vector);

  helper_bind_optional_buffer(run, "ctrl_param", p_ctrl_param);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     u_param,
                     v_param,
                     seed,
                     octaves,
                     weight,
                     persistence,
                     lacunarity,
                     p_ctrl_param ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array voronoi_edge_distance(glm::ivec2    shape,
                            glm::vec2     kw,
                            std::uint32_t seed,
                            glm::vec2     jitter,
                            float         bias,
                            const Array  *p_ctrl_param,
                            const Array  *p_noise_x,
                            const Array  *p_noise_y,
                            glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("voronoi_edge_distance");

  run.bind_buffer<float>("array", array.vector);

  helper_bind_optional_buffer(run, "ctrl_param", p_ctrl_param);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     jitter,
                     bias,
                     p_ctrl_param ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array vororand(glm::ivec2        shape,
               float             density,
               float             variability,
               std::uint32_t     seed,
               float             k_smoothing,
               float             exp_sigma,
               VoronoiReturnType return_type,
               const Array      *p_noise_x,
               const Array      *p_noise_y,
               glm::vec4         bbox,
               glm::vec4         bbox_points)
{
  if (!validate_shape(shape)) return Array();
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  // --- generate random set of points

  // TODO adjust extension with respect to the density?

  // take a bounding box a bit larger to reduce border effects
  float     lx = variability * (bbox_points.y - bbox_points.x);
  float     ly = variability * (bbox_points.w - bbox_points.z);
  glm::vec4 bbox_points_mod = adjust(bbox_points, -lx, lx, -ly, ly);

  // density is the number of pts per unit surface
  int npoints = static_cast<int>(density *
                                 (bbox_points_mod.y - bbox_points_mod.x) *
                                 (bbox_points_mod.w - bbox_points_mod.z));
  npoints = std::max(1, npoints);
  Cloud cloud = Cloud(npoints, seed, bbox_points_mod);

  std::vector<float> xp = cloud.get_x();
  std::vector<float> yp = cloud.get_y();

  // --- generate noise

  Array array = vororand(shape,
                         xp,
                         yp,
                         k_smoothing,
                         exp_sigma,
                         return_type,
                         p_noise_x,
                         p_noise_y,
                         bbox);

  return array;
}

Array vororand(glm::ivec2                shape,
               const std::vector<float> &xp,
               const std::vector<float> &yp,
               float                     k_smoothing,
               float                     exp_sigma,
               VoronoiReturnType         return_type,
               const Array              *p_noise_x,
               const Array              *p_noise_y,
               glm::vec4                 bbox)
{
  if (!validate_shape(shape)) return Array();
  if (!validate_non_empty(xp, "xp") || !validate_non_empty(yp, "yp") ||
      xp.size() != yp.size())
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("vororand");

  run.bind_buffer<float>("array", array.vector);

  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_buffer<float>("xp", xp);
  run.bind_buffer<float>("yp", yp);

  run.write_buffer("xp");
  run.write_buffer("yp");

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     static_cast<int>(xp.size()),
                     k_smoothing,
                     exp_sigma,
                     (int)return_type,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

Array wavelet_noise(glm::ivec2    shape,
                    glm::vec2     kw,
                    std::uint32_t seed,
                    float         kw_multiplier,
                    float         vorticity,
                    float         density,
                    int           octaves,
                    float         weight,
                    float         persistence,
                    float         lacunarity,
                    const Array  *p_ctrl_param,
                    const Array  *p_noise_x,
                    const Array  *p_noise_y,
                    glm::vec4     bbox)
{
  if (!validate_shape(shape)) return Array();
  if (p_ctrl_param && !validate_same_shape(shape, *p_ctrl_param))
    return Array(shape);
  if (p_noise_x && !validate_same_shape(shape, *p_noise_x)) return Array(shape);
  if (p_noise_y && !validate_same_shape(shape, *p_noise_y)) return Array(shape);

  Array array(shape);

  auto run = clwrapper::Run("wavelet_noise");

  run.bind_buffer<float>("array", array.vector);
  helper_bind_optional_buffer(run, "ctrl_param", p_ctrl_param);
  helper_bind_optional_buffer(run, "noise_x", p_noise_x);
  helper_bind_optional_buffer(run, "noise_y", p_noise_y);

  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kw.x,
                     kw.y,
                     seed,
                     kw_multiplier,
                     vorticity,
                     density,
                     octaves,
                     weight,
                     persistence,
                     lacunarity,
                     p_ctrl_param ? 1 : 0,
                     p_noise_x ? 1 : 0,
                     p_noise_y ? 1 : 0,
                     bbox);

  run.execute({array.shape.x, array.shape.y});
  run.read_buffer("array");

  return array;
}

} // namespace hmap::gpu
