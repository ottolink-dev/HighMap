/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file coherent_noise.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2023
 */
#pragma once

#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/interpolate/interpolate2d.hpp"
#include "highmap/math/profiles.hpp"

#define HMAP_GRADIENT_OFFSET 0.001f

namespace hmap
{

/**
 * @enum VoronoiReturnType
 * @brief Selects the value returned by the Voronoi evaluation.
 *
 * F1 is the distance to the nearest site, F2 to the second nearest.
 */
enum VoronoiReturnType : int
{
  F1_SQUARED,            ///< Returns F1^2
  F2_SQUARED,            ///< Returns F2^2
  F1TF2_SQUARED,         ///< Returns (F1 * F2)^2
  F1DF2_SQUARED,         ///< Returns (F1 / F2)^2
  F2MF1_SQUARED,         ///< Returns (F2 - F1)^2
  EDGE_DISTANCE_EXP,     ///< Exponential edge distance
  EDGE_DISTANCE_SQUARED, ///< Squared edge distance
  CONSTANT,              ///< Constant value
  CONSTANT_F2MF1_SQUARED ///< Constant × (F2 - F1)^2
};

/**
 * @brief Dendry is a locally computable procedural function that generates
 * branching patterns at various scales (see @cite Gaillard2019).
 *
 * @param  shape                       Array shape.
 * @param  kw                          Noise wavenumber with respect to a unit
 *                                     domain.
 * @param  seed                        Random seed number.
 * @param  control_array               Control array (can be of any shape,
 *                                     different from
 * `shape`).
 * @param  eps                         Epsilon used to bias the area where
 *                                     points are generated in cells.
 * @param  resolution                  Number of resolutions in the noise
 *                                     function.
 * @param  displacement                Maximum displacement of segments.
 * @param  primitives_resolution_steps Additional resolution steps in the
 *                                     ComputeColorPrimitives function.
 * @param  slope_power                 Additional parameter to control the
 *                                     variation of slope on terrains.
 * @param  noise_amplitude_proportion  Proportion of the amplitude of the
 *                                     control function as noise.
 * @param  add_control_function        Add control function to the output.
 * @param  control_function_overlap    Extent of the extension added at the
 *                                     domain frontiers of the control array.
 * @param  p_noise_x, p_noise_y        Reference to the input noise arrays.
 * @param  p_stretching                Local wavenumber multiplier.
 * @param  bbox                        Domain bounding box.
 * @return                             Array New array.
 *
 * **Example**
 * @include ex_dendry.cpp
 *
 * **Result**
 * @image html ex_dendry.png
 */
Array dendry(glm::ivec2    shape,
             glm::vec2     kw,
             std::uint32_t seed,
             Array        &control_function,
             float         eps = 0.05,
             int           resolution = 1,
             float         displacement = 0.075,
             int           primitives_resolution_steps = 3,
             float         slope_power = 2.f,
             float         noise_amplitude_proportion = 0.01,
             bool          add_control_function = true,
             float         control_function_overlap = 0.5f,
             const Array  *p_noise_x = nullptr,
             const Array  *p_noise_y = nullptr,
             glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f},
             int           subsampling = 1);

Array dendry(glm::ivec2     shape,
             glm::vec2      kw,
             std::uint32_t  seed,
             NoiseFunction &noise_function,
             float          noise_function_offset = 0.f,
             float          noise_function_scaling = 1.f,
             float          eps = 0.05,
             int            resolution = 1,
             float          displacement = 0.075,
             int            primitives_resolution_steps = 3,
             float          slope_power = 2.f,
             float          noise_amplitude_proportion = 0.01,
             bool           add_control_function = true,
             float          control_function_overlap = 0.5f,
             const Array   *p_noise_x = nullptr,
             const Array   *p_noise_y = nullptr,
             glm::vec4      bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a diffusion-limited aggregation (DLA) pattern.
 *
 * This function simulates the diffusion-limited aggregation process to generate
 * a pattern within a grid of specified dimensions. The DLA process models
 * particles that undergo a random walk until they stick to a seed, gradually
 * forming complex fractal structures.
 *
 * @param  shape                      The dimensions of the grid where the DLA
 *                                    pattern will be generated. It is
 *                                    represented as a `glm::ivec2` object,
 *                                    where the first element is the width and
 *                                    the second element is the height.
 * @param  scale                      A scaling factor that influences the
 *                                    density of the particles in the DLA
 *                                    pattern.
 * @param  seeding_radius             The radius within which initial seeding of
 *                                    particles occurs. This radius defines the
 *                                    area where the first particles are placed.
 * @param  seeding_outer_radius_ratio The ratio between the outer seeding radius
 *                                    and the initial seeding radius. It
 *                                    determines the outer boundary for particle
 *                                    seeding.
 * @param  slope                      Slope of the talus added to the DLA
 *                                    pattern.
 * @param  noise_ratio                A parameter that controls the amount of
 *                                    randomness or noise introduced in the
 *                                    talus formation process.
 * @param  seed                       The seed for the random number generator,
 *                                    ensuring reproducibility of the pattern.
 *                                    The same seed will generate the same
 *                                    pattern.
 *
 * @return                            A 2D array representing the generated DLA
 *                                    pattern. The array is of the same size as
 *                                    specified by `shape`.
 *
 * **Example**
 * @include ex_diffusion_limited_aggregation.cpp
 *
 * **Result**
 * @image html ex_diffusion_limited_aggregation.png
 */
Array diffusion_limited_aggregation(glm::ivec2    shape,
                                    float         scale,
                                    std::uint32_t seed,
                                    float         seeding_radius = 0.4f,
                                    float seeding_outer_radius_ratio = 0.2f,
                                    float slope = 8.f,
                                    float noise_ratio = 0.2f);

Array diffusion_limited_aggregation_trimesh(
    glm::ivec2            shape,
    std::uint32_t         seed,
    size_t                control_points_count = 5000,
    glm::vec2             seed_position = {0.5f, 0.5f},
    float                 ratio = 0.98f,
    float                 stop_proba = 1.f,
    float                 slope = 16.f,
    InterpolationMethod2D interpolation_method =
        InterpolationMethod2D::ITP2D_DELAUNAY_GRADIENT,
    const Array *p_noise_x = nullptr,
    const Array *p_noise_y = nullptr);

/**
 * @brief Return a sparse Gabor noise.
 *
 * @param  shape   Array shape.
 * @param  kw      Kernel wavenumber, with respect to a unit domain.
 * @param  angle   Kernel angle (in degree).
 * @param  width   Kernel width (in pixels).
 * @param  density Spot noise density.
 * @param  seed    Random seed number.
 * @return         Array New array.
 *
 * **Example**
 * @include ex_gabor_noise.cpp
 *
 * **Result**
 * @image html ex_gabor_noise.png
 */
Array gabor_noise(glm::ivec2    shape,
                  float         kw,
                  float         angle,
                  int           width,
                  float         density,
                  std::uint32_t seed);

/**
 * @brief Return an array filled with coherence noise.
 *
 * @param  noise_type           Noise type.
 * @param  shape                Array shape.
 * @param  kw                   Noise wavenumbers {kx, ky} for each directions.
 * @param  seed                 Random seed number.
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  bbox                 Domain bounding box.
 * @return                      Array Resulting array.
 *
 * **Example**
 * @include ex_noise.cpp
 *
 * **Result**
 * @image html ex_noise.png
 */
Array noise(NoiseType     noise_type,
            glm::ivec2    shape,
            glm::vec2     kw,
            std::uint32_t seed,
            const Array  *p_noise_x = nullptr,
            const Array  *p_noise_y = nullptr,
            glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return an array filled with coherence fbm noise.
 *
 * @param  noise_type           Noise type.
 * @param  shape                Array shape.
 * @param  kw                   Noise wavenumbers {kx, ky} for each directions.
 * @param  seed                 Random seed number.
 * @param  octaves              Number of octaves.
 * @param  weigth               Octave weighting.
 * @param  persistence          Octave persistence.
 * @param  lacunarity           Defines the wavenumber ratio between each
 *                              octaves.
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the weight parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  bbox                 Domain bounding box.
 * @return                      Array Fractal noise.
 *
 * **Example**
 * @include ex_noise_fbm.cpp
 *
 * **Result**
 * @image html ex_noise_fbm0.png
 * @image html ex_noise_fbm1.png
 * @image html ex_noise_fbm2.png
 * @image html ex_noise_fbm3.png
 * @image html ex_noise_fbm4.png
 * @image html ex_noise_fbm5.png
 * @image html ex_noise_fbm6.png
 */
Array noise_fbm(NoiseType     noise_type,
                glm::ivec2    shape,
                glm::vec2     kw,
                std::uint32_t seed,
                int           octaves = 8,
                float         weight = 0.7f,
                float         persistence = 0.5f,
                float         lacunarity = 2.f,
                const Array  *p_ctrl_param = nullptr,
                const Array  *p_noise_x = nullptr,
                const Array  *p_noise_y = nullptr,
                glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return an array filled with coherence fbm noise.
 *
 * @param  noise_type           Noise type.
 * @param  shape                Array shape.
 * @param  kw                   Noise wavenumbers {kx, ky} for each directions.
 * @param  seed                 Random seed number.
 * @param  octaves              Number of octaves.
 * @param  weigth               Octave weighting.
 * @param  persistence          Octave persistence.
 * @param  lacunarity           Defines the wavenumber ratio between each
 *                              octaves.
 * @param  gradient_scale       Gradient scale.
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the weight parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  bbox                 Domain bounding box.
 * @return                      Array Fractal noise.
 *
 * **Example**
 * @include ex_noise_fbm.cpp
 *
 * **Result**
 * @image html ex_noise_fbm0.png
 * @image html ex_noise_fbm1.png
 * @image html ex_noise_fbm2.png
 * @image html ex_noise_fbm3.png
 * @image html ex_noise_fbm4.png
 * @image html ex_noise_fbm5.png
 * @image html ex_noise_fbm6.png
 */
Array noise_iq(NoiseType     noise_type,
               glm::ivec2    shape,
               glm::vec2     kw,
               std::uint32_t seed,
               int           octaves = 8,
               float         weight = 0.7f,
               float         persistence = 0.5f,
               float         lacunarity = 2.f,
               float         gradient_scale = 0.05f,
               const Array  *p_ctrl_param = nullptr,
               const Array  *p_noise_x = nullptr,
               const Array  *p_noise_y = nullptr,
               glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return an array filled with coherence fbm noise.
 *
 * @param  noise_type           Noise type.
 * @param  shape                Array shape.
 * @param  kw                   Noise wavenumbers {kx, ky} for each directions.
 * @param  seed                 Random seed number.
 * @param  octaves              Number of octaves.
 * @param  weigth               Octave weighting.
 * @param  persistence          Octave persistence.
 * @param  lacunarity           Defines the wavenumber ratio between each
 *                              octaves.
 * @param  warp0                Initial warp scale.
 * @param  damp0                Initial damp scale.
 * @param  warp_scale           Warp scale.
 * @param  damp_scale           Damp scale.
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the weight parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  bbox                 Domain bounding box.
 * @return                      Array Fractal noise.
 *
 * **Example**
 * @include ex_noise_fbm.cpp
 *
 * **Result**
 * @image html ex_noise_fbm0.png
 * @image html ex_noise_fbm1.png
 * @image html ex_noise_fbm2.png
 * @image html ex_noise_fbm3.png
 * @image html ex_noise_fbm4.png
 * @image html ex_noise_fbm5.png
 * @image html ex_noise_fbm6.png
 */
Array noise_jordan(NoiseType     noise_type,
                   glm::ivec2    shape,
                   glm::vec2     kw,
                   std::uint32_t seed,
                   int           octaves = 8,
                   float         weight = 0.7f,
                   float         persistence = 0.5f,
                   float         lacunarity = 2.f,
                   float         warp0 = 0.4f,
                   float         damp0 = 1.f,
                   float         warp_scale = 0.4f,
                   float         damp_scale = 1.f,
                   const Array  *p_ctrl_param = nullptr,
                   const Array  *p_noise_x = nullptr,
                   const Array  *p_noise_y = nullptr,
                   glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return an array filled with coherent fbm Parberry variant of Perlin
 * noise.
 *
 * @param  shape                Array shape.
 * @param  kw                   Noise wavenumbers {kx, ky} for each directions.
 * @param  seed                 Random seed number.
 * @param  octaves              Number of octaves.
 * @param  weigth               Octave weighting.
 * @param  persistence          Octave persistence.
 * @param  lacunarity           Defines the wavenumber ratio between each
 *                              octaves.
 * @param  mu                   Gradient magnitude exponent.
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the weight parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  bbox                 Domain bounding box.
 * @return                      Array Fractal noise.
 *
 * **Example**
 * @include ex_noise_fbm.cpp
 *
 * **Result**
 * @image html ex_noise_fbm0.png
 * @image html ex_noise_fbm1.png
 * @image html ex_noise_fbm2.png
 * @image html ex_noise_fbm3.png
 * @image html ex_noise_fbm4.png
 * @image html ex_noise_fbm5.png
 * @image html ex_noise_fbm6.png
 */
Array noise_parberry(glm::ivec2    shape,
                     glm::vec2     kw,
                     std::uint32_t seed,
                     int           octaves = 8,
                     float         weight = 0.7f,
                     float         persistence = 0.5f,
                     float         lacunarity = 2.f,
                     float         mu = 1.02f,
                     const Array  *p_ctrl_param = nullptr,
                     const Array  *p_noise_x = nullptr,
                     const Array  *p_noise_y = nullptr,
                     glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return an array filled with coherence fbm pingpong noise.
 *
 * @param  noise_type           Noise type.
 * @param  shape                Array shape.
 * @param  kw                   Noise wavenumbers {kx, ky} for each directions.
 * @param  seed                 Random seed number.
 * @param  octaves              Number of octaves.
 * @param  weigth               Octave weighting.
 * @param  persistence          Octave persistence.
 * @param  lacunarity           Defines the wavenumber ratio between each
 *                              octaves.
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the weight parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  bbox                 Domain bounding box.
 * @return                      Array Fractal noise.
 *
 * **Example**
 * @include ex_noise_fbm.cpp
 *
 * **Result**
 * @image html ex_noise_fbm0.png
 * @image html ex_noise_fbm1.png
 * @image html ex_noise_fbm2.png
 * @image html ex_noise_fbm3.png
 * @image html ex_noise_fbm4.png
 * @image html ex_noise_fbm5.png
 * @image html ex_noise_fbm6.png
 */
Array noise_pingpong(NoiseType     noise_type,
                     glm::ivec2    shape,
                     glm::vec2     kw,
                     std::uint32_t seed,
                     int           octaves = 8,
                     float         weight = 0.7f,
                     float         persistence = 0.5f,
                     float         lacunarity = 2.f,
                     const Array  *p_ctrl_param = nullptr,
                     const Array  *p_noise_x = nullptr,
                     const Array  *p_noise_y = nullptr,
                     glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return an array filled with coherence fbm ridged noise.
 *
 * @param  noise_type           Noise type.
 * @param  shape                Array shape.
 * @param  kw                   Noise wavenumbers {kx, ky} for each directions.
 * @param  seed                 Random seed number.
 * @param  octaves              Number of octaves.
 * @param  weigth               Octave weighting.
 * @param  persistence          Octave persistence.
 * @param  lacunarity           Defines the wavenumber ratio between each
 *                              octaves.
 * @param  k_smoothing          Smoothing parameter.
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the weight parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  bbox                 Domain bounding box.
 * @return                      Array Fractal noise.
 *
 * **Example**
 * @include ex_noise_fbm.cpp
 *
 * **Result**
 * @image html ex_noise_fbm0.png
 * @image html ex_noise_fbm1.png
 * @image html ex_noise_fbm2.png
 * @image html ex_noise_fbm3.png
 * @image html ex_noise_fbm4.png
 * @image html ex_noise_fbm5.png
 * @image html ex_noise_fbm6.png
 */
Array noise_ridged(NoiseType     noise_type,
                   glm::ivec2    shape,
                   glm::vec2     kw,
                   std::uint32_t seed,
                   int           octaves = 8,
                   float         weight = 0.7f,
                   float         persistence = 0.5f,
                   float         lacunarity = 2.f,
                   float         k_smoothing = 0.1f,
                   const Array  *p_ctrl_param = nullptr,
                   const Array  *p_noise_x = nullptr,
                   const Array  *p_noise_y = nullptr,
                   glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return an array filled with coherence fbm swiss noise.
 *
 * @param  noise_type           Noise type.
 * @param  shape                Array shape.
 * @param  kw                   Noise wavenumbers {kx, ky} for each directions.
 * @param  seed                 Random seed number.
 * @param  octaves              Number of octaves.
 * @param  weigth               Octave weighting.
 * @param  persistence          Octave persistence.
 * @param  lacunarity           Defines the wavenumber ratio between each
 *                              octaves.
 * @param  warp_scale           Warp scale.
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the weight parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  bbox                 Domain bounding box.
 * @return                      Array Fractal noise.
 *
 * **Example**
 * @include ex_noise_fbm.cpp
 *
 * **Result**
 * @image html ex_noise_fbm0.png
 * @image html ex_noise_fbm1.png
 * @image html ex_noise_fbm2.png
 * @image html ex_noise_fbm3.png
 * @image html ex_noise_fbm4.png
 * @image html ex_noise_fbm5.png
 * @image html ex_noise_fbm6.png
 */
Array noise_swiss(NoiseType     noise_type,
                  glm::ivec2    shape,
                  glm::vec2     kw,
                  std::uint32_t seed,
                  int           octaves = 8,
                  float         weight = 0.7f,
                  float         persistence = 0.5f,
                  float         lacunarity = 2.f,
                  float         warp_scale = 0.1f,
                  const Array  *p_ctrl_param = nullptr,
                  const Array  *p_noise_x = nullptr,
                  const Array  *p_noise_y = nullptr,
                  glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return an array filled with the maximum of two Worley (cellular)
 * noises.
 *
 * @param  shape                Array shape.
 * @param  kw                   Noise wavenumbers {kx, ky} for each directions,
 *                              with respect to a unit domain.
 * @param  seed                 Random seed number.
 * @param  ratio                Amplitude ratio between each Worley noise.
 * @param  k                    Transition smoothing parameter.
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the ratio parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  p_stretching         Local wavenumber multiplier.
 * @param  bbox                 Domain bounding box.
 * @return                      Array Noise.
 *
 * **Example**
 * @include ex_worley_double.cpp
 *
 * **Result**
 * @image html ex_worley_double.png
 */
Array worley_double(glm::ivec2    shape,
                    glm::vec2     kw,
                    std::uint32_t seed,
                    float         ratio = 0.5f,
                    float         k = 0.f,
                    const Array  *p_ctrl_param = nullptr,
                    const Array  *p_noise_x = nullptr,
                    const Array  *p_noise_y = nullptr,
                    glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

} // namespace hmap

namespace hmap::gpu
{

/**
 * @brief Return an array filled with coherence Gabor noise.
 *
 * @param  shape              Array shape.
 * @param  kw                 Noise wavenumbers {kx, ky} for each directions.
 * @param  seed               Random seed number.
 * @param  angle              Base orientation angle for the Gabor wavelets (in
 *                            radians). Defaults to 0.
 * @param  angle_spread_ratio Ratio that controls the spread of wave
 *                            orientations around the base angle. Defaults to 1.
 * @param  bbox               Domain bounding box.
 * @return                    Array Fractal noise.
 *
 * @note Taken from https://www.shadertoy.com/view/clGyWm
 *
 * @note Only available if OpenCL is enabled.
 *
 * **Example**
 * @include ex_gabor_wave.cpp
 *
 * **Result**
 * @image html ex_gabor_wave.png
 */
Array gabor_wave(glm::ivec2    shape,
                 glm::vec2     kw,
                 std::uint32_t seed,
                 const Array  &angle,
                 float         angle_spread_ratio = 1.f,
                 glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

Array gabor_wave(glm::ivec2    shape,
                 glm::vec2     kw,
                 std::uint32_t seed,
                 float         angle = 0.f,
                 float         angle_spread_ratio = 1.f,
                 glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return an array filled with coherence Gabor noise.
 *
 * @param  shape                Array shape.
 * @param  kw                   Noise wavenumbers {kx, ky} for each directions.
 * @param  seed                 Random seed number.
 * @param  angle                Base orientation angle for the Gabor wavelets
 *                              (in radians). Defaults to 0.
 * @param  angle_spread_ratio   Ratio that controls the spread of wave
 *                              orientations around the base angle. Defaults to
 *                              1.
 * @param  octaves              Number of octaves.
 * @param  weigth               Octave weighting.
 * @param  persistence          Octave persistence.
 * @param  lacunarity           Defines the wavenumber ratio between each
 *                              octaves.
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the weight parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  bbox                 Domain bounding box.
 * @return                      Array Fractal noise.
 *
 * @note Taken from https://www.shadertoy.com/view/clGyWm
 *
 * @note Only available if OpenCL is enabled.
 *
 * **Example**
 * @include ex_gabor_wave.cpp
 *
 * **Result**
 * @image html ex_gabor_wave.png
 */
Array gabor_wave_fbm(glm::ivec2    shape,
                     glm::vec2     kw,
                     std::uint32_t seed,
                     const Array  &angle,
                     float         angle_spread_ratio = 1.f,
                     int           octaves = 8,
                     float         weight = 0.7f,
                     float         persistence = 0.5f,
                     float         lacunarity = 2.f,
                     const Array  *p_ctrl_param = nullptr,
                     const Array  *p_noise_x = nullptr,
                     const Array  *p_noise_y = nullptr,
                     glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

Array gabor_wave_fbm(glm::ivec2    shape,
                     glm::vec2     kw,
                     std::uint32_t seed,
                     float         angle = 0.f,
                     float         angle_spread_ratio = 1.f,
                     int           octaves = 8,
                     float         weight = 0.7f,
                     float         persistence = 0.5f,
                     float         lacunarity = 2.f,
                     const Array  *p_ctrl_param = nullptr,
                     const Array  *p_noise_x = nullptr,
                     const Array  *p_noise_y = nullptr,
                     glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a 2D array using the GavoroNoise algorithm, which is a
 * procedural noise technique for terrain generation and other applications.
 *
 * @param  shape           Dimensions of the output array.
 * @param  kw              Wave number vector controlling the noise frequency.
 * @param  seed            Seed value for random number generation.
 * @param  amplitude       Amplitude of the noise.
 * @param  kw_multiplier   Multiplier for wave numbers in the noise function.
 * @param  slope_strength  Strength of slope-based directional erosion in the
 *                         noise.
 * @param  branch_strength Strength of branch-like structures in the generated
 *                         noise.
 * @param  z_cut_min       Minimum cutoff for Z-value in the noise.
 * @param  z_cut_max       Maximum cutoff for Z-value in the noise.
 * @param  octaves         Number of octaves for fractal Brownian motion (fBm).
 * @param  persistence     Amplitude scaling factor between noise octaves.
 * @param  lacunarity      Frequency scaling factor between noise octaves.
 * @param  p_ctrl_param    Optional array for control parameters, can modify the
 *                         Z cutoff dynamically.
 * @param  p_noise_x       Optional array for X-axis noise perturbation.
 * @param  p_noise_y       Optional array for Y-axis noise perturbation.
 * @param  bbox            Bounding box for mapping grid coordinates to world
 *                         space.
 *
 * @return                 A 2D array containing the generated GavoroNoise
 *                         values.
 *
 * @note Taken from https://www.shadertoy.com/view/MtGcWh
 *
 * @note Only available if OpenCL is enabled.
 *
 * This function leverages an OpenCL kernel to compute the GavoroNoise values on
 * the GPU, allowing for efficient large-scale generation. The kernel applies a
 * combination of fractal Brownian motion (fBm), directional erosion, and other
 * procedural techniques to generate intricate noise patterns.
 *
 * The optional `p_ctrl_param`, `p_noise_x`, and `p_noise_y` buffers provide
 * additional flexibility for dynamically adjusting noise parameters and
 * perturbations.
 *
 * **Example**
 * @include ex_gavoronoise.cpp
 *
 * **Result**
 * @image html ex_gavoronoise.png
 */
Array gavoronoise(glm::ivec2    shape,
                  glm::vec2     kw,
                  std::uint32_t seed,
                  const Array  &angle,
                  float         amplitude = 0.05f,
                  float         angle_spread_ratio = 1.f,
                  glm::vec2     kw_multiplier = {4.f, 4.f},
                  float         slope_strength = 1.f,
                  float         branch_strength = 2.f,
                  float         z_cut_min = 0.2f,
                  float         z_cut_max = 1.f,
                  int           octaves = 8,
                  float         persistence = 0.4f,
                  float         lacunarity = 2.f,
                  const Array  *p_ctrl_param = nullptr,
                  const Array  *p_noise_x = nullptr,
                  const Array  *p_noise_y = nullptr,
                  glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

Array gavoronoise(glm::ivec2    shape,
                  glm::vec2     kw,
                  std::uint32_t seed,
                  float         angle = 0.f,
                  float         amplitude = 0.05f,
                  float         angle_spread_ratio = 1.f,
                  glm::vec2     kw_multiplier = {4.f, 4.f},
                  float         slope_strength = 1.f,
                  float         branch_strength = 2.f,
                  float         z_cut_min = 0.2f,
                  float         z_cut_max = 1.f,
                  int           octaves = 8,
                  float         persistence = 0.4f,
                  float         lacunarity = 2.f,
                  const Array  *p_ctrl_param = nullptr,
                  const Array  *p_noise_x = nullptr,
                  const Array  *p_noise_y = nullptr,
                  glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

Array gavoronoise(const Array  &base,
                  glm::vec2     kw,
                  std::uint32_t seed,
                  float         amplitude = 0.05f,
                  glm::vec2     kw_multiplier = {4.f, 4.f},
                  float         slope_strength = 1.f,
                  float         branch_strength = 2.f,
                  float         z_cut_min = 0.2f,
                  float         z_cut_max = 1.f,
                  int           octaves = 8,
                  float         persistence = 0.4f,
                  float         lacunarity = 2.f,
                  const Array  *p_ctrl_param = nullptr,
                  const Array  *p_noise_x = nullptr,
                  const Array  *p_noise_y = nullptr,
                  glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/*! @brief See hmap::noise. @p period sets the tiling period in lattice cells.
   For the heightmap to be periodic with respect to the full domain, @p period
   must be equal to the wavenumber @p kw, i.e. period = {round(kw.x),
   round(kw.y)}; a smaller period makes the pattern repeat within the domain. A
   component <= 0 disables wrapping on that axis (default {0, 0}:
 * non-periodic). Tiling is exact only for lattice noise types with integer kw;
 * simplex is never wrapped. */
Array noise(NoiseType     noise_type,
            glm::ivec2    shape,
            glm::vec2     kw,
            std::uint32_t seed,
            const Array  *p_noise_x = nullptr,
            const Array  *p_noise_y = nullptr,
            glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f},
            glm::ivec2    period = {0, 0});

/*! @brief See hmap::noise_fbm. @p period sets the tiling period in lattice
   cells for the base octave. For the heightmap to be periodic with respect to
   the full domain, @p period must be equal to the wavenumber @p kw, i.e.
 * period = {round(kw.x), round(kw.y)}; a smaller period makes the pattern
 * repeat within the domain. A component <= 0 disables wrapping on that axis
 * (default {0, 0}: non-periodic). Seamless tiling additionally requires an
 * integer lacunarity (the default 2); simplex is never wrapped. */
Array noise_fbm(NoiseType     noise_type,
                glm::ivec2    shape,
                glm::vec2     kw,
                std::uint32_t seed,
                int           octaves = 8,
                float         weight = 0.7f,
                float         persistence = 0.5f,
                float         lacunarity = 2.f,
                const Array  *p_ctrl_param = nullptr,
                const Array  *p_noise_x = nullptr,
                const Array  *p_noise_y = nullptr,
                glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f},
                glm::ivec2    period = {0, 0});

/**
 * @brief Generates a procedural phasor-based pattern.
 *
 * Produces a directional, phase-driven field using a selected profile, optional
 * angle control, and spatial distortion noise.
 *
 * @param  phasor_profile   Phasor shaping profile.
 * @param  shape            Output array dimensions.
 * @param  kp_global        Global frequency factor.
 * @param  seed             Random seed.
 * @param  angle_shift      Global phase angle offset (degrees).
 * @param  n_kernel_samples Number of kernel samples.
 * @param  jitter           Sampling jitter (x,y).
 * @param  delta            Finite difference step.
 * @param  phase_smoothing  Phase smoothing factor.
 * @param  p_angle          Optional external angle field.
 * @param  p_noise_x        Optional X distortion field.
 * @param  p_noise_y        Optional Y distortion field.
 * @param  bbox             Domain bounding box.
 * @return                  Array Generated phasor field.
 *
 * **Example**
 * @include ex_phasor.cpp
 *
 * **Result**
 * @image html ex_phasor.png
 */
Array phasor(PhasorProfile   phasor_profile,
             glm::ivec2      shape,
             float           kp_global,
             std::uint32_t   seed,
             float           angle_shift = 0.f,
             float           normalization = 1.f,
             const glm::vec2 jitter = {1.f, 1.f},
             float           delta = 0.01f,
             float           phase_smoothing = 10.f,
             const Array    *p_angle = nullptr,
             const Array    *p_noise_x = nullptr,
             const Array    *p_noise_y = nullptr,
             glm::vec4       bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a multi-octave (fBm) phasor-based pattern.
 *
 * Accumulates several scaled `phasor` layers using persistence and lacunarity
 * for richer multi-scale structure.
 *
 * @param  phasor_profile   Phasor shaping profile.
 * @param  shape            Output array dimensions.
 * @param  kp_global        Base frequency factor.
 * @param  seed             Random seed.
 * @param  angle_shift      Global phase angle offset (degrees).
 * @param  octaves          Number of fBm layers.
 * @param  weight           Initial octave weight.
 * @param  persistence      Amplitude multiplier per octave.
 * @param  lacunarity       Frequency multiplier per octave.
 * @param  normalization    Degree of normalization applied [0, 1].
 * @param  jitter           Sampling jitter (x,y).
 * @param  delta            Finite difference step.
 * @param  phase_smoothing  Phase smoothing factor.
 * @param  p_angle          Optional external angle field.
 * @param  p_noise_x        Optional X distortion field.
 * @param  p_noise_y        Optional Y distortion field.
 * @param  bbox             Domain bounding box.
 * @return                  Array Generated multi-scale phasor field.
 *
 * **Example**
 * @include ex_phasor.cpp
 *
 * **Result**
 * @image html ex_phasor.png
 */
Array phasor_fbm(PhasorProfile   phasor_profile,
                 glm::ivec2      shape,
                 float           kp_global,
                 std::uint32_t   seed,
                 float           angle_shift = 0.f,
                 int             octaves = 8,
                 float           weight = 0.7f,
                 float           persistence = 0.5f,
                 float           lacunarity = 2.f,
                 float           normalization = 1.f,
                 const glm::vec2 jitter = {1.f, 1.f},
                 float           delta = 0.01f,
                 float           phase_smoothing = 10.f,
                 const Array    *p_angle = nullptr,
                 const Array    *p_noise_x = nullptr,
                 const Array    *p_noise_y = nullptr,
                 glm::vec4       bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a Voronoi-based pattern where cells are defined by proximity
 * to random lines.
 *
 * This function generates an OpenCL-accelerated Voronoi-like pattern based on
 * the distance from each pixel to a set of randomly oriented lines. Each line
 * is defined by a random point and a direction sampled from a uniform
 * distribution around a given angle.
 *
 * @param  shape       The resolution of the resulting 2D array (width, height).
 * @param  density     Number of base points per unit area used to define lines.
 * @param  seed        Seed for the random number generator used to generate
 *                     base points and directions.
 * @param  k_smoothing Kernel smoothing factor; controls how sharp or soft the
 *                     distance fields are.
 * @param  exp_sigma   Exponential smoothing parameter applied to the computed
 *                     distance field.
 * @param  alpha       Base angle (in radians) used to orient the generated
 *                     lines.
 * @param  alpha_span  Maximum angular deviation from `alpha`; controls line
 *                     orientation variability.
 * @param  return_type Type of Voronoi output to return (e.g., F1, F2, edge
 *                     distance, smoothed field, etc.).
 * @param  p_noise_x   Optional pointer to an input noise field applied to the X
 *                     coordinates (can be nullptr).
 * @param  p_noise_y   Optional pointer to an input noise field applied to the Y
 *                     coordinates (can be nullptr).
 * @param  bbox        Bounding box in normalized coordinates (min_x, max_x,
 *                     min_y, max_y) of the final array.
 * @param  bbox_points Bounding box within which random base points are sampled.
 *
 * @return             A 2D array (of type Array) containing the computed
 *                     distance field based on line proximity.
 *
 * @note Each line is defined from a point (x, y) to a direction offset using
 * angle `theta = alpha + rand * alpha_span`.
 * @note The OpenCL kernel "vorolines" must be defined and compiled beforehand.
 *
 * **Example**
 * @include ex_vorolines.cpp
 *
 * **Result**
 * @image html ex_vorolines.png
 * @image html ex_vorolines_fbm.png
 */
Array vorolines(glm::ivec2        shape,
                float             density,
                std::uint32_t     seed,
                float             k_smoothing = 0.f,
                float             exp_sigma = 0.f,
                float             alpha = 0.f,
                float             alpha_span = M_PI,
                VoronoiReturnType return_type = VoronoiReturnType::F1_SQUARED,
                const Array      *p_noise_x = nullptr,
                const Array      *p_noise_y = nullptr,
                glm::vec4         bbox = {0.f, 1.f, 0.f, 1.f},
                glm::vec4         bbox_points = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a Voronoi-based pattern using distances to lines defined by
 * random points and angles, with additional fractal Brownian motion (fBm) noise
 * modulation.
 *
 * This function extends the standard `vorolines` generation by introducing
 * fBm-based warping of the coordinate space, resulting in more organic and
 * fractal-like structures. It creates a Voronoi distance field based on
 * proximity to oriented line segments and distorts the result using
 * multi-octave procedural noise.
 *
 * @param  shape       Output resolution of the 2D array (width, height).
 * @param  density     Number of base points per unit area used to define lines.
 * @param  seed        Seed value for the random number generator.
 * @param  k_smoothing Kernel smoothing coefficient to soften distance values
 *                     (e.g., for blending).
 * @param  exp_sigma   Sigma value for optional exponential smoothing on the
 *                     final field.
 * @param  alpha       Base orientation angle (in radians) of lines generated
 *                     from random points.
 * @param  alpha_span  Maximum angle deviation from `alpha`, determining
 *                     directional randomness of lines.
 * @param  return_type Type of output to return (e.g., F1, F2, distance to edge,
 *                     smoothed version).
 * @param  octaves     Number of noise octaves used in the fBm modulation.
 * @param  weight      Weight of each octave's contribution to the total noise.
 * @param  persistence Amplitude decay factor for each successive octave
 *                     (commonly 0.5–0.8).
 * @param  lacunarity  Frequency multiplier for each successive octave (commonly
 *                     2.0).
 * @param  p_noise_x   Optional pointer to an external noise field applied to X
 *                     coordinates (can be nullptr).
 * @param  p_noise_y   Optional pointer to an external noise field applied to Y
 *                     coordinates (can be nullptr).
 * @param  bbox        Bounding box for the final image domain (min_x, max_x,
 *                     min_y, max_y).
 * @param  bbox_points Bounding box from which the initial set of points are
 *                     sampled.
 *
 * @return             A 2D `Array` representing the Voronoi-fBm field,
 *                     distorted by noise and influenced by distance to random
 *                     lines.
 *
 * @note This version uses internally computed fBm noise unless external fields
 * (`p_noise_x`, `p_noise_y`) are provided.
 * @note This function requires an OpenCL kernel named "vorolines_fbm" to be
 * compiled and accessible.
 *
 * **Example**
 * @include ex_vorolines.cpp
 *
 * **Result**
 * @image html ex_vorolines.png
 * @image html ex_vorolines_fbm.png
 */
Array vorolines_fbm(
    glm::ivec2        shape,
    float             density,
    std::uint32_t     seed,
    float             k_smoothing = 0.f,
    float             exp_sigma = 0.f,
    float             alpha = 0.f,
    float             alpha_span = M_PI,
    VoronoiReturnType return_type = VoronoiReturnType::F1_SQUARED,
    int               octaves = 8,
    float             weight = 0.7f,
    float             persistence = 0.5f,
    float             lacunarity = 2.f,
    const Array      *p_noise_x = nullptr,
    const Array      *p_noise_y = nullptr,
    glm::vec4         bbox = {0.f, 1.f, 0.f, 1.f},
    glm::vec4         bbox_points = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a Voronoi diagram in a 2D array with configurable
 * properties.
 *
 * @param  shape        The dimensions of the output array as a 2D vector of
 *                      integers.
 * @param  kw           The frequency scale factors for the Voronoi cells, given
 *                      as a 2D vector of floats.
 * @param  seed         A seed value for random number generation, ensuring
 *                      reproducibility.
 * @param  jitter       (Optional) The amount of random variation in the
 *                      positions of Voronoi cell sites, given as a 2D vector of
 *                      floats. Defaults to {0.5f, 0.5f}.
 * @param  return_type  (Optional) The type of value to compute for the Voronoi
 *                      diagram. Defaults to `VoronoiReturnType::F1_SQUARED`.
 * @param  p_ctrl_param (Optional) A pointer to an `Array` used to control the
 *                      Voronoi computation. Used here as a multiplier for the
 *                      jitter. If nullptr, no control is applied.
 * @param  p_noise_x    (Optional) A pointer to an `Array` providing additional
 *                      noise in the x-direction for cell positions. If nullptr,
 *                      no x-noise is applied.
 * @param  p_noise_y    (Optional) A pointer to an `Array` providing additional
 *                      noise in the y-direction for cell positions. If nullptr,
 *                      no y-noise is applied.
 * @param  bbox         (Optional) The bounding box for the Voronoi computation,
 *                      given as a 4D vector of floats representing {min_x,
 *                      max_x, min_y, max_y}. Defaults to
 * {0.f, 1.f, 0.f, 1.f}.
 *
 * @return              A 2D array representing the generated Voronoi diagram.
 *
 * @note Only available if OpenCL is enabled.
 *
 * **Example**
 * @include ex_voronoi.cpp
 *
 * **Result**
 * @image html ex_voronoi.png
 */
Array voronoi(glm::ivec2        shape,
              glm::vec2         kw,
              std::uint32_t     seed,
              glm::vec2         jitter = {0.5f, 0.5f},
              float             k_smoothing = 0.f,
              float             exp_sigma = 0.f,
              VoronoiReturnType return_type = VoronoiReturnType::F1_SQUARED,
              const Array      *p_ctrl_param = nullptr,
              const Array      *p_noise_x = nullptr,
              const Array      *p_noise_y = nullptr,
              glm::vec4         bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a Voronoi diagram in a 2D array with configurable
 * properties.
 *
 * @param  shape        The dimensions of the output array as a 2D vector of
 *                      integers.
 * @param  kw           The frequency scale factors for the base Voronoi cells,
 *                      given as a 2D vector of floats.
 * @param  seed         A seed value for random number generation, ensuring
 *                      reproducibility.
 * @param  jitter       (Optional) The amount of random variation in the
 *                      positions of Voronoi cell sites, given as a 2D vector of
 *                      floats. Defaults to {0.5f, 0.5f}.
 * @param  return_type  (Optional) The type of value to compute for the Voronoi
 *                      diagram. Defaults to `VoronoiReturnType::F1_SQUARED`.
 * @param  octaves      (Optional) The number of layers (octaves) in the fractal
 *                      Brownian motion. Defaults to 8.
 * @param  weight       (Optional) The initial weight of the base layer in the
 *                      FBM computation. Defaults to 0.7f.
 * @param  persistence  (Optional) The persistence factor that controls the
 *                      amplitude reduction between octaves. Defaults to 0.5f.
 * @param  lacunarity   (Optional) The lacunarity factor that controls the
 *                      frequency increase between octaves. Defaults to 2.f.
 * @param  p_ctrl_param (Optional) A pointer to an `Array` used to control the
 *                      Voronoi computation. If nullptr, no control is applied.
 * @param  p_noise_x    (Optional) A pointer to an `Array` providing additional
 *                      noise in the x-direction for cell positions. If nullptr,
 *                      no x-noise is applied.
 * @param  p_noise_y    (Optional) A pointer to an `Array` providing additional
 *                      noise in the y-direction for cell positions. If nullptr,
 *                      no y-noise is applied.
 * @param  bbox         (Optional) The bounding box for the Voronoi computation,
 *                      given as a 4D vector of floats representing {min_x,
 *                      max_x, min_y, max_y}. Defaults to
 * {0.f, 1.f, 0.f, 1.f}.
 *
 * @return              A 2D array representing the generated Voronoi diagram.
 *
 * @note Only available if OpenCL is enabled.
 *
 * **Example**
 * @include ex_voronoi.cpp
 *
 * **Result**
 * @image html ex_voronoi.png
 */
Array voronoi_fbm(glm::ivec2        shape,
                  glm::vec2         kw,
                  std::uint32_t     seed,
                  glm::vec2         jitter = {0.5f, 0.5f},
                  float             k_smoothing = 0.f,
                  float             exp_sigma = 0.f,
                  VoronoiReturnType return_type = VoronoiReturnType::F1_SQUARED,
                  int               octaves = 8,
                  float             weight = 0.7f,
                  float             persistence = 0.5f,
                  float             lacunarity = 2.f,
                  const Array      *p_ctrl_param = nullptr,
                  const Array      *p_noise_x = nullptr,
                  const Array      *p_noise_y = nullptr,
                  glm::vec4         bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Computes the Voronoi edge distance.
 *
 * @param shape                The shape of the grid as a 2D vector (width,
 *                             height).
 * @param kw                   The weights for the Voronoi kernel as a 2D
 *                             vector.
 * @param seed                 The random seed used for generating Voronoi
 *                             points.
 * @param jitter               Optional parameter for controlling jitter in
 *                             Voronoi point placement (default is {0.5f,
 *                             0.5f}).
 * @param p_ctrl_param         Optional pointer to an Array specifying control
 *                             parameters for Voronoi grid jitter (default is
 *                             nullptr).
 * @param p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param bbox                 The bounding box for the Voronoi diagram as
 *                            {x_min, x_max, y_min, y_max} (default is {0.f,
 * 1.f, 0.f, 1.f}).
 *
 * @note Taken from https://www.shadertoy.com/view/llG3zy
 *
 * @note Only available if OpenCL is enabled.
 *
 * @note The resulting Array has the same dimensions as the input shape.
 */
Array voronoi_edge_distance(glm::ivec2    shape,
                            glm::vec2     kw,
                            std::uint32_t seed,
                            glm::vec2     jitter = {0.5f, 0.5f},
                            const Array  *p_ctrl_param = nullptr,
                            const Array  *p_noise_x = nullptr,
                            const Array  *p_noise_y = nullptr,
                            glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a 2D Voronoi noise array.
 *
 * This function computes a Voronoi noise pattern based on the input parameters
 * and returns it as a 2D array. The noise is calculated in the OpenCL kernel
 * `noise_voronoise`, which uses a combination of hashing and smoothstep
 * functions to generate a weighted Voronoi noise field.
 *
 * @note Taken from https://www.shadertoy.com/view/Xd23Dh
 *
 * @note Only available if OpenCL is enabled.
 *
 * @param  shape                The dimensions of the 2D output array as a
 *                              vector (width and height).
 * @param  kw                   Wave numbers for scaling the noise pattern,
 *                              represented as a 2D vector.
 * @param  u_param              A control parameter for the noise, adjusting the
 *                              contribution of random offsets.
 * @param  v_param              A control parameter for the noise, affecting the
 *                              smoothness of the pattern.
 * @param  p_ctrl_param         Optional pointer to an Array specifying control
 *                              parameters for Voronoi grid jitter (default is
 *                              nullptr).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  seed                 A seed value for random number generation,
 *                              ensuring reproducibility.
 *
 * @return                      An `Array` object containing the generated 2D
 *                              Voronoi noise values.
 *
 * **Example**
 * @include ex_voronoise.cpp
 *
 * **Result**
 * @image html ex_voronoise.png
 */
Array voronoise(glm::ivec2    shape,
                glm::vec2     kw,
                float         u_param,
                float         v_param,
                std::uint32_t seed,
                const Array  *p_noise_x = nullptr,
                const Array  *p_noise_y = nullptr,
                glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Return an array filled with coherence Voronoise.
 *
 * @param  shape                Array shape.
 * @param  kw                   Noise wavenumbers {kx, ky} for each directions.
 * @param  seed                 Random seed number.
 * @param  octaves              Number of octaves.
 * @param  weigth               Octave weighting.
 * @param  persistence          Octave persistence.
 * @param  lacunarity           Defines the wavenumber ratio between each
 *                              octaves.
 * @param  p_ctrl_param         Reference to the control parameter array (acts
 *                              as a multiplier for the weight parameter).
 * @param  p_noise_x, p_noise_y Reference to the input noise arrays.
 * @param  bbox                 Domain bounding box.
 * @return                      Array Fractal noise.
 *
 * @note Taken from https://www.shadertoy.com/view/clGyWm
 *
 * @note Only available if OpenCL is enabled.
 *
 * **Example**
 * @include ex_voronoise.cpp
 *
 * **Result**
 * @image html ex_voronoise.png
 */
Array voronoise_fbm(glm::ivec2    shape,
                    glm::vec2     kw,
                    float         u_param,
                    float         v_param,
                    std::uint32_t seed,
                    int           octaves = 8,
                    float         weight = 0.7f,
                    float         persistence = 0.5f,
                    float         lacunarity = 2.f,
                    const Array  *p_ctrl_param = nullptr,
                    const Array  *p_noise_x = nullptr,
                    const Array  *p_noise_y = nullptr,
                    glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a 2D Voronoi-based scalar field using OpenCL.
 *
 * This function computes a Voronoi diagram or derived metric (such as F1, F2,
 * or edge distances) on a grid of given shape. A set of random points is
 * generated within an extended bounding box, based on the desired density and
 * variability, to reduce edge artifacts. Optionally, per-pixel displacement can
 * be applied through noise fields.
 *
 * The result is stored in an `Array` object representing a 2D scalar field.
 *
 * @param  shape       Dimensions of the output array (width x height).
 * @param  density     Number of random points per unit area for Voronoi
 *                     diagram.
 * @param  variability Amount of randomness added to the point generation
 *                     bounding box.
 * @param  seed        Seed for random number generation used in point sampling.
 * @param  k_smoothing Smoothing factor used in soft minimum/maximum Voronoi
 *                     distance computations.
 * @param  exp_sigma   Standard deviation used in exponential falloff for edge
 *                     distance computation.
 * @param  return_type Type of Voronoi computation to perform (e.g., F1, F2,
 *                     F2-F1, edge distance).
 * @param  p_noise_x   Optional pointer to a noise field applied to X
 *                     coordinates of grid points.
 * @param  p_noise_y   Optional pointer to a noise field applied to Y
 *                     coordinates of grid points.
 * @param  bbox        Bounding box of the domain in which the field is
 *                     computed: {xmin, xmax, ymin, ymax}.
 * @param  bbox_points Bounding box for point generation, usually larger than
 * `bbox` to avoid edge effects.
 *
 * @return             An `Array` object containing the computed scalar field.
 *
 * @note
 * - The kernel `"vororand"` must be compiled and available in the OpenCL
 * context.
 * - If `p_noise_x` or `p_noise_y` are provided, they must match the shape of
 * the output array.
 * - The generated point cloud will be larger than `bbox` to reduce border
 * artifacts.
 *
 * **Example**
 * @include ex_vororand.cpp
 *
 * **Result**
 * @image html ex_vororand.png
 */
Array vororand(glm::ivec2        shape,
               float             density,
               float             variability,
               std::uint32_t     seed,
               float             k_smoothing = 0.f,
               float             exp_sigma = 0.f,
               VoronoiReturnType return_type = VoronoiReturnType::F1_SQUARED,
               const Array      *p_noise_x = nullptr,
               const Array      *p_noise_y = nullptr,
               glm::vec4         bbox = {0.f, 1.f, 0.f, 1.f},
               glm::vec4         bbox_points = {0.f, 1.f, 0.f, 1.f});

Array vororand(glm::ivec2                shape,
               const std::vector<float> &xp,
               const std::vector<float> &yp,
               float                     k_smoothing = 0.f,
               float                     exp_sigma = 0.f,
               VoronoiReturnType return_type = VoronoiReturnType::F1_SQUARED,
               const Array      *p_noise_x = nullptr,
               const Array      *p_noise_y = nullptr,
               glm::vec4         bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates 2D wavelet noise using an OpenCL kernel.
 *
 * This function creates a 2D noise field of the given @p shape using
 * multi-octave wavelet turbulence. The noise is computed on the GPU via the
 * "wavelet_noise" OpenCL kernel. Various parameters control the frequency
 * evolution, amplitude decay, vorticity injection, and optional modulation
 * using additional input arrays.
 *
 * @param  shape         Resolution of the output noise array (width, height).
 * @param  kw            Base wave number (frequency) in each axis. Controls the
 *                       initial spatial frequency of the wavelet field.
 * @param  seed          Integer seed used by the random hashing functions
 *                       inside the kernel.
 * @param  kw_multiplier Multiplier applied to the wave number per octave
 *                       (similar to gain). Higher values increase frequency
 *                       growth across octaves.
 * @param  vorticity     Amount of rotational distortion injected into the
 *                       noise. When > 0, the kernel applies additional angular
 *                       perturbations to produce turbulent flow-like patterns.
 * @param  density       Global scaling factor controlling the overall contrast
 *                       or amplitude of the generated noise.
 * @param  octaves       Number of wavelet noise layers combined. Higher octaves
 *                       yield more detail but increase computational cost.
 * @param  weight        Weight applied to the base octave, used as a global
 *                       intensity multiplier before persistence attenuation is
 *                       applied.
 * @param  persistence   Amplitude decay factor per octave. Values in (0,1)
 *                       produce the classic fractal noise falloff; higher
 *                       values retain more high-frequency detail.
 * @param  lacunarity    Frequency multiplier per octave. Controls how rapidly
 *                       frequency increases at each octave, with typical values
 *                       in the range [1.5, 3.0].
 * @param  p_ctrl_param  Optional pointer to a control-parameter array. When
 *                       provided, values inside this array modulate the
 *                       generated noise spatially. Pass nullptr to disable this
 *                       feature.
 * @param  p_noise_x     Optional pointer to an auxiliary noise array used to
 *                       perturb sampling positions horizontally. Pass nullptr
 *                       to disable.
 * @param  p_noise_y     Optional pointer to an auxiliary noise array used to
 *                       perturb sampling positions vertically. Pass nullptr to
 *                       disable.
 * @param  bbox          Bounding-box mapping (xmin, ymin, xmax, ymax). Converts
 *                       pixel-space indices into world-space coordinates for
 *                       spatially consistent noise.
 *
 * @return               Array A 2D floating-point array containing the
 *                       generated wavelet noise.
 *
 * **Example**
 * @include ex_wavelet_noise.cpp
 *
 * **Result**
 * @image html ex_wavelet_noise.png
 */
Array wavelet_noise(glm::ivec2    shape,
                    glm::vec2     kw,
                    std::uint32_t seed,
                    float         kw_multiplier = 2.f,
                    float         vorticity = 0.f,
                    float         density = 1.f,
                    int           octaves = 8,
                    float         weight = 0.7f,
                    float         persistence = 0.5f,
                    float         lacunarity = 2.f,
                    const Array  *p_ctrl_param = nullptr,
                    const Array  *p_noise_x = nullptr,
                    const Array  *p_noise_y = nullptr,
                    glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

} // namespace hmap::gpu
