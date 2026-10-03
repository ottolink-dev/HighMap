/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once

#include <cstdint>

#include "highmap/array.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <cfloat>

namespace hmap::gpu
{

/**
 * @brief Simulate a mudslide (landslide-driven material redistribution) on a
 * height field.
 *
 * @param[in,out] z                   Input/output height field modified in
 *                                    place.
 * @param[in]     landslide_mask      Mask defining affected areas (non-zero =
 *                                    active).
 * @param[in]     depth               Maximum erosion/deposition depth.
 * @param[in]     iterations          Number of simulation iterations.
 * @param[in]     depth_map_exponent  Exponent applied to depth influence
 *                                    (default 0.5f).
 * @param[in]     viscosity_law_power Power exponent controlling viscosity
 *                                    response (default 1.5f).
 * @param[out]    p_depth_end         Optional pointer to store resulting depth
 *                                    map.
 * @param[out]    p_depth_init        Optional pointer to store initial depth
 *                                    map.
 *
 * **Example**
 * @include ex_mudslide.cpp
 *
 * **Result**
 * @image html ex_mudslide.png
 */
void mudslide(Array       &z,
              const Array &landslide_mask,
              float        depth,
              int          iterations,
              float        depth_map_exponent = 0.5f,
              float        viscosity_law_power = 1.5f,
              Array       *p_depth_end = nullptr,
              Array       *p_depth_init = nullptr);

/*! @brief See hmap::gpu::mudslide */
void mudslide(Array &z,
              float  talus_limit,
              float  depth,
              int    iterations,
              float  depth_map_exponent = 0.5f,
              float  viscosity_law_power = 1.5f,
              Array *p_depth_end = nullptr,
              Array *p_depth_init = nullptr);

/**
 * @brief Applies a "rift" deformation effect to a heightmap array.
 *
 * @param z                     Reference to the heightmap array to be modified
 *                              in-place.
 * @param kw                    Frequency vector (kx, ky) scaling deformation.
 * @param angle                 Orientation of rift in degrees.
 * @param amplitude             Strength of rift deformation.
 * @param seed                  Random seed.
 * @param elevation_noise_shift Vertical offset applied to base noise.
 * @param k_smooth_bottom       Lower smoothing factor for Voronoi noise.
 * @param k_smooth_top          Upper smoothing factor for Voronoi noise.
 * @param radial_spread_amp     Amplitude controlling radial spreading.
 * @param elevation_noise_amp   Amplitude scaling elevation influence.
 * @param clamp_vmin            Minimum clamp value.
 * @param remap_vmin            Minimum remap value.
 * @param apply_mask            If true, applies power-based blending mask.
 * @param reverse_mask          Reverse mask blending.
 * @param mask_gamma            Gamma exponent for blending mask.
 * @param p_noise_x             Optional X noise array.
 * @param p_noise_y             Optional Y noise array.
 * @param p_mask                Optional mask array.
 * @param center                2D center point.
 * @param bbox                  Bounding box (xmin, xmax, ymin, ymax).
 *
 * **Example**
 * @include ex_rifts.cpp
 *
 * **Result**
 * @image html ex_rifts.png
 */
void rifts(Array           &z,
           const glm::vec2 &kw,
           float            angle,
           float            amplitude,
           std::uint32_t    seed,
           float            elevation_noise_shift = 0.f,
           float            k_smooth_bottom = 0.05f,
           float            k_smooth_top = 0.05f,
           float            radial_spread_amp = 0.2f,
           float            elevation_noise_amp = 0.1f,
           float            clamp_vmin = 0.f,
           float            remap_vmin = 0.f,
           bool             apply_mask = true,
           bool             reverse_mask = false,
           float            mask_gamma = 1.f,
           const Array     *p_noise_x = nullptr,
           const Array     *p_noise_y = nullptr,
           const Array     *p_mask = nullptr,
           const glm::vec2 &center = {0.5f, 0.5f},
           const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Perform sediment deposition combined with thermal erosion.
 *
 * @param z                     Input array.
 * @param p_mask                Intensity mask, expected in [0, 1].
 * @param talus                 Talus limit.
 * @param p_deposition_map      Reference to deposition map output.
 * @param max_deposition        Maximum height of sediment deposition.
 * @param iterations            Number of iterations.
 * @param thermal_subiterations Number of thermal erosion iterations per pass.
 *
 * **Example**
 * @include ex_sediment_deposition.cpp
 *
 * **Result**
 * @image html ex_sediment_deposition.png
 */
void sediment_deposition(Array       &z,
                         const Array *p_mask,
                         const Array &talus,
                         Array       *p_deposition_map = nullptr,
                         float        max_deposition = 0.01,
                         int          iterations = 5,
                         int          thermal_subiterations = 10);

void sediment_deposition(Array       &z,
                         const Array &talus,
                         Array       *p_deposition_map = nullptr,
                         float        max_deposition = 0.01,
                         int          iterations = 5,
                         int          thermal_subiterations = 10);

/**
 * @brief Applies a talus-based sediment deposition layer.
 *
 * @param[in,out] z                 Heightmap to modify.
 * @param[in]     talus_layer       Base talus values for erosion.
 * @param[in]     talus_upper_limit Upper talus threshold used for masking.
 * @param[in]     iterations        Number of erosion iterations.
 * @param[in]     apply_post_filter Enable post smoothing and blending.
 * @param[out]    p_deposition_map  Optional output deposition map.
 *
 * **Example**
 * @include ex_sediment_layer.cpp
 *
 * **Result**
 * @image html ex_sediment_layer.png
 */
void sediment_layer(Array       &z,
                    const Array &talus_layer,
                    const Array &talus_upper_limit,
                    int          iterations,
                    bool         apply_post_filter = true,
                    Array       *p_deposition_map = nullptr);

/**
 * @brief Applies stratification to a heightfield using directional noise and
 * multiscale gamma transformations.
 *
 * @param z                    Reference to heightfield array, MUST BE
 *                             NORMALIZED in [0, 1].
 * @param angle                Horizontal orientation of strata in degrees.
 * @param slope                Vertical slope of strata.
 * @param gamma                Gamma exponent for non-linear remapping.
 * @param seed                 Seed for deterministic noise generation.
 * @param linear_gamma         Linear or smooth gamma mapping.
 * @param kz                   Base scaling factor for frequency.
 * @param octaves              Number of iterative passes.
 * @param lacunarity           Frequency multiplier per octave.
 * @param gamma_noise_ratio    Noise influence ratio on gamma.
 * @param noise_amp            Amplitude of base Perlin noise.
 * @param noise_kw             Frequency vector for Perlin noise.
 * @param enable_ridge_noise   Enable Voronoi ridge noise.
 * @param ridge_noise_kw       Frequency vector for ridge noise.
 * @param ridge_angle_shift    Angular shift for ridge direction.
 * @param ridge_noise_amp      Amplitude of ridge noise.
 * @param ridge_clamp_vmin     Minimum clamp value for ridge noise.
 * @param ridge_remap_vmin     Minimum remap value for ridge noise.
 * @param apply_elevation_mask Apply elevation mask.
 * @param apply_ridge_mask     Apply ridge mask.
 * @param mask_gamma           Mask gamma exponent.
 * @param p_mask               Optional filter mask.
 * @param bbox                 Domain bounding box.
 *
 * **Example**
 * @include ex_strata.cpp
 *
 * **Result**
 * @image html ex_strata.png
 */
void strata(Array           &z,
            float            angle,
            float            slope,
            float            gamma,
            std::uint32_t    seed,
            bool             linear_gamma = true,
            float            kz = 1.f,
            int              octaves = 4,
            float            lacunarity = 2.f,
            float            gamma_noise_ratio = 0.5f,
            float            noise_amp = 0.4f,
            const glm::vec2 &noise_kw = {4.f, 4.f},
            bool             enable_ridge_noise = true,
            const glm::vec2 &ridge_noise_kw = {4.f, 1.2f},
            float            ridge_angle_shift = 45.f,
            float            ridge_noise_amp = 0.5f,
            float            ridge_clamp_vmin = 0.f,
            float            ridge_remap_vmin = 0.f,
            bool             apply_elevation_mask = true,
            bool             apply_ridge_mask = true,
            float            mask_gamma = 0.4f,
            const Array     *p_mask = nullptr,
            const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Applies procedural stratified cell displacement to a heightmap.
 *
 * **Example**
 * @include ex_strata_cells.cpp
 *
 * **Result**
 * @image html ex_strata_cells.png
 */
void strata_cells(Array        &z,
                  glm::vec2     kw,
                  float         amp,
                  std::uint32_t seed,
                  float         gamma = 0.5f,
                  float         gamma_lateral = 0.4f,
                  float         angle = 0.f,
                  float         noise_amp = 0.5f,
                  bool          absolute_displacement = true,
                  float         occurence_probability = 0.5f,
                  const Array  *p_noise_x = nullptr,
                  const Array  *p_noise_y = nullptr,
                  glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

void strata_cells(Array        &z,
                  glm::vec2     kw,
                  float         amp,
                  std::uint32_t seed,
                  const Array  *p_mask,
                  float         gamma = 0.5f,
                  float         gamma_lateral = 0.4f,
                  float         angle = 0.f,
                  float         noise_amp = 0.5f,
                  bool          absolute_displacement = true,
                  float         occurence_probability = 0.5f,
                  const Array  *p_noise_x = nullptr,
                  const Array  *p_noise_y = nullptr,
                  glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Applies multi-octave (fBm) stratified cell displacement.
 *
 * **Example**
 * @include ex_strata_cells.cpp
 *
 * **Result**
 * @image html ex_strata_cells.png
 */
void strata_cells_fbm(Array        &z,
                      glm::vec2     kw,
                      float         amp,
                      std::uint32_t seed,
                      float         gamma = 0.5f,
                      float         gamma_lateral = 0.4f,
                      float         angle = 0.f,
                      bool          enable_default_noise = true,
                      float         default_noise_amp = 0.05f,
                      bool          absolute_displacement = true,
                      float         occurence_probability = 0.5f,
                      int           octaves = 8,
                      float         persistence = 0.4f,
                      float         lacunarity = 2.2f,
                      const Array  *p_noise_x = nullptr,
                      const Array  *p_noise_y = nullptr,
                      glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

void strata_cells_fbm(Array        &z,
                      glm::vec2     kw,
                      float         amp,
                      std::uint32_t seed,
                      const Array  *p_mask,
                      float         gamma = 0.5f,
                      float         gamma_lateral = 0.4f,
                      float         angle = 0.f,
                      bool          enable_default_noise = true,
                      float         default_noise_amp = 0.05f,
                      bool          absolute_displacement = true,
                      float         occurence_probability = 0.5f,
                      int           octaves = 8,
                      float         persistence = 0.4f,
                      float         lacunarity = 2.2f,
                      const Array  *p_noise_x = nullptr,
                      const Array  *p_noise_y = nullptr,
                      glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Apply stratified talus projection along multiple directions.
 *
 * **Example**
 * @include ex_strata_plates.cpp
 *
 * **Result**
 * @image html ex_strata_plates.png
 */
void strata_plates(Array        &z,
                   const Array  &talus,
                   int           direction_offset = 0,
                   int           direction_count = 3,
                   bool          random_directions = false,
                   std::uint32_t seed = 0,
                   float         vmin = -FLT_MAX,
                   float         skew = 0.f,
                   float         mix_ratio = 0.9f,
                   const Array  *p_mask = nullptr,
                   const Array  *p_dx = nullptr,
                   const Array  *p_dy = nullptr);

/**
 * @brief Applies a terrace (stratification) filter to a heightmap.
 *
 * **Example**
 * @include ex_strata_terrace.cpp
 *
 * **Result**
 * @image html ex_strata_terrace.png
 */
void strata_terrace(Array        &z,
                    float         gamma,
                    std::uint32_t seed,
                    float         kz = 4.f,
                    bool          linear_gamma = true,
                    float         gamma_noise_ratio = 0.5f,
                    float         slope = 0.f,
                    float         angle = 0.f,
                    const Array  *p_noise = nullptr,
                    glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

void strata_terrace(Array        &z,
                    float         gamma,
                    std::uint32_t seed,
                    const Array  *p_mask,
                    float         kz = 4.f,
                    bool          linear_gamma = true,
                    float         gamma_noise_ratio = 0.5f,
                    float         slope = 0.f,
                    float         angle = 0.f,
                    const Array  *p_noise = nullptr,
                    glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Fill valleys using thermal scree deposition and height-based blending.
 *
 * **Example**
 * @include ex_valley_fill.cpp
 *
 * **Result**
 * @image html ex_valley_fill.png
 */
void valley_fill(Array       &z,
                 const Array &talus,
                 int          iterations = 100,
                 float        gamma = 2.f,
                 float        ratio = 0.8f,
                 float        zmin = 0.f,
                 float        zmax = 0.f,
                 float        elevation_max_ratio = 1.f,
                 bool         preserve_elevation_range = true,
                 const Array *p_noise = nullptr,
                 Array       *p_deposition_map = nullptr);

void valley_fill(Array       &z,
                 const Array *p_mask,
                 const Array &talus,
                 int          iterations = 100,
                 float        gamma = 2.f,
                 float        ratio = 0.8f,
                 float        zmin = 0.f,
                 float        zmax = 0.f,
                 float        elevation_max_ratio = 1.f,
                 bool         preserve_elevation_range = true,
                 const Array *p_noise = nullptr,
                 Array       *p_deposition_map = nullptr);

/**
 * @brief Carves watershed ridges using basin-wise distance transforms.
 *
 * **Example**
 * @include ex_watershed_ridge.cpp
 *
 * **Result**
 * @image html ex_watershed_ridge.png
 */
Array watershed_ridge(
    const Array        &z,
    float               amplitude = 0.2f,
    float               width = 32.f,
    float               edt_exponent = 0.5f,
    int                 prefilter_ir = 0,
    FlowDirectionMethod fd_method = FlowDirectionMethod::FDM_D8,
    const Array        *p_noise_x = nullptr,
    const Array        *p_noise_y = nullptr,
    const Array        *p_scaling = nullptr);

Array watershed_ridge(
    const Array        &z,
    const Array        *p_mask,
    float               amplitude = 0.2f,
    float               width = 32.f,
    float               edt_exponent = 0.5f,
    int                 prefilter_ir = 0,
    FlowDirectionMethod fd_method = FlowDirectionMethod::FDM_D8,
    const Array        *p_noise_x = nullptr,
    const Array        *p_noise_y = nullptr,
    const Array        *p_scaling = nullptr);

} // namespace hmap::gpu
