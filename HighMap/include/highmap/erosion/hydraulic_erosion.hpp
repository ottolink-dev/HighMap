/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once

#include <cmath>
#include <cstdint>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/erosion/erosion_parameters.hpp"
#include "highmap/interpolate/interpolate2d.hpp"
#include "highmap/math/profiles.hpp"
#include "highmap/terrain_tri_mesh.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

// neighbor pattern search based on Moore pattern and define diagonal
// weight coefficients ('c' corresponds to a weight coefficient
// applied to take into account the longer distance for diagonal
// comparison between cells)

// clang-format off
// 6 2 8
// 1 . 4
// 5 3 7
#define HMAP_DI {-1, 0, 0, 1, -1, -1, 1, 1}
#define HMAP_DJ {0, 1, -1, 0, -1, 1, -1, 1}
#define HMAP_CD  {1.f, 1.f, 1.f, 1.f, (float)M_SQRT2, (float)M_SQRT2, (float)M_SQRT2, (float)M_SQRT2}
#define HMAP_CD_INV  {1.f, 1.f, 1.f, 1.f, 1.f / (float)M_SQRT2, 1.f / (float)M_SQRT2, 1.f / (float)M_SQRT2, 1.f / (float)M_SQRT2}
// clang-format on

namespace hmap
{

/**
 * @brief Apply an algerbic formula based on the local gradient to perform
 * erosion/deposition.
 *
 * @param z                  Input array.
 * @param p_mask             Intensity mask, expected in [0, 1] (applied as a
 *                           post-processing).
 * @param talus_ref          Reference talus.
 * @param ir                 Smoothing prefilter radius.
 * @param p_bedrock          Reference to the bedrock heightmap.
 * @param p_erosion_map[out] Reference to the erosion map, provided as an output
 *                           field.
 * @param p_deposition_map   [out] Reference to the deposition map, provided as
 *                           an output field.
 * @param c_erosion          Erosion coefficient.
 * @param c_deposition       Deposition coefficient.
 * @param iterations         Number of iterations.
 *
 * **Example**
 * @include ex_hydraulic_algebric.cpp
 *
 * **Result**
 * @image html ex_hydraulic_algebric.png
 */
void hydraulic_algebric(Array &z,
                        Array *p_mask,
                        float  talus_ref,
                        int    ir,
                        Array *p_bedrock = nullptr,
                        Array *p_erosion_map = nullptr,
                        Array *p_deposition_map = nullptr,
                        float  c_erosion = 0.05f,
                        float  c_deposition = 0.05f,
                        int    iterations = 1);

/// @overload
void hydraulic_algebric(Array &z,
                        float  talus_ref,
                        int    ir,
                        Array *p_bedrock = nullptr,
                        Array *p_erosion_map = nullptr,
                        Array *p_deposition_map = nullptr,
                        float  c_erosion = 0.05f,
                        float  c_deposition = 0.05f,
                        int    iterations = 1);

/**
 * @brief Apply cell-based hydraulic erosion/deposition based on Benes et al.
 * procedure.
 *
 * See @cite Benes2002 and @cite Olsen2004.
 *
 * @param z                  Input array.
 * @param p_mask             Intensity mask, expected in [0, 1] (applied as a
 *                           post-processing).
 * @param iterations         Number of iterations.
 * @param p_bedrock          Reference to the bedrock heightmap.
 * @param p_moisture_map     Reference to the moisture map (quantity of rain),
 *                           expected to be in [0, 1].
 * @param p_erosion_map[out] Reference to the erosion map, provided as an output
 *                           field.
 * @param p_deposition_map   [out] Reference to the deposition map, provided as
 *                           an output field.
 * @param c_capacity         Sediment capacity.
 * @param c_deposition       Deposition coefficient.
 * @param c_erosion          Erosion coefficient.
 * @param water_level        Water level.
 * @param evap_rate          Water evaporation rate.
 * @param rain_rate          Rain relaxation rate.
 *
 * **Example**
 * @include ex_hydraulic_benes.cpp
 *
 * **Result**
 * @image html ex_hydraulic_benes.png
 */
void hydraulic_benes(Array &z,
                     Array *p_mask,
                     int    iterations = 50,
                     Array *p_bedrock = nullptr,
                     Array *p_moisture_map = nullptr,
                     Array *p_erosion_map = nullptr,
                     Array *p_deposition_map = nullptr,
                     float  c_capacity = 40.f,
                     float  c_erosion = 0.2f,
                     float  c_deposition = 0.8f,
                     float  water_level = 0.005f,
                     float  evap_rate = 0.01f,
                     float  rain_rate = 0.5f);

/// @overload
void hydraulic_benes(Array &z,
                     int    iterations = 50,
                     Array *p_bedrock = nullptr,
                     Array *p_moisture_map = nullptr,
                     Array *p_erosion_map = nullptr,
                     Array *p_deposition_map = nullptr,
                     float  c_capacity = 40.f,
                     float  c_erosion = 0.2f,
                     float  c_deposition = 0.8f,
                     float  water_level = 0.005f,
                     float  evap_rate = 0.01f,
                     float  rain_rate = 0.5f);

/**
 * @brief Apply cell-based hydraulic erosion using a nonlinear diffusion model.
 * @param z           Input array.
 * @param radius      Gaussian filter radius (with respect to a unit domain).
 * @param vmax        Maximum elevation for the details.
 * @param k_smoothing Smoothing factor, if any.
 *
 * **Example**
 * @include ex_hydraulic_blur.cpp
 *
 * **Result**
 * @image html ex_hydraulic_blur.png
 */
void hydraulic_blur(Array &z,
                    float  radius,
                    float  vmax,
                    float  k_smoothing = 0.1f);

/**
 * @brief Apply cell-based hydraulic erosion using a nonlinear diffusion model.
 *
 * See @cite Roering2001.
 *
 * @param z           Input array.
 * @param c_diffusion Diffusion coefficient.
 * @param talus       Reference talus (must be higher than the maximum talus of
 *                    the map).
 * @param iterations  Number of iterations.
 *
 * **Example**
 * @include ex_hydraulic_diffusion.cpp
 *
 * **Result**
 * @image html ex_hydraulic_diffusion.png
 */
void hydraulic_diffusion(Array &z,
                         float  c_diffusion,
                         float  talus,
                         int    iterations);

/**
 * @brief Simulates multiscale implicit stream-power hydraulic erosion (MISE).
 *
 * Adapted from MISE by Leonhard (https://github.com/Leonhardmaster2).
 *
 * @param z                Heightmap array to be eroded in-place.
 * @param params           Solver parameters.
 * @param p_bedrock        Optional bedrock array bounding maximum erosion
 *                         depth.
 * @param p_erodibility    Optional spatial erodibility multiplier array.
 * @param p_moisture_map   Optional precipitation / moisture weight map.
 * @param p_outlet         Optional outlet mask (>0 flags fixed base-level
 *                         cells).
 * @param p_sediment       Optional output array receiving loose sediment
 *                         thickness.
 * @param p_flow           Optional output array receiving normalized drainage
 *                         flow area.
 * @param p_erosion_map    Optional output array receiving erosion delta
 *                         (z_before - z, >= 0).
 * @param p_deposition_map Optional output array receiving deposition delta (z -
 *                         z_before, >= 0).
 *
 * **Example**
 * @include ex_hydraulic_mise.cpp
 *
 * **Result**
 * @image html ex_hydraulic_mise.png
 */
void hydraulic_mise(Array            &z,
                    const MiseParams &params = MiseParams(),
                    const Array      *p_bedrock = nullptr,
                    const Array      *p_erodibility = nullptr,
                    const Array      *p_moisture_map = nullptr,
                    const Array      *p_outlet = nullptr,
                    Array            *p_sediment = nullptr,
                    Array            *p_flow = nullptr,
                    Array            *p_erosion_map = nullptr,
                    Array            *p_deposition_map = nullptr);

/**
 * @brief Simulates multiscale implicit stream-power hydraulic erosion (MISE)
 * with an intensity mask.
 *
 * @param z                Heightmap array to be eroded in-place.
 * @param p_mask           Intensity mask in [0, 1] used for blending
 *                         post-erosion.
 * @param params           Solver parameters.
 * @param p_bedrock        Optional bedrock array bounding maximum erosion
 *                         depth.
 * @param p_erodibility    Optional spatial erodibility multiplier array.
 * @param p_moisture_map   Optional precipitation / moisture weight map.
 * @param p_outlet         Optional outlet mask (>0 flags fixed base-level
 *                         cells).
 * @param p_sediment       Optional output array receiving loose sediment
 *                         thickness.
 * @param p_flow           Optional output array receiving normalized drainage
 *                         flow area.
 * @param p_erosion_map    Optional output array receiving erosion delta
 *                         (z_before - z, >= 0).
 * @param p_deposition_map Optional output array receiving deposition delta (z -
 *                         z_before, >= 0).
 */
void hydraulic_mise(Array            &z,
                    const Array      *p_mask,
                    const MiseParams &params = MiseParams(),
                    const Array      *p_bedrock = nullptr,
                    const Array      *p_erodibility = nullptr,
                    const Array      *p_moisture_map = nullptr,
                    const Array      *p_outlet = nullptr,
                    Array            *p_sediment = nullptr,
                    Array            *p_flow = nullptr,
                    Array            *p_erosion_map = nullptr,
                    Array            *p_deposition_map = nullptr);

/**
 * @brief Apply cell-based hydraulic erosion/deposition of Musgrave et al.
 * (1989).
 *
 * A simple grid-based erosion technique was published by Musgrave, Kolb, and
 * Mace in 1989 @cite Musgrave1989.
 *
 * @param z            Input array.
 * @param moisture_map Moisture map (quantity of rain), expected to be in [0,
 *                     1].
 * @param iterations   Number of iterations.
 * @param c_capacity   Sediment capacity.
 * @param c_deposition Deposition coefficient.
 * @param c_erosion    Erosion coefficient.
 * @param water_level  Water level.
 * @param evap_rate    Water evaporation rate.
 *
 * **Example**
 * @include ex_hydraulic_musgrave.cpp
 *
 * **Result**
 * @image html ex_hydraulic_musgrave.png
 */
void hydraulic_musgrave(Array &z,
                        Array &moisture_map,
                        int    iterations = 100,
                        float  c_capacity = 1.f,
                        float  c_erosion = 0.1f,
                        float  c_deposition = 0.1f,
                        float  water_level = 0.01f,
                        float  evap_rate = 0.01f);

/// @overload
void hydraulic_musgrave(Array &z,
                        int    iterations = 100,
                        float  c_capacity = 1.f,
                        float  c_erosion = 0.1f,
                        float  c_deposition = 0.1f,
                        float  water_level = 0.01f,
                        float  evap_rate = 0.01f);

/**
 * @brief Perform hydraulic erosion on a triangulated terrain mesh.
 *
 * Iteratively updates elevations using a drainage model with uplift, slope
 * constraints, and spatially varying erodibility.
 *
 * @param mesh           Input/output terrain mesh.
 * @param erodibility    Per-vertex erodibility coefficients.
 * @param max_slope      Per-vertex maximum slope constraint.
 * @param m_exp          Exponent controlling flow response.
 * @param uplift_rate    Constant uplift applied each iteration.
 * @param tolerance      Convergence threshold on elevation updates.
 * @param max_iterations Maximum number of erosion iterations.
 *
 * **Example**
 * @include ex_hydraulic_saleve.cpp
 *
 * **Result**
 * @image html ex_hydraulic_saleve.png
 */
void hydraulic_saleve(TerrainTriMesh           &mesh,
                      const std::vector<float> &erodibility,
                      const std::vector<float> &max_slope,
                      float                     m_exp = 0.8f,
                      float                     uplift_rate = 1.f,
                      float                     tolerance = 1e-3f,
                      int                       max_iterations = 200,
                      float                     noise_strength = 0.f,
                      std::uint32_t             seed = 0,
                      bool  enable_post_slope_limiter = false,
                      float post_slope_limit = 0.f,
                      bool  enable_post_smoothing = false);

/**
 * @brief Apply hydraulic erosion to a heightmap using an adaptive mesh.
 *
 * Converts the input heightmap to a triangulated mesh, performs erosion, and
 * interpolates the result back to a grid.
 *
 * @param  z                        Input heightmap.
 * @param  seed                     Random seed for point sampling.
 * @param  control_points_count     Number of mesh control points.
 * @param  m_exp                    Flow response exponent.
 * @param  uplift_rate              Constant uplift per iteration.
 * @param  tolerance                Convergence threshold.
 * @param  max_iterations           Maximum number of iterations.
 * @param  smin                     Minimum slope constraint.
 * @param  smax                     Maximum slope constraint.
 * @param  strength                 Blending factor between original and eroded
 *                                  terrain.
 * @param  scale_erodibility_with_z Modulate erodibility with elevation.
 * @param  erodibility_distrib_exp  Exponent for erodibility distribution.
 * @param  p_noise_x                Optional X displacement field.
 * @param  p_noise_y                Optional Y displacement field.
 *
 * @return                          Eroded heightmap.
 *
 * **Example**
 * @include ex_hydraulic_saleve.cpp
 *
 * **Result**
 * @image html ex_hydraulic_saleve.png
 */
Array hydraulic_saleve(const Array          &z,
                       std::uint32_t         seed,
                       size_t                control_points_count = 10000,
                       float                 m_exp = 0.8f,
                       float                 uplift_rate = 1.f,
                       float                 tolerance = 1e-3f,
                       int                   max_iterations = 200,
                       float                 smin = 0.f,
                       float                 smax = 6.f,
                       float                 strength = 0.5f,
                       bool                  scale_erodibility_with_z = true,
                       float                 erodibility_distrib_exp = 1.f,
                       float                 noise_strength = 0.f,
                       bool                  enable_post_slope_limiter = false,
                       float                 post_slope_limit = 0.f,
                       bool                  enable_post_smoothing = false,
                       InterpolationMethod2D interpolation_method =
                           InterpolationMethod2D::ITP2D_DELAUNAY_GRADIENT,
                       const Array *p_noise_x = nullptr,
                       const Array *p_noise_y = nullptr);

Array hydraulic_saleve(const Array          &z,
                       const Array          *p_mask,
                       std::uint32_t         seed,
                       size_t                control_points_count = 10000,
                       float                 m_exp = 0.8f,
                       float                 uplift_rate = 1.f,
                       float                 tolerance = 1e-3f,
                       int                   max_iterations = 200,
                       float                 smin = 0.f,
                       float                 smax = 6.f,
                       float                 strength = 0.5f,
                       bool                  scale_erodibility_with_z = true,
                       float                 erodibility_distrib_exp = 1.f,
                       float                 noise_strength = 0.f,
                       bool                  enable_post_slope_limiter = false,
                       float                 post_slope_limit = 0.f,
                       bool                  enable_post_smoothing = false,
                       InterpolationMethod2D interpolation_method =
                           InterpolationMethod2D::ITP2D_DELAUNAY_GRADIENT,
                       const Array *p_noise_x = nullptr,
                       const Array *p_noise_y = nullptr);

/**
 * @brief Apply hydraulic erosion based on a flow accumulation map.
 *
 * @param z                  Input array.
 * @param p_mask             Intensity mask, expected in [0, 1] (applied as a
 *                           post-processing).
 * @param c_erosion          Erosion coefficient.
 * @param talus_ref          Reference talus used to localy define the
 *                           flow-partition exponent (small values of
 *                           `talus_ref` will lead to thinner flow streams, see
 *                           {@link flow_accumulation_dinf}).
 * @param p_bedrock          Lower elevation limit.
 * @param p_moisture_map     Reference to the moisture map (quantity of rain),
 *                           expected to be in [0, 1].
 * @param p_erosion_map[out] Reference to the erosion map, provided as an output
 *                           field.
 * @param ir                 Kernel radius. If `ir > 1`, a cone kernel is used
 *                           to carv channel flow erosion.
 * @param clipping_ratio     Flow accumulation clipping ratio.
 *
 * **Example**
 * @include ex_hydraulic_stream.cpp
 *
 * **Result**
 * @image html ex_hydraulic_stream0.png
 * @image html ex_hydraulic_stream1.png
 */
void hydraulic_stream(Array       &z,
                      float        c_erosion,
                      float        talus_ref,
                      const Array *p_bedrock = nullptr,
                      const Array *p_moisture_map = nullptr,
                      Array       *p_erosion_map = nullptr, // -> out
                      int          ir = 1,
                      float        clipping_ratio = 10.f);

/// @overload
void hydraulic_stream(Array       &z,
                      const Array *p_mask,
                      float        c_erosion,
                      float        talus_ref,
                      const Array *p_bedrock = nullptr,
                      const Array *p_moisture_map = nullptr,
                      Array       *p_erosion_map = nullptr,
                      int          ir = 1,
                      float        clipping_ratio = 10.f);

/**
 * @brief Apply hydraulic erosion based on a flow accumulation map, alternative
 * formulation.
 *
 * @param z                      Input array representing the terrain elevation.
 * @param c_erosion              Erosion coefficient controlling the intensity
 *                               of erosion.
 * @param talus_ref              Reference talus used to locally define the
 *                               flow-partition exponent. Small values lead to
 *                               thinner flow streams (see
 *                               {@link flow_accumulation_dinf}).
 * @param deposition_ir          Kernel radius for sediment deposition. If
 *                               greater than 1, a smoothing effect is applied.
 * @param deposition_scale_ratio Scaling factor for sediment deposition.
 * @param gradient_power         Exponent applied to the terrain gradient to
 *                               control erosion intensity.
 * @param gradient_scaling_ratio Scaling factor for gradient-based erosion.
 * @param gradient_prefilter_ir  Kernel radius for pre-filtering the terrain
 *                               gradient.
 * @param saturation_ratio       Ratio controlling the water saturation
 *                               threshold for erosion processes.
 * @param p_bedrock              Pointer to an optional lower elevation limit.
 * @param p_moisture_map         Pointer to the moisture map (rainfall
 *                               quantity), expected to be in [0, 1].
 * @param p_erosion_map[out]     Pointer to the erosion map, provided as an
 *                               output field.
 * @param p_flow_map[out]        Pointer to the flow accumulation map, provided
 *                               as an output field.
 * @param ir                     Kernel radius. If `ir > 1`, a cone kernel is
 *                               used to carve channel flow erosion.
 *
 * **Example**
 * @include ex_hydraulic_stream.cpp
 *
 * **Result**
 * @image html ex_hydraulic_stream0.png
 * @image html ex_hydraulic_stream1.png
 */
void hydraulic_stream_log(Array       &z,
                          float        c_erosion,
                          float        talus_ref,
                          int          deposition_ir = 32,
                          float        deposition_scale_ratio = 1.f,
                          float        gradient_power = 0.8f,
                          float        gradient_scaling_ratio = 1.f,
                          int          gradient_prefilter_ir = 16,
                          float        saturation_ratio = 1.f,
                          const Array *p_bedrock = nullptr,
                          const Array *p_moisture_map = nullptr,
                          Array       *p_erosion_map = nullptr,
                          Array       *p_deposition_map = nullptr,
                          Array       *p_flow_map = nullptr);

/// @overload
void hydraulic_stream_log(Array       &z,
                          float        c_erosion,
                          float        talus_ref,
                          const Array *p_mask,
                          int          deposition_ir = 32,
                          float        deposition_scale_ratio = 1.f,
                          float        gradient_power = 0.8f,
                          float        gradient_scaling_ratio = 1.f,
                          int          gradient_prefilter_ir = 16,
                          float        saturation_ratio = 1.f,
                          const Array *p_bedrock = nullptr,
                          const Array *p_moisture_map = nullptr,
                          Array       *p_erosion_map = nullptr,
                          Array       *p_deposition_map = nullptr,
                          Array       *p_flow_map = nullptr);

/**
 * @brief Multiscale stream-power hydraulic erosion cascade.
 *
 * Runs stream-power erosion over a geometric resolution ladder (coarsest first,
 * up to z.shape). Coarse levels carve broad valleys and major drainage
 * networks, while fine levels carve detailed tributaries.
 *
 * @param z                      Heightmap array to modify.
 * @param c_erosion              Erosion coefficient.
 * @param talus_ref              Reference talus.
 * @param levels                 Number of resolution levels.
 * @param deposition_ir          Deposition smoothing radius at full resolution.
 * @param deposition_scale_ratio Deposition blend factor (applied at last
 *                               scale).
 * @param gradient_power         Power applied to normalized gradient.
 * @param gradient_scaling_ratio Gradient scaling weight.
 * @param gradient_prefilter_ir  Gradient prefilter smoothing radius at full
 *                               resolution.
 * @param saturation_ratio       Flow saturation ratio.
 * @param p_bedrock              Optional bedrock array bounding maximum erosion
 *                               depth.
 * @param p_moisture_map         Optional moisture map.
 * @param p_erosion_map          Optional output erosion map.
 * @param p_deposition_map       Optional output deposition map.
 * @param p_flow_map             Optional output flow accumulation map.
 * @param mix                    Interpolation blend factor between levels.
 * @param warp                   Domain warp amplitude during upsampling.
 * @param seed                   Random seed for domain warping noise.
 */
void hydraulic_stream_log_multiscale(Array        &z,
                                     float         c_erosion,
                                     float         talus_ref,
                                     int           levels = 3,
                                     int           deposition_ir = 32,
                                     float         deposition_scale_ratio = 1.f,
                                     float         gradient_power = 0.8f,
                                     float         gradient_scaling_ratio = 1.f,
                                     int           gradient_prefilter_ir = 16,
                                     float         saturation_ratio = 1.f,
                                     const Array  *p_bedrock = nullptr,
                                     const Array  *p_moisture_map = nullptr,
                                     Array        *p_erosion_map = nullptr,
                                     Array        *p_deposition_map = nullptr,
                                     Array        *p_flow_map = nullptr,
                                     float         mix = 1.f,
                                     float         warp = 0.35f,
                                     std::uint32_t seed = 1);

/// @overload
void hydraulic_stream_log_multiscale(Array        &z,
                                     float         c_erosion,
                                     float         talus_ref,
                                     const Array  *p_mask,
                                     int           levels = 3,
                                     int           deposition_ir = 32,
                                     float         deposition_scale_ratio = 1.f,
                                     float         gradient_power = 0.8f,
                                     float         gradient_scaling_ratio = 1.f,
                                     int           gradient_prefilter_ir = 16,
                                     float         saturation_ratio = 1.f,
                                     const Array  *p_bedrock = nullptr,
                                     const Array  *p_moisture_map = nullptr,
                                     Array        *p_erosion_map = nullptr,
                                     Array        *p_deposition_map = nullptr,
                                     Array        *p_flow_map = nullptr,
                                     float         mix = 1.f,
                                     float         warp = 0.35f,
                                     std::uint32_t seed = 1);

/**
 * @brief Applies hydraulic erosion with upscaling amplification.
 *
 * This function progressively upscales the input array `z` by powers of 2 and
 * applies hydraulic erosion based on flow accumulation at each level of
 * upscaling. After all upscaling levels are processed, the array is resampled
 * back to its original resolution using bilinear interpolation.
 *
 * @param z                Input array representing elevation data.
 * @param c_erosion        Erosion coefficient.
 * @param talus_ref        Reference talus used to locally define the
 *                         flow-partition exponent. Smaller values lead to
 *                         thinner flow streams.
 * @param upscaling_levels Number of upscaling levels to apply. The function
 *                         will resample the array at each level.
 * @param persistence      A scaling factor applied at each level to adjust the
 *                         impact of the unary operation. Higher persistence
 *                         values will amplify the effects at each level.
 * @param ir               Kernel radius. If `ir > 1`, a cone kernel is used to
 *                         carve channel flow erosion.
 * @param clipping_ratio   Flow accumulation clipping ratio.
 *
 * @note The function first applies upscaling using bicubic resampling, performs
 * hydraulic erosion at each level, and finally resamples the array back to its
 * initial resolution using bilinear interpolation.
 *
 * **Example**
 * @include ex_hydraulic_stream_upscale_amplification.cpp
 *
 * **Result**
 * @image html ex_hydraulic_stream_upscale_amplification.png
 */
void hydraulic_stream_upscale_amplification(Array &z,
                                            float  c_erosion,
                                            float  talus_ref,
                                            int    upscaling_levels = 1,
                                            float  persistence = 1.f,
                                            int    ir = 1,
                                            float  clipping_ratio = 10.f);

/**
 * @brief Applies hydraulic erosion with upscaling amplification, with a
 * post-processing intensity mask.
 *
 * Similar to the overloaded version, this function progressively upscales the
 * input array `z` and applies hydraulic erosion. Additionally, an intensity
 * mask `p_mask` is applied as a post-processing step.
 *
 * @param z                Input array representing elevation data.
 * @param p_mask           Intensity mask, expected in [0, 1], which is applied
 *                         as a post-processing step.
 * @param c_erosion        Erosion coefficient.
 * @param talus_ref        Reference talus used to locally define the
 *                         flow-partition exponent. Smaller values lead to
 *                         thinner flow streams.
 * @param upscaling_levels Number of upscaling levels to apply. The function
 *                         will resample the array at each level.
 * @param persistence      A scaling factor applied at each level to adjust the
 *                         impact of the unary operation. Higher persistence
 *                         values will amplify the effects at each level.
 * @param ir               Kernel radius. If `ir > 1`, a cone kernel is used to
 *                         carve channel flow erosion.
 * @param clipping_ratio   Flow accumulation clipping ratio.
 *
 * @note This version of the function applies an additional intensity mask as
 * part of the upscaling amplification process.
 *
 * **Example**
 * @include ex_hydraulic_stream_upscale_amplification.cpp
 *
 * **Result**
 * @image html ex_hydraulic_stream_upscale_amplification.png
 */

/// @overload
void hydraulic_stream_upscale_amplification(Array       &z,
                                            const Array *p_mask,
                                            float        c_erosion,
                                            float        talus_ref,
                                            int          upscaling_levels = 1,
                                            float        persistence = 1.f,
                                            int          ir = 1,
                                            float        clipping_ratio = 10.f);
} // namespace hmap

namespace hmap::gpu
{

/**
 * @brief Performs iterative particle-based convolution erosion on a heightmap.
 *
 * This function simulates erosion by spawning particles, tracing their paths
 * over the heightmap, accumulating a mask and size field, and applying a
 * kernel-based convolution to compute erosion deltas.
 *
 * The process is repeated for a number of iterations to refine the result.
 *
 * @param z                Heightmap to be eroded (modified in-place).
 * @param seed             RNG seed for deterministic behavior.
 * @param iterations       Number of erosion iterations.
 * @param particle_count   Number of particles per iteration.
 * @param ir_min           Minimum kernel radius scale.
 * @param ir_max           Maximum kernel radius scale.
 * @param size_distrib_exp Exponent controlling particle size bias.
 * @param erosion_strength Global strength of erosion applied per iteration.
 * @param randomness       Controls randomness in particle trajectories.
 * @param exit_forcing     Forcing strength towards the domain frontier.
 *
 *  **Example**
 * @include ex_conv_erosion.cpp
 *
 * **Result**
 * @image html ex_conv_erosion.png
 */
void conv_erosion(Array        &z,
                  std::uint32_t seed,
                  int           iterations = 20,
                  int           particle_count = 1000,
                  int           ir_min = 8,
                  int           ir_max = 64,
                  float         size_distrib_exp = 1.f,
                  float         erosion_strength = 0.02f,
                  float         randomness = 0.01f,
                  float         exit_forcing = 0.05f,
                  int           gradient_ir = 16,
                  float         gradient_exp = 0.5f,
                  float         gradient_strength_min = 0.f);

/**
 * @brief Fill holes using Gaussian-based deposition.
 *
 * Applies a smoothing/deposition pass that fills local depressions while
 * preserving overall terrain shape through gradient-aware blending.
 *
 * @param z                   Heightmap to modify (in-place).
 * @param deposition_ir       Influence radius of the Gaussian filter.
 * @param deposition_strength Blending factor controlling deposition intensity.
 * @param iterations          Number of successive deposition passes.
 */
void deposition_fill_holes(Array &z,
                           int    deposition_ir,
                           float  deposition_strength,
                           int    iterations = 1);

/// @overload
void deposition_fill_holes(Array       &z,
                           int          deposition_ir,
                           float        deposition_strength,
                           const Array *p_mask,
                           int          iterations = 1);

/**
 * @brief Simulates hydraulic erosion on a heightmap using particle-based flow.
 *
 * Particles traverse the heightmap `z`, eroding and depositing material
 * according to local slope, capacity, inertia, and optional directional bias.
 *
 * @param z                 Heightmap array to modify.
 * @param nparticles        Number of erosion particles to simulate.
 * @param seed              Random seed for particle initialization.
 * @param p_bedrock         Optional bedrock array to limit erosion.
 * @param p_moisture_map    Optional moisture map affecting erosion/deposition.
 * @param p_elevation_shift Optional elevation shift map.
 * @param p_erosion_map     Optional output array recording total erosion.
 * @param p_deposition_map  Optional output array recording deposition.
 * @param c_capacity        Sediment capacity of each particle.
 * @param c_erosion         Erosion rate coefficient.
 * @param c_deposition      Deposition rate coefficient.
 * @param c_inertia         Particle inertia factor.
 * @param c_gravity         Gravity effect on particle movement.
 * @param drag_rate         Particle velocity damping per step.
 * @param evap_rate         Sediment evaporation rate.
 * @param talus_slope       Talus slope angle.
 * @param collapse_rate     Collapse rate.
 * @param iterations        Number of multi-pass iterations.
 *
 *  **Example**
 * @include ex_hydraulic_particle.cpp
 *
 * **Result**
 * @image html ex_hydraulic_particle.png
 */
void hydraulic_particle(Array        &z,
                        int           nparticles,
                        std::uint32_t seed,
                        const Array  *p_bedrock = nullptr,
                        const Array  *p_moisture_map = nullptr,
                        const Array  *p_elevation_shift = nullptr,
                        Array        *p_erosion_map = nullptr,
                        Array        *p_deposition_map = nullptr,
                        float         c_capacity = 10.f,
                        float         c_erosion = 0.05f,
                        float         c_deposition = 0.05f,
                        float         c_inertia = 0.01f,
                        float         c_gravity = 1.f,
                        float         drag_rate = 0.001f,
                        float         evap_rate = 0.001f,
                        float         talus_slope = 2.f,
                        float         collapse_rate = 0.1f,
                        int           iterations = 1);

void hydraulic_particle(Array        &z,
                        const Array  *p_mask,
                        int           nparticles,
                        std::uint32_t seed,
                        const Array  *p_bedrock = nullptr,
                        const Array  *p_moisture_map = nullptr,
                        const Array  *p_elevation_shift = nullptr,
                        Array        *p_erosion_map = nullptr,
                        Array        *p_deposition_map = nullptr,
                        float         c_capacity = 10.f,
                        float         c_erosion = 0.05f,
                        float         c_deposition = 0.05f,
                        float         c_inertia = 0.01f,
                        float         c_gravity = 1.f,
                        float         drag_rate = 0.001f,
                        float         evap_rate = 0.001f,
                        float         talus_slope = 2.f,
                        float         collapse_rate = 0.1f,
                        int           iterations = 1);

/**
 * @brief Multiscale particle-based hydraulic erosion cascade.
 *
 * Runs particle erosion over a geometric resolution ladder (coarsest first, up
 * to z.shape). Coarse levels carve broad continental valleys and major river
 * corridors, while fine levels carve detailed tributaries and gullies.
 *
 * @param z                 Heightmap array to modify.
 * @param seed              Random seed.
 * @param steps_per_level   Number of iteration steps per level.
 * @param p_bedrock         Optional bedrock array.
 * @param p_moisture_map    Optional moisture map.
 * @param p_elevation_shift Optional elevation shift map.
 * @param p_erosion_map     Optional output erosion map.
 * @param p_deposition_map  Optional output deposition map.
 * @param particles_ratio   Particles per level as a multiple of resolution.
 * @param c_capacity        Sediment capacity factor.
 * @param c_erosion         Erosion rate coefficient.
 * @param c_deposition      Deposition rate coefficient.
 * @param c_inertia         Inertia factor.
 * @param c_gravity         Gravity factor.
 * @param drag_rate         Velocity damping.
 * @param evap_rate         Evaporation rate.
 * @param talus_slope       Domain-normalized talus slope threshold.
 * @param collapse_rate     Bank collapse rate.
 * @param mix               Interpolation blend factor.
 */
void hydraulic_particle_multiscale(
    Array                  &z,
    std::uint32_t           seed,
    const std::vector<int> &steps_per_level = {4, 2, 1},
    const Array            *p_bedrock = nullptr,
    const Array            *p_moisture_map = nullptr,
    const Array            *p_elevation_shift = nullptr,
    Array                  *p_erosion_map = nullptr,
    Array                  *p_deposition_map = nullptr,
    float                   particles_ratio = 0.5f,
    float                   c_capacity = 10.f,
    float                   c_erosion = 0.05f,
    float                   c_deposition = 0.05f,
    float                   c_inertia = 0.01f,
    float                   c_gravity = 1.f,
    float                   drag_rate = 0.001f,
    float                   evap_rate = 0.001f,
    float                   talus_slope = 2.f,
    float                   collapse_rate = 0.1f,
    float                   mix = 1.f);

void hydraulic_particle_multiscale(
    Array                  &z,
    const Array            *p_mask,
    std::uint32_t           seed,
    const std::vector<int> &steps_per_level = {4, 2, 1},
    const Array            *p_bedrock = nullptr,
    const Array            *p_moisture_map = nullptr,
    const Array            *p_elevation_shift = nullptr,
    Array                  *p_erosion_map = nullptr,
    Array                  *p_deposition_map = nullptr,
    float                   particles_ratio = 0.5f,
    float                   c_capacity = 10.f,
    float                   c_erosion = 0.05f,
    float                   c_deposition = 0.05f,
    float                   c_inertia = 0.01f,
    float                   c_gravity = 1.f,
    float                   drag_rate = 0.001f,
    float                   evap_rate = 0.001f,
    float                   talus_slope = 2.f,
    float                   collapse_rate = 0.1f,
    float                   mix = 1.f);

/**
 * @brief Particle-based hydraulic erosion with flow-field coupling (McDonald's
 * model).
 *
 * @param z               Input/output heightmap.
 * @param steps           Number of erosion iterations.
 * @param seed            Random seed number.
 * @param params          Parameters structure.
 * @param p_moisture_map  Optional moisture map.
 * @param p_sediment_map  Optional output: final sediment layer.
 * @param p_discharge_map Optional output: water discharge field.
 *
 * **Example**
 * @include ex_hydraulic_mcdonald.cpp
 *
 * **Result**
 * @image html ex_hydraulic_mcdonald0.png
 * @image html ex_hydraulic_mcdonald1.png
 * @image html ex_hydraulic_mcdonald2.png
 */
void hydraulic_mcdonald(Array                &z,
                        int                   steps,
                        std::uint32_t         seed,
                        const McDonaldParams &params,
                        const Array          *p_moisture_map = nullptr,
                        Array                *p_sediment_map = nullptr,
                        Array                *p_discharge_map = nullptr);

void hydraulic_mcdonald(Array        &z,
                        int           steps,
                        std::uint32_t seed,
                        const Array  *p_moisture_map = nullptr,
                        Array        *p_sediment_map = nullptr,
                        Array        *p_discharge_map = nullptr,
                        float         world_extent_km = 40.f,
                        float         z_scale_km = 4.f,
                        int           samples = 8192,
                        int           maxage = 512,
                        float         lrate = 0.2f,
                        float         time_step = 10.f,
                        float         rainfall = 1.f,
                        float         evap_rate = 1e-4f,
                        float         gravity = 9.81f,
                        float         viscosity = 0.025f,
                        float         bed_shear = 0.01f,
                        float         crit_slope = 0.57f,
                        float         settle_rate = 0.1f,
                        float         thermal_rate = 2.5e-3f,
                        float         deposition_rate = 5e-3f,
                        float         suspension_rate = 2.5e-4f,
                        float         exit_slope = 0.01f);

/**
 * @brief Multiscale driver for hydraulic_mcdonald.
 *
 * **Example**
 * @include ex_hydraulic_mcdonald.cpp
 *
 * **Result**
 * @image html ex_hydraulic_mcdonald3.png
 */
void hydraulic_mcdonald_multiscale(Array                  &z,
                                   std::uint32_t           seed,
                                   const std::vector<int> &steps_per_level,
                                   const McDonaldParams   &params,
                                   const Array *p_moisture_map = nullptr,
                                   Array       *p_sediment_map = nullptr,
                                   Array       *p_discharge_map = nullptr);

void hydraulic_mcdonald_multiscale(
    Array                  &z,
    std::uint32_t           seed,
    const std::vector<int> &steps_per_level = {512, 256, 128},
    const Array            *p_moisture_map = nullptr,
    Array                  *p_sediment_map = nullptr,
    Array                  *p_discharge_map = nullptr,
    float                   world_extent_km = 40.f,
    float                   z_scale_km = 4.f,
    int                     samples = 8192,
    int                     maxage = 512,
    float                   lrate = 0.2f,
    float                   time_step = 10.f,
    float                   rainfall = 1.f,
    float                   evap_rate = 1e-4f,
    float                   gravity = 9.81f,
    float                   viscosity = 0.025f,
    float                   bed_shear = 0.01f,
    float                   crit_slope = 0.57f,
    float                   settle_rate = 0.1f,
    float                   thermal_rate = 2.5e-3f,
    float                   deposition_rate = 5e-3f,
    float                   suspension_rate = 2.5e-4f,
    float                   exit_slope = 0.01f);

/**
 * @brief Apply cell-based hydraulic erosion/deposition of Musgrave et al.
 * (1989) on GPU.
 *
 * **Example**
 * @include ex_hydraulic_musgrave.cpp
 *
 * **Result**
 * @image html ex_hydraulic_musgrave.png
 */
void hydraulic_musgrave(Array &z,
                        Array &moisture_map,
                        int    iterations = 100,
                        float  c_capacity = 1.f,
                        float  c_erosion = 0.1f,
                        float  c_deposition = 0.1f,
                        float  water_level = 0.01f,
                        float  evap_rate = 0.01f);

/// @overload
void hydraulic_musgrave(Array &z,
                        int    iterations = 100,
                        float  c_capacity = 1.f,
                        float  c_erosion = 0.1f,
                        float  c_deposition = 0.1f,
                        float  water_level = 0.01f,
                        float  evap_rate = 0.01f);

/**
 * @brief Apply phase-guided hydraulic procedural erosion to a heightmap.
 *
 * **Example**
 * @include ex_hydraulic_procedural.cpp
 *
 * **Result**
 * @image html ex_hydraulic_procedural.png
 */
void hydraulic_procedural(
    Array         &z,
    float          kp_global,
    float          c_erosion,
    std::uint32_t  seed,
    ErosionProfile erosion_profile = ErosionProfile::EP_TRIANGLE_GRENIER,
    float          erosion_profile_parameter = 0.01f,
    float          angle_shift = 0.f, // degs
    float          phase_smoothing = 0.1f,
    float          talus_ref = 0.001f,
    float          gradient_scaling_ratio = 1.f,
    float          gradient_power = 0.8f,
    bool           exclude_ridges = true,
    bool           apply_deposition = false,
    float          deposition_strength = 1.f,
    bool           enable_default_noise = true,
    float          noise_amp = 0.01f,
    const Array   *p_kp_multiplier = nullptr,
    const Array   *p_angle_shift = nullptr,
    const Array   *p_noise_x = nullptr,
    const Array   *p_noise_y = nullptr,
    Array         *p_ridge_mask = nullptr, // ouptput
    glm::vec4      bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Multi-octave (fBm) variant of hydraulic_procedural().
 *
 * **Example**
 * @include ex_hydraulic_procedural.cpp
 *
 * **Result**
 * @image html ex_hydraulic_procedural.png
 */
void hydraulic_procedural_fbm(
    Array         &z,
    float          kp_global,
    float          c_erosion,
    std::uint32_t  seed,
    ErosionProfile erosion_profile = ErosionProfile::EP_TRIANGLE_GRENIER,
    int            octaves = 3,
    float          persistence = 0.5f,
    float          lacunarity = 2.f,
    float          erosion_profile_parameter = 0.01f,
    float          angle_shift = 0.f, // degs
    float          phase_smoothing = 0.1f,
    float          talus_ref = 0.001f,
    float          gradient_scaling_ratio = 1.f,
    float          gradient_power = 0.8f,
    bool           exclude_ridges = true,
    bool           apply_deposition = false,
    float          deposition_strength = 1.f,
    bool           enable_default_noise = true,
    float          noise_amp = 0.01f,
    const Array   *p_kp_multiplier = nullptr,
    const Array   *p_angle_shift = nullptr,
    const Array   *p_noise_x = nullptr,
    const Array   *p_noise_y = nullptr,
    Array         *p_ridge_mask = nullptr, // ouptput
    glm::vec4      bbox = {0.f, 1.f, 0.f, 1.f});

void hydraulic_procedural_fbm(
    Array         &z,
    float          kp_global,
    float          c_erosion,
    std::uint32_t  seed,
    const Array   *p_mask,
    ErosionProfile erosion_profile = ErosionProfile::EP_TRIANGLE_GRENIER,
    int            octaves = 3,
    float          persistence = 0.5f,
    float          lacunarity = 2.f,
    float          erosion_profile_parameter = 0.01f,
    float          angle_shift = 0.f, // degs
    float          phase_smoothing = 0.1f,
    float          talus_ref = 0.001f,
    float          gradient_scaling_ratio = 1.f,
    float          gradient_power = 0.8f,
    bool           exclude_ridges = true,
    bool           apply_deposition = false,
    float          deposition_strength = 1.f,
    bool           enable_default_noise = true,
    float          noise_amp = 0.01f,
    const Array   *p_kp_multiplier = nullptr,
    const Array   *p_angle_shift = nullptr,
    const Array   *p_noise_x = nullptr,
    const Array   *p_noise_y = nullptr,
    Array         *p_ridge_mask = nullptr, // ouptput
    glm::vec4      bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Simulates hydraulic erosion and deposition on a heightmap using Schott
 * method.
 *
 * **Example**
 * @include ex_hydraulic_schott.cpp
 *
 * **Result**
 * @image html ex_hydraulic_schott.png
 */
void hydraulic_schott(Array       &z,
                      int          iterations,
                      const Array &talus,
                      float        c_erosion = 1.f,
                      float        c_thermal = 0.1f,
                      float        c_deposition = 0.2f,
                      float        flow_acc_exponent = 0.8f,
                      float        flow_acc_exponent_depo = 0.8f,
                      float        flow_routing_exponent = 1.3f,
                      float        thermal_weight = 1.5f,
                      float        deposition_weight = 2.5f,
                      Array       *p_flow = nullptr);

/// @overload
void hydraulic_schott(Array       &z,
                      const Array *p_mask,
                      int          iterations,
                      const Array &talus,
                      float        c_erosion = 1.f,
                      float        c_thermal = 0.1f,
                      float        c_deposition = 0.2f,
                      float        flow_acc_exponent = 0.8f,
                      float        flow_acc_exponent_depo = 0.8f,
                      float        flow_routing_exponent = 1.3f,
                      float        thermal_weight = 1.5f,
                      float        deposition_weight = 2.5f,
                      Array       *p_flow = nullptr);

/*! @brief See hmap::gpu::hydraulic_schott */
void hydraulic_schott_erosion(Array       &z,
                              int          iterations,
                              float        c_erosion = 1.f,
                              float        flow_acc_exponent = 0.8f,
                              float        flow_routing_exponent = 1.3f,
                              const Array *p_moisture_map = nullptr,
                              Array       *p_flow = nullptr);

/*! @brief See hmap::hydraulic_stream_log */
void hydraulic_stream_log(Array &z,
                          float  c_erosion,
                          float  talus_ref,
                          int    deposition_ir = 32,
                          float  deposition_scale_ratio = 1.f,
                          float  gradient_power = 0.8f,
                          float  gradient_scaling_ratio = 1.f,
                          int    gradient_prefilter_ir = 16,
                          float  saturation_ratio = 1.f,
                          Array *p_bedrock = nullptr,
                          Array *p_moisture_map = nullptr,
                          Array *p_erosion_map = nullptr,
                          Array *p_deposition_map = nullptr,
                          Array *p_flow_map = nullptr);

/// @overload
void hydraulic_stream_log(Array       &z,
                          float        c_erosion,
                          float        talus_ref,
                          const Array *p_mask,
                          int          deposition_ir = 32,
                          float        deposition_scale_ratio = 1.f,
                          float        gradient_power = 0.8f,
                          float        gradient_scaling_ratio = 1.f,
                          int          gradient_prefilter_ir = 16,
                          float        saturation_ratio = 1.f,
                          Array       *p_bedrock = nullptr,
                          Array       *p_moisture_map = nullptr,
                          Array       *p_erosion_map = nullptr,
                          Array       *p_deposition_map = nullptr,
                          Array       *p_flow_map = nullptr);

/*! @brief See hmap::hydraulic_stream_log_multiscale */
void hydraulic_stream_log_multiscale(Array        &z,
                                     float         c_erosion,
                                     float         talus_ref,
                                     int           levels = 3,
                                     int           deposition_ir = 32,
                                     float         deposition_scale_ratio = 1.f,
                                     float         gradient_power = 0.8f,
                                     float         gradient_scaling_ratio = 1.f,
                                     int           gradient_prefilter_ir = 16,
                                     float         saturation_ratio = 1.f,
                                     const Array  *p_bedrock = nullptr,
                                     const Array  *p_moisture_map = nullptr,
                                     Array        *p_erosion_map = nullptr,
                                     Array        *p_deposition_map = nullptr,
                                     Array        *p_flow_map = nullptr,
                                     float         mix = 1.f,
                                     float         warp = 0.35f,
                                     std::uint32_t seed = 1);

/// @overload
void hydraulic_stream_log_multiscale(Array        &z,
                                     float         c_erosion,
                                     float         talus_ref,
                                     const Array  *p_mask,
                                     int           levels = 3,
                                     int           deposition_ir = 32,
                                     float         deposition_scale_ratio = 1.f,
                                     float         gradient_power = 0.8f,
                                     float         gradient_scaling_ratio = 1.f,
                                     int           gradient_prefilter_ir = 16,
                                     float         saturation_ratio = 1.f,
                                     const Array  *p_bedrock = nullptr,
                                     const Array  *p_moisture_map = nullptr,
                                     Array        *p_erosion_map = nullptr,
                                     Array        *p_deposition_map = nullptr,
                                     Array        *p_flow_map = nullptr,
                                     float         mix = 1.f,
                                     float         warp = 0.35f,
                                     std::uint32_t seed = 1);

/*! @brief See hmap::gpu::hydraulic_vpipes */
void hydraulic_vpipes(Array &z,
                      float  water_height = 1e-2f,
                      bool   maintain_water_volume = true,
                      float  evap_rate = 0.1f,
                      int    iterations = 50,
                      float  dt = 0.5f,
                      float  k_capacity = 0.5f,
                      float  k_erode = 0.001f,
                      float  k_depose = 0.01f,
                      float  k_discharge_exp = 1.f,
                      float  downcutting_max_depth_ratio = 10.f,
                      bool   flux_diffusion = true,
                      float  flux_diffusion_strength = 0.01f,
                      Array *p_rain_map = nullptr,
                      Array *p_water_depth = nullptr,
                      Array *p_sediment = nullptr,
                      Array *p_vel_u = nullptr,
                      Array *p_vel_v = nullptr);

} // namespace hmap::gpu

namespace hmap::va
{

VirtualArray hydraulic_saleve(
    const ComputeMode    &cm,
    const VirtualArray   &z,
    std::uint32_t         seed,
    size_t                control_points_count = 10000,
    float                 m_exp = 0.8f,
    float                 uplift_rate = 1.f,
    float                 tolerance = 1e-3f,
    int                   max_iterations = 200,
    float                 smin = 0.f,
    float                 smax = 6.f,
    float                 strength = 0.5f,
    bool                  scale_erodibility_with_z = true,
    float                 erodibility_distrib_exp = 1.f,
    float                 noise_strength = 0.f,
    bool                  enable_post_slope_limiter = false,
    float                 post_slope_limit = 0.f,
    bool                  enable_post_smoothing = false,
    InterpolationMethod2D interpolation_method =
        InterpolationMethod2D::ITP2D_DELAUNAY_GRADIENT,
    const VirtualArray *p_noise_x = nullptr,
    const VirtualArray *p_noise_y = nullptr,
    const VirtualArray *p_mask = nullptr);

} // namespace hmap::va
