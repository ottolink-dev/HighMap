/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once

#include <cstdint>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/interpolate/interpolate2d.hpp"
#include "highmap/terrain_tri_mesh.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

// --- Deprecated Hydraulic Functions

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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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
 * @param c_erosion          Erosion coefficient.
 * @param c_deposition       Deposition coefficient.
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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
void hydraulic_diffusion(Array &z,
                         float  c_diffusion,
                         float  talus,
                         int    iterations);

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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
void hydraulic_musgrave(Array &z,
                        Array &moisture_map,
                        int    iterations = 100,
                        float  c_capacity = 1.f,
                        float  c_erosion = 0.1f,
                        float  c_deposition = 0.1f,
                        float  water_level = 0.01f,
                        float  evap_rate = 0.01f);

/// @overload
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
void hydraulic_musgrave(Array &z,
                        int    iterations = 100,
                        float  c_capacity = 1.f,
                        float  c_erosion = 0.1f,
                        float  c_deposition = 0.1f,
                        float  water_level = 0.01f,
                        float  evap_rate = 0.01f);

/**
 * @brief Apply hydraulic erosion based on a flow accumulation map.
 *
 * @param z                  Input array.
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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
void hydraulic_stream(Array       &z,
                      float        c_erosion,
                      float        talus_ref,
                      const Array *p_bedrock = nullptr,
                      const Array *p_moisture_map = nullptr,
                      Array       *p_erosion_map = nullptr, // -> out
                      int          ir = 1,
                      float        clipping_ratio = 10.f);

/// @overload
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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
                          const Array *p_erodibility = nullptr,
                          const Array *p_moisture_map = nullptr,
                          Array       *p_erosion_map = nullptr,
                          Array       *p_deposition_map = nullptr,
                          Array       *p_flow_map = nullptr);

/// @overload
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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
                          const Array *p_erodibility = nullptr,
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
 * @param p_erodibility          Optional spatial erodibility multiplier array.
 * @param p_moisture_map         Optional moisture map.
 * @param p_erosion_map          Optional output erosion map.
 * @param p_deposition_map       Optional output deposition map.
 * @param p_flow_map             Optional output flow accumulation map.
 * @param mix                    Interpolation blend factor between levels.
 * @param warp                   Domain warp amplitude during upsampling.
 * @param seed                   Random seed for domain warping noise.
 */
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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
                                     const Array  *p_erodibility = nullptr,
                                     const Array  *p_moisture_map = nullptr,
                                     Array        *p_erosion_map = nullptr,
                                     Array        *p_deposition_map = nullptr,
                                     Array        *p_flow_map = nullptr,
                                     float         mix = 1.f,
                                     float         warp = 0.35f,
                                     std::uint32_t seed = 1);

/// @overload
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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
                                     const Array  *p_erodibility = nullptr,
                                     const Array  *p_moisture_map = nullptr,
                                     Array        *p_erosion_map = nullptr,
                                     Array        *p_deposition_map = nullptr,
                                     Array        *p_flow_map = nullptr,
                                     float         mix = 1.f,
                                     float         warp = 0.35f,
                                     std::uint32_t seed = 1);

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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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

/// @overload
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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

} // namespace hmap

namespace hmap::gpu
{

// --- Deprecated GPU Hydraulic Functions

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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
void hydraulic_musgrave(Array &z,
                        Array &moisture_map,
                        int    iterations = 100,
                        float  c_capacity = 1.f,
                        float  c_erosion = 0.1f,
                        float  c_deposition = 0.1f,
                        float  water_level = 0.01f,
                        float  evap_rate = 0.01f);

/// @overload
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
void hydraulic_musgrave(Array &z,
                        int    iterations = 100,
                        float  c_capacity = 1.f,
                        float  c_erosion = 0.1f,
                        float  c_deposition = 0.1f,
                        float  water_level = 0.01f,
                        float  evap_rate = 0.01f);

} // namespace hmap::gpu

namespace hmap::va
{

// --- Deprecated VirtualArray Hydraulic Functions

// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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
