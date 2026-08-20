/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once
#include <cmath>
#include <cstdint>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/interpolate/interpolate2d.hpp"
#include "highmap/math/profiles.hpp"
#include "highmap/terrain_tri_mesh.hpp"

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
#define HMAP_CD  {1.f, 1.f, 1.f, 1.f, M_SQRT2, M_SQRT2, M_SQRT2, M_SQRT2}
#define HMAP_CD_INV  {1.f, 1.f, 1.f, 1.f, M_SQRT2, M_SQRT2, M_SQRT2, M_SQRT2}
// clang-format on

namespace hmap
{

/**
 * @brief Simulates terrain diffusion due to coastal erosion.
 *
 * This function applies an iterative coastal erosion diffusion process on a
 * terrain elevation array (`z`), taking into account the presence and depth of
 * water. The erosion model smooths the terrain near shorelines while
 * maintaining constant water-surface height.
 *
 * At each iteration:
 *  - A local water mask is computed using `water_mask(water_depth, z,
 * additional_depth)`.
 *  - The terrain elevations are smoothed (diffused) using a masked Laplacian
 * filter, applied only near the water boundary.
 *  - The water depth is adjusted to preserve the total water surface height
 * (i.e., `z + water_depth` remains constant).
 *
 * This results in a realistic simulation of coastal erosion processes where the
 * terrain near the waterline is progressively smoothed and redistributed.
 *
 * @param z                Reference to the terrain elevation array (modified in
 *                         place).
 * @param water_depth      Reference to the array representing water depth
 *                         values (modified to preserve water surface level).
 * @param additional_depth Additional virtual water depth used to estimate the
 *                         influence region of the shoreline during mask
 *                         computation.
 * @param iterations       Number of erosion–diffusion iterations to apply.
 * @param p_water_mask     Optional output pointer. If non-null, it receives the
 *                         last computed water mask used during the final
 *                         iteration.
 *
 * **Example**
 * @include ex_coastal_erosion_diffusion.cpp
 *
 * **Result**
 * @image html ex_coastal_erosion_diffusion.png
 */
void coastal_erosion_diffusion(Array       &z,
                               Array       &water_depth,
                               float        additional_depth,
                               int          iterations = 10,
                               const Array *p_mask = nullptr,
                               Array       *p_water_mask = nullptr);

/**
 * @brief Applies a coastal erosion profile to a terrain elevation field.
 *
 * This function modifies the elevation array @p z by carving a coastal profile
 * at the interface between ground and water. It uses distance transforms both
 * from ground regions and from water regions to determine how far each point is
 * from the shoreline. According to this distance, the function applies:
 *
 * - A ground-side shore slope with an optional scarp transition.
 * - A water-side underwater slope ensuring continuity with the ground slope.
 * - An optional post-filtering step (Laplace smoothing) restricted to the
 * shoreline region.
 *
 * Water depth is adjusted after filtering to preserve the original water
 * surface height.
 *
 * Optionally, a shoreline mask can be returned describing the region influenced
 * by the coastal transformation.
 *
 * @param z                   Array of terrain elevations. Modified in-place.
 * @param water_depth         Array of water depths. Modified to preserve water
 *                            surface height after terrain changes.
 * @param shore_ground_extent Horizontal extent (in grid units) over which the
 *                            ground-side shore profile is applied.
 * @param shore_water_extent  Horizontal extent (in grid units) over which the
 *                            underwater profile is applied.
 * @param slope_shore         Ground-side slope magnitude of the coastal
 *                            profile. Expressed in elevation units per domain
 *                            width.
 * @param slope_shore_water   Water-side slope magnitude of the underwater
 *                            profile. Expressed similarly to @p slope_shore.
 * @param scarp_extent_ratio  Ratio defining the relative extent of the scarp
 *                            region. A value in [0,1]:
 *                             - 0 → no scarp, only slope
 *                             - 1 → all scarp
 * @param apply_post_filter   If true, applies Laplacian smoothing to @p z
 *                            restricted to shoreline areas.
 * @param p_shore_mask        Optional output pointer. If non-null, receives the
 *                            shoreline mask (values in [0,1]) indicating where
 *                            the coastal transformation was applied.
 *
 * **Example**
 * @include ex_coastal_erosion_profile.cpp
 *
 * **Result**
 * @image html ex_coastal_erosion_profile.png
 */
void coastal_erosion_profile(Array &z,
                             Array &water_depth,
                             float  shore_ground_extent, // pixels
                             float  shore_water_extent,
                             float  slope_shore = 0.5f,
                             float  slope_shore_water = 0.5f,
                             float  scarp_extent_ratio = 0.5f, // in [0, 1]
                             bool   apply_post_filter = true,
                             int    post_filter_iterations = 3,
                             bool   solid_shore_mask = true,
                             float  scarp_mask_transition_ratio = 0.2f,
                             const Array *p_noise = nullptr,
                             Array       *p_shore_mask = nullptr,
                             Array       *p_scarp_mask = nullptr);

void coastal_erosion_profile(Array       &z,
                             const Array *p_mask,
                             Array       &water_depth,
                             float        shore_ground_extent, // pixels
                             float        shore_water_extent,
                             float        slope_shore = 0.5f,
                             float        slope_shore_water = 0.5f,
                             float scarp_extent_ratio = 0.5f, // in [0, 1]
                             bool  apply_post_filter = true,
                             int   post_filter_iterations = 3,
                             bool  solid_shore_mask = true,
                             float scarp_mask_transition_ratio = 0.2f,
                             const Array *p_noise = nullptr,
                             Array       *p_shore_mask = nullptr,
                             Array       *p_scarp_mask = nullptr);

/**
 * @brief Fill the depressions of the heightmap using the Planchon-Darboux
 * algorithm.
 *
 * Fill heightmap depressions to ensure that every cell can be connected to the
 * boundaries following a downward slope @cite Planchon2002.
 *
 * @param z          Input array.
 * @param iterations Number of iterations.
 * @param epsilon
 *
 * **Example**
 * @include ex_depression_filling.cpp
 *
 * **Result**
 * @image html ex_depression_filling.png
 */
void depression_filling(Array &z, int iterations = 1000, float epsilon = 1e-4f);

void depression_filling_priority_flood(Array &z,
                                       bool   apply_post_filter = false);

/**
 * @brief
 *
 * @param z_before       Input array (before erosion).
 * @param z_after        Input array (after erosion).
 * @param erosion_map    Erosion map.
 * @param deposition_map Deposition map.
 * @param tolerance      Tolerance for erosion / deposition definition.
 *
 * **Example**
 * @include ex_erosions_maps.cpp
 *
 * **Result**
 * @image html ex_erosions_maps0.png
 * @image html ex_erosions_maps1.png
 */
void erosion_maps(Array &z_before,
                  Array &z_after,
                  Array &erosion_map,
                  Array &deposition_map,
                  float  tolerance = 0.f);

/**
 * @brief Generates a modified bedrock heightmap from an input elevation array.
 *
 * This function adjusts elevations based on overall height and local slope:
 * - Reduces values relative to the elevation range using `elevation_strength`.
 * - Reduces values in steep areas using `slope_strength` and `slope_limit`.
 *
 * @param  z                  Input elevation array.
 * @param  elevation_strength Strength of elevation-based adjustment.
 * @param  slope_strength     Strength of slope-based adjustment.
 * @param  slope_limit        Slope threshold for slope-based adjustment.
 * @param  zmin               Minimum elevation to consider (computed from z if
 *                            zmin > zmax).
 * @param  zmax               Maximum elevation to consider (computed from z if
 *                            zmin > zmax).
 * @return                    Array Modified bedrock heightmap.
 *
 *  **Example**
 * @include ex_hydraulic_particle.cpp
 *
 * **Result**
 * @image html ex_hydraulic_particle.png
 */
Array generate_bedrock(const Array &z,
                       float        elevation_strength,
                       float        slope_strength,
                       float        slope_limit,
                       float        zmin = 0.f,
                       float        zmax = -1.f);

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

void hydraulic_algebric(Array &z,
                        float  talus_ref,
                        int    ir,
                        Array *p_bedrock = nullptr,
                        Array *p_erosion_map = nullptr,
                        Array *p_deposition_map = nullptr,
                        float  c_erosion = 0.05f,
                        float  c_deposition = 0.05f,
                        int    iterations = 1); ///< @overload

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
                     float  rain_rate = 0.5f); ///< @overload

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

void hydraulic_musgrave(Array &z,
                        int    iterations = 100,
                        float  c_capacity = 1.f,
                        float  c_erosion = 0.1f,
                        float  c_deposition = 0.1f,
                        float  water_level = 0.01f,
                        float  evap_rate = 0.01f); ///< @overload

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
 * {@link flow_accumulation_dinf}).
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
void hydraulic_stream(Array &z,
                      float  c_erosion,
                      float  talus_ref,
                      Array *p_bedrock = nullptr,
                      Array *p_moisture_map = nullptr,
                      Array *p_erosion_map = nullptr, // -> out
                      int    ir = 1,
                      float  clipping_ratio = 10.f);

void hydraulic_stream(Array &z,
                      Array *p_mask,
                      float  c_erosion,
                      float  talus_ref,
                      Array *p_bedrock = nullptr,
                      Array *p_moisture_map = nullptr,
                      Array *p_erosion_map = nullptr, // -> out
                      int    ir = 1,
                      float  clipping_ratio = 10.f); ///< @overload

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
 * {@link flow_accumulation_dinf}).
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
                          Array       *p_flow_map = nullptr); ///< @overload

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
void hydraulic_stream_upscale_amplification(
    Array &z,
    Array *p_mask,
    float  c_erosion,
    float  talus_ref,
    int    upscaling_levels = 1,
    float  persistence = 1.f,
    int    ir = 1,
    float  clipping_ratio = 10.f); ///< @overload

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

void deposition_fill_holes(Array       &z,
                           int          deposition_ir,
                           float        deposition_strength,
                           const Array *p_mask,
                           int          iterations = 1); ///< @overload

/**
 * @brief Simulates hydraulic erosion on a heightmap using particle-based flow.
 *
 * Particles traverse the heightmap `z`, eroding and depositing material
 * according to local slope, capacity, inertia, and optional directional bias.
 *
 * @param z                       Heightmap array to modify.
 * @param nparticles              Number of erosion particles to simulate.
 * @param seed                    Random seed for particle initialization.
 * @param p_bedrock               Optional bedrock array to limit erosion.
 * @param p_moisture_map          Optional moisture map affecting
 *                                erosion/deposition.
 * @param p_elevation_shift       Optional elevation shift map.
 * @param p_erosion_map           Optional output array recording total erosion.
 * @param p_deposition_map        Optional output array recording deposition.
 * @param c_capacity              Sediment capacity of each particle.
 * @param c_erosion               Erosion rate coefficient.
 * @param c_deposition            Deposition rate coefficient.
 * @param c_inertia               Particle inertia factor.
 * @param c_gravity               Gravity effect on particle movement.
 * @param drag_rate               Particle velocity damping.
 * @param evap_rate               Sediment evaporation rate.
 * @param enable_directional_bias Enable bias along a preferred slope direction.
 * @param angle_bias              Bias angle in degrees (if directional bias
 *                                enabled).
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
                        float         c_inertia = 0.1f,
                        float         c_gravity = 1.f,
                        float         drag_rate = 0.001f,
                        float         evap_rate = 0.001f,
                        bool          enable_directional_bias = false,
                        float         angle_bias = 30.f);

/*! @brief See hmap::hydraulic_particle */
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
                        float         c_inertia = 0.1f,
                        float         c_gravity = 1.f,
                        float         drag_rate = 0.001f,
                        float         evap_rate = 0.001f,
                        bool          enable_directional_bias = false,
                        float         angle_bias = 30.f);

/**
 * @brief Particle-based hydraulic erosion with flow-field coupling (McDonald's
 * model): persistent per-cell discharge and momentum fields, exponentially
 * filtered across iterations, couple particles through the mean local flow,
 * producing coherent drainage networks; a bank-stability debris flow runs in
 * the same solver loop against a separate sediment layer. Clean-room port of
 * erosiv/soillib (LGPL-3, Nicholas McDonald),
 * https://github.com/erosiv/soillib.
 *
 * Parameters are physical: the terrain is interpreted as a world_extent_km x
 * world_extent_km domain with height range z_scale_km. This makes behavior
 * consistent across resolutions.
 *
 * @warning Single-resolution runs numerically diverge at high resolutions (>=
 * 1024^2) beyond a few hundred steps with default parameters — use
 * hydraulic_mcdonald_multiscale for high-resolution terrain.
 *
 * Results are reproducible up to atomic scheduling order (particles read fields
 * that other particles concurrently modify via atomic adds — a property shared
 * with the reference implementation).
 *
 * @param z               Input/output heightmap. In: bedrock. Out: bedrock plus
 *                        sediment (total surface).
 * @param steps           Number of erosion iterations.
 * @param seed            Random seed number.
 * @param p_sediment_map  Optional output: final sediment layer.
 * @param p_discharge_map Optional output: water discharge field (usable as a
 *                        river / water mask).
 * @param world_extent_km Physical domain edge length [km].
 * @param z_scale_km      Physical height range of z's [0, 1] span [km].
 * @param samples         Particles per iteration.
 * @param maxage          Maximum particle age (steps along a trajectory).
 * @param lrate           Exponential filter rate for the flow fields.
 * @param time_step       Geological timestep [y].
 * @param rainfall        Rainfall rate [m/y].
 * @param evap_rate       Evaporation rate.
 * @param gravity         Gravity [m/s^2].
 * @param viscosity       Kinematic viscosity (flow-coupling strength).
 * @param bed_shear       Bed shear coefficient.
 * @param crit_slope      Critical (stable bank) slope [m/m].
 * @param settle_rate     Debris settling rate.
 * @param thermal_rate    Thermal erosion rate.
 * @param deposition_rate Fluvial deposition rate.
 * @param suspension_rate Fluvial suspension rate.
 * @param exit_slope      Boundary exit slope [m/m].
 *
 * **Example**
 * @include ex_hydraulic_mcdonald.cpp
 *
 * **Result**
 * @image html ex_hydraulic_mcdonald0.png
 * @image html ex_hydraulic_mcdonald1.png
 * @image html ex_hydraulic_mcdonald2.png
 */
void hydraulic_mcdonald(Array        &z,
                        int           steps,
                        std::uint32_t seed,
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
 * @brief Multiscale driver for hydraulic_mcdonald: erodes on a halving
 * resolution ladder derived from z.shape (coarsest first), resampling the full
 * model state (bedrock, sediment, discharge, momentum) between levels. The
 * numerically stable route to high-resolution erosion with this model.
 * {512, 256, 128} on a 1024^2 input runs 256^2 (512 steps), 512^2 (256), then
 * 1024^2 (128). See hydraulic_mcdonald for the model and parameters.
 *
 * **Example**
 * @include ex_hydraulic_mcdonald.cpp
 *
 * **Result**
 * @image html ex_hydraulic_mcdonald3.png
 */
void hydraulic_mcdonald_multiscale(
    Array                  &z,
    std::uint32_t           seed,
    const std::vector<int> &steps_per_level = {512, 256, 128},
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
 * @brief Apply phase-guided hydraulic procedural erosion to a heightmap.
 *
 * Computes a phase field from the input terrain and generates ridge-aligned
 * erosion patterns modulated by flow accumulation and local gradient.
 * Optionally applies procedural noise, deposition, and outputs a ridge mask.
 *
 * @param[in,out] z                         Heightmap to erode (modified in
 *                                          place).
 * @param         kp_global                 Global kernel size controlling
 *                                          feature scale.
 * @param         c_erosion                 Global erosion intensity.
 * @param         seed                      Random seed for phase/noise
 *                                          generation.
 * @param         erosion_profile           Ridge shaping profile type.
 * @param         erosion_profile_parameter Parameter controlling the profile
 *                                          shape.
 * @param         angle_shift               Global phase angle offset (degrees).
 * @param         phase_smoothing           Controls smoothing between phase
 *                                          transitions.
 * @param         talus_ref                 Reference talus angle for flow
 *                                          accumulation.
 * @param         gradient_scaling_ratio    Blend factor for gradient-based
 *                                          modulation.
 * @param         gradient_power            Exponent applied to normalized
 *                                          gradient.
 * @param         apply_deposition          If true, applies post-erosion
 *                                          deposition pass.
 * @param         enable_default_noise      If true, generates internal FBM
 *                                          noise fields.
 * @param         noise_amp                 Amplitude of procedural noise.
 * @param         p_kp_multiplier           Optional spatial modulation of
 *                                          kernel size.
 * @param         p_angle_shift             Optional spatial phase shift.
 * @param         p_noise_x                 Optional external noise field (x
 *                                          component).
 * @param         p_noise_y                 Optional external noise field (y
 *                                          component).
 * @param         p_ridge_mask              Optional output ridge mask (cosine
 *                                          of phase).
 * @param         bbox                      World-space bounding box.
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
 * Applies the erosion operator across multiple octaves with
 * lacunarity/persistence scaling, accumulating a normalized ridge mask.
 *
 * @param[in,out] z                         Heightmap to erode (modified in
 *                                          place).
 * @param         kp_global                 Base kernel size.
 * @param         c_erosion                 Base erosion intensity.
 * @param         seed                      Random seed.
 * @param         erosion_profile           Ridge shaping profile.
 * @param         octaves                   Number of fBm octaves.
 * @param         persistence               Amplitude decay per octave.
 * @param         lacunarity                Frequency growth per octave.
 * @param         erosion_profile_parameter Profile parameter.
 * @param         angle_shift               Global phase shift (degrees).
 * @param         phase_smoothing           Phase smoothing factor.
 * @param         talus_ref                 Reference talus angle.
 * @param         gradient_scaling_ratio    Gradient modulation ratio.
 * @param         gradient_power            Gradient exponent.
 * @param         apply_deposition          Enable deposition pass.
 * @param         enable_default_noise      Enable internal noise generation.
 * @param         noise_amp                 Noise amplitude.
 * @param         p_kp_multiplier           Optional spatial kernel modulation.
 * @param         p_angle_shift             Optional spatial phase shift.
 * @param         p_noise_x                 Optional external noise (x).
 * @param         p_noise_y                 Optional external noise (y).
 * @param         p_ridge_mask              Optional accumulated ridge mask
 *                                          output.
 * @param         bbox                      World-space bounding box.
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
 * @brief Simulates hydraulic erosion and deposition on a heightmap using the
 * Schott method.
 *
 * This function performs hydraulic erosion on the given heightmap `z` over a
 * specified number of iterations. It includes parameters for controlling
 * erosion, deposition, and flow routing. Optional flow accumulation can also be
 * computed and stored in the `p_flow` array.
 *
 * @param[in,out] z                      The heightmap array to be modified.
 *                                       Heights are updated in-place.
 * @param[in]     iterations             The number of iterations for the
 *                                       hydraulic erosion process.
 * @param[in]     talus                  An array defining the slope threshold
 *                                       for erosion.
 * @param[in]     c_erosion              Erosion coefficient (default: 1.0f).
 * @param[in]     c_thermal              Thermal erosion coefficient (default:
 *                                       0.1f).
 * @param[in]     c_deposition           Deposition coefficient (default:
 * 0.2f).
 * @param[in]     flow_acc_exponent      Exponent controlling the influence of
 *                                       flow accumulation on erosion (default:
 *                                       0.8f).
 * @param[in]     flow_acc_exponent_depo Exponent controlling the influence of
 *                                       flow accumulation on deposition
 *                                       (default: 0.8f).
 * @param[in]     flow_routing_exponent  Exponent controlling flow routing
 *                                       behavior (default: 1.3f).
 * @param[in]     thermal_weight         Weight of thermal erosion effects
 *                                       (default: 1.5f).
 * @param[in]     deposition_weight      Weight of deposition effects (default:
 *                                       2.5f).
 * @param[out]    p_flow                 Optional pointer to an array for
 *                                       storing flow accumulation data. If
 *                                       null, flow data is not returned
 *                                       (default: nullptr).
 *
 * @note Taken from https://hal.science/hal-04565030v1/document
 *
 * @note Only available if OpenCL is enabled.
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

void hydraulic_schott(Array       &z,
                      int          iterations,
                      const Array &talus,
                      Array       *p_mask,
                      float        c_erosion = 1.f,
                      float        c_thermal = 0.1f,
                      float        c_deposition = 0.2f,
                      float        flow_acc_exponent = 0.8f,
                      float        flow_acc_exponent_depo = 0.8f,
                      float        flow_routing_exponent = 1.3f,
                      float        thermal_weight = 1.5f,
                      float        deposition_weight = 2.5f,
                      Array       *p_flow = nullptr); ///< @overload

/*! @brief See hmap::hydraulic_schott */
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
                          Array       *p_flow_map = nullptr); ///< @overload

/*! @brief See hmap::hydraulic_vpipes */
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

/**
 * @brief Simulate a mudslide (landslide-driven material redistribution) on a
 * height field.
 *
 * Applies an iterative viscous flow process to the height array @p z,
 * constrained by @p landslide_mask and limited by the specified @p depth. The
 * flow behavior is controlled by a depth exponent and a viscosity power law.
 * Optionally outputs the transported depth map.
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
 * @param[out]    p_depth             Optional pointer to store resulting depth
 *                                    map (nullable).
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

/*! @brief See hmap::mudslide */
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
 * This function modifies the given heightmap by introducing linear or radial
 * "rift-like" noise patterns. The deformation can be controlled by several
 * parameters such as direction, amplitude, noise shifts, and optional external
 * noise arrays. Optionally, the effect can be masked using a power-based
 * blending function.
 *
 * @param z                     Reference to the heightmap array to be modified
 *                              in-place.
 * @param kw                    Frequency vector (kx, ky) scaling the
 *                              deformation in X and Y directions.
 * @param angle                 Orientation of the rift in degrees (0° =
 *                              horizontal, 90° = vertical).
 * @param amplitude             Strength of the rift deformation applied to the
 *                              heightmap.
 * @param seed                  Random seed used for deterministic noise
 *                              generation.
 * @param elevation_noise_shift Vertical offset applied to the base noise to
 *                              shift elevation influence.
 * @param k_smooth_bottom       Lower smoothing factor for Voronoi-based noise
 *                              computation.
 * @param k_smooth_top          Upper smoothing factor for Voronoi-based noise
 *                              computation.
 * @param radial_spread_amp     Amplitude controlling radial spreading away from
 *                              the rift axis.
 * @param elevation_noise_amp   Amplitude scaling the influence of the
 *                              heightmap's initial values as noise input.
 * @param clamp_vmin            Minimum clamp value for the Voronoi noise before
 *                              remapping.
 * @param remap_vmin            Minimum remap value for scaling noise output.
 * @param apply_mask            If true, applies a power-based blending mask
 *                              instead of a direct overwrite.
 * @param mask_gamma            Gamma exponent used when applying the mask to
 *                              control blending.
 * @param p_noise_x             Optional pointer to an external noise array for
 *                              X-offset perturbation (nullptr if unused).
 * @param p_noise_y             Optional pointer to an external noise array for
 *                              Y-offset perturbation (nullptr if unused).
 * @param center                2D vector specifying the central point around
 *                              which the rift effect is computed.
 * @param bbox                  Bounding box (xmin, xmax, ymin, ymax) defining
 *                              the spatial domain of the heightmap.
 *
 * **Example**
 * @include ex_rifts.cpp
 *
 * **Result**
 * @image html ex_rifts.png
 */
void rifts(Array           &z,
           const glm::vec2 &kw,    //  = {4.f, 1.2f},
           float            angle, // degs
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
 * @todo deposition map
 *
 * @param z                             Input array.
 * @param p_mask                        Intensity mask, expected in [0, 1]
 *                                      (applied as a post-processing).
 * @param talus                         Talus limit.
 * @param p_deposition_map              [out] Reference to the deposition map,
 *                                      provided as an output field.
 * @param max_deposition                Maximum height of sediment deposition.
 * @param iterations                    Number of iterations.
 * @param thermal_erosion_subiterations Number of thermal erosion iterations for
 *                                      each pass.
 *
 * **Example**
 * @include ex_sediment_deposition.cpp
 *
 * **Result**
 * @image html ex_sediment_deposition.png
 */
void sediment_deposition(Array       &z,
                         Array       *p_mask,
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
 * Computes a sediment (deposition) map using thermal erosion constrained by
 * per-cell talus limits. Optionally smooths the result and blends it back into
 * the heightmap.
 *
 * @param[in,out] z                 Heightmap to modify.
 * @param[in]     talus_layer       Base talus values for erosion.
 * @param[in]     talus_upper_limit Upper talus threshold used for masking.
 * @param[in]     iterations        Number of erosion iterations.
 * @param[in]     apply_post_filter Enable post smoothing and blending.
 * @param[out]    p_deposition_map  Optional output deposition map (may be
 *                                  null).
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
 * This function modifies the input heightfield `z` by simulating geological
 * strata patterns. It combines directional shifts, fractal noise, and
 * ridge-based perturbations to produce layered structures in the data. The MUST
 * BE NORMALIZED in [0, 1].
 *
 * @param z                 Reference to the heightfield array to modify, MUST
 *                          BE NORMALIZED in [0, 1].
 * @param angle             Horizontal orientation of the strata in degrees.
 * @param slope             Vertical slope of the strata.
 * @param gamma             Gamma exponent for the non-linear remapping (e.g.,
 *                          0.5 for smoothing, 1.5 for sharpening).
 * @param seed              Seed for deterministic noise generation.
 * @param linear_gamma      If true, applies sharp linear gamma mapping; if
 *                          false, uses smooth gamma mapping.
 * @param kz                Base scaling factor for the stratification
 *                          frequency.
 * @param octaves           Number of iterative stratification passes
 *                          (multiscale detail).
 * @param lacunarity        Frequency multiplier applied at each octave for
 *                          fractal scaling.
 * @param gamma_noise_ratio Ratio controlling how noise influences gamma
 *                          variation (0 = no noise, 1 = full influence).
 * @param noise_amp         Amplitude of the base Perlin noise used to modulate
 *                          the strata.
 * @param noise_kw          Frequency vector for the base Perlin noise along X
 *                          and Y axes.
 * @param ridge_noise_kw    Frequency vector for the Voronoi ridge noise (x =
 *                          main frequency, y = directional frequency).
 * @param ridge_angle_shift Additional angular shift (in degrees) for the ridge
 *                          direction, relative to `angle`.
 * @param ridge_noise_amp   Amplitude of the ridge noise modulation.
 * @param ridge_clamp_vmin  Minimum clamp value for ridge noise response.
 * @param ridge_remap_vmin  Minimum remap value for ridge modulation (used for
 *                          reverse remapping).
 * @param apply_mask        If true, applied an elevation mask on the effect.
 * @param mask_gamma        Gamma applied to the mask used for blending original
 *                          and stratified values for the elevation mask.
 * @param p_mask            Optional filter mask, expected in the range [0, 1].
 * @param bbox              Bounding box of the domain as `{xmin, xmax, ymin,
 *                          ymax}`.
 *
 * @note
 * - Setting `linear_gamma` to `false` produces smoother transitions, while
 * `true` creates sharper layer boundaries.
 * - Increasing `octaves` adds multiscale detail but also increases
 * computational cost.
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
            float            gamma, // e.g 0.5f or 1.5f
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
 * Generates layered cell-like structures with optional directional distortion
 * and noise-based perturbation.
 *
 * @param z                     Heightmap array to modify.
 * @param kw                    Cell frequency (x,y).
 * @param amp                   Displacement amplitude.
 * @param seed                  Random seed.
 * @param gamma                 Vertical shaping exponent.
 * @param gamma_lateral         Lateral shaping exponent.
 * @param angle                 Layer orientation angle (degrees).
 * @param noise_amp             Amplitude of input distortion noise.
 * @param absolute_displacement Use absolute instead of relative displacement.
 * @param occurence_probability Probability of cell activation.
 * @param p_noise_x             Optional X distortion field.
 * @param p_noise_y             Optional Y distortion field.
 * @param bbox                  Domain bounding box.
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

/**
 * @brief Applies stratified cell displacement with optional masking.
 *
 * If `p_mask` is provided, the effect is blended with the original heightmap
 * using the mask as interpolation weight.
 *
 * @param z                     Heightmap array to modify.
 * @param kw                    Cell frequency (x,y).
 * @param amp                   Displacement amplitude.
 * @param seed                  Random seed.
 * @param p_mask                Optional mask controlling spatial blending.
 * @param gamma                 Vertical shaping exponent.
 * @param gamma_lateral         Lateral shaping exponent.
 * @param angle                 Layer orientation angle (degrees).
 * @param noise_amp             Amplitude of input distortion noise.
 * @param absolute_displacement Use absolute instead of relative displacement.
 * @param occurence_probability Probability of cell activation.
 * @param p_noise_x             Optional X distortion field.
 * @param p_noise_y             Optional Y distortion field.
 * @param bbox                  Domain bounding box.
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
 * Accumulates several scaled `strata_cells` layers using persistence and
 * lacunarity. Can optionally generate internal distortion noise.
 *
 * @param z                     Heightmap array to modify.
 * @param kw                    Base cell frequency (x,y).
 * @param amp                   Base displacement amplitude.
 * @param seed                  Random seed.
 * @param gamma                 Vertical shaping exponent.
 * @param gamma_lateral         Lateral shaping exponent.
 * @param angle                 Layer orientation angle (degrees).
 * @param enable_default_noise  Generate internal distortion noise.
 * @param default_noise_amp     Amplitude of generated noise.
 * @param absolute_displacement Use absolute instead of relative displacement.
 * @param occurence_probability Probability of cell activation.
 * @param octaves               Number of fBm layers.
 * @param persistence           Amplitude multiplier per octave.
 * @param lacunarity            Frequency multiplier per octave.
 * @param p_noise_x             Optional X distortion field.
 * @param p_noise_y             Optional Y distortion field.
 * @param bbox                  Domain bounding box.
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

/**
 * @brief Applies masked multi-octave stratified cell displacement.
 *
 * Works like `strata_cells_fbm` but blends the result with the original
 * heightmap using `p_mask`.
 *
 * @param z                     Heightmap array to modify.
 * @param kw                    Base cell frequency (x,y).
 * @param amp                   Base displacement amplitude.
 * @param seed                  Random seed.
 * @param p_mask                Optional mask controlling spatial blending.
 * @param gamma                 Vertical shaping exponent.
 * @param gamma_lateral         Lateral shaping exponent.
 * @param angle                 Layer orientation angle (degrees).
 * @param enable_default_noise  Generate internal distortion noise.
 * @param default_noise_amp     Amplitude of generated noise.
 * @param absolute_displacement Use absolute instead of relative displacement.
 * @param occurence_probability Probability of cell activation.
 * @param octaves               Number of fBm layers.
 * @param persistence           Amplitude multiplier per octave.
 * @param lacunarity            Frequency multiplier per octave.
 * @param p_noise_x             Optional X distortion field.
 * @param p_noise_y             Optional Y distortion field.
 * @param bbox                  Domain bounding box.
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
 * Projects the input array along several directions (sequential or random),
 * blends the results, and optionally applies a spatial warp.
 *
 * @param[in,out] z                 Input/output height array.
 * @param[in]     talus             Talus values used for projection.
 * @param[in]     direction_offset  Offset for sequential direction selection.
 * @param[in]     direction_count   Number of directions to use (max 8).
 * @param[in]     random_directions Whether to use random directions.
 * @param[in]     seed              Seed for random direction selection.
 * @param[in]     mix_ratio         Blending factor between iterations.
 * @param[in]     p_mask            Optional mask array.
 * @param[in]     p_dx              Optional X displacement map for warp.
 * @param[in]     p_dy              Optional Y displacement map for warp.
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
 * Quantizes elevation into kz strata levels and shapes terrace transitions
 * using a gamma profile. Optional noise can modulate the profile.
 *
 * @param z                 Heightmap modified in place.
 * @param gamma             Terrace profile exponent (e.g. 0.5 = sharper, 1.5 =
 *                          smoother).
 * @param seed              Random seed used for profile variation.
 * @param kz                Number of strata levels.
 * @param linear_gamma      If true, uses linear terrace transitions.
 * @param gamma_noise_ratio Influence of noise on the terrace profile.
 * @param slope             Width of the inclined terrace face in [0, 1] (0 =
 *                          vertical step, 1 = fully sloped transition).
 * @param angle             Direction of the terrace inclination in radians.
 * @param p_noise           Optional noise map for profile modulation.
 * @param bbox              Bounding box used to evaluate the terrace direction.
 *
 * **Example**
 * @include ex_strata_terrace.cpp
 *
 * **Result**
 * @image html ex_strata_terrace.png
 */
void strata_terrace(Array        &z,
                    float         gamma, // e.g 0.5f or 1.5f
                    std::uint32_t seed,
                    float         kz = 4.f, // 4-layers
                    bool          linear_gamma = true,
                    float         gamma_noise_ratio = 0.5f,
                    float         slope = 0.f,
                    float         angle = 0.f,
                    const Array  *p_noise = nullptr,
                    glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Applies a masked terrace (stratification) filter to a heightmap.
 *
 * Same as strata_terrace(), but only affects regions defined by p_mask.
 *
 * @param z                 Heightmap modified in place.
 * @param gamma             Terrace profile exponent.
 * @param seed              Random seed used for profile variation.
 * @param p_mask            Optional mask controlling where terraces are
 *                          applied.
 * @param kz                Number of strata levels.
 * @param linear_gamma      If true, uses linear terrace transitions.
 * @param gamma_noise_ratio Influence of noise on the terrace profile.
 * @param slope             Width of the inclined terrace face in [0, 1] (0 =
 *                          vertical step, 1 = fully sloped transition).
 * @param angle             Direction of the terrace inclination in radians.
 * @param p_noise           Optional noise map for profile modulation.
 * @param bbox              Bounding box used to evaluate the terrace direction.
 *
 * **Example**
 * @include ex_strata_terrace.cpp
 *
 * **Result**
 * @image html ex_strata_terrace.png
 */
void strata_terrace(Array        &z,
                    float         gamma, // e.g 0.5f or 1.5f
                    std::uint32_t seed,
                    const Array  *p_mask,
                    float         kz = 4.f, // 4-layers
                    bool          linear_gamma = true,
                    float         gamma_noise_ratio = 0.5f,
                    float         slope = 0.f,
                    float         angle = 0.f,
                    const Array  *p_noise = nullptr,
                    glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Apply thermal weathering erosion.
 *
 * Based on https://www.shadertoy.com/view/XtKSWh
 *
 * @param z                Input array.
 * @param p_mask           Filter mask, expected in [0, 1].
 * @param talus            Talus limit.
 * @param p_bedrock        Lower elevation limit.
 * @param p_deposition_map [out] Reference to the deposition map, provided as an
 *                         output field.
 * @param iterations       Number of iterations.
 *
 * **Example**
 * @include ex_thermal.cpp
 *
 * **Result**
 * @image html ex_thermal.png
 */
void thermal(Array       &z,
             const Array &talus,
             int          iterations = 10,
             Array       *p_bedrock = nullptr,
             Array       *p_deposition_map = nullptr);

void thermal(Array       &z,
             const Array *p_mask,
             const Array &talus,
             int          iterations = 10,
             Array       *p_bedrock = nullptr,
             Array       *p_deposition_map = nullptr); ///< @overload

void thermal(Array &z,
             float  talus,
             int    iterations = 10,
             Array *p_bedrock = nullptr,
             Array *p_deposition_map = nullptr); ///< @overload

/**
 * @brief Apply thermal weathering erosion with automatic determination of the
 * bedrock.
 *
 * @todo more comprehensive documentation on auto-bedrock algo.
 * @todo fix hard-coded parameters.
 *
 * @see                    {@link thermal}
 *
 * @param z                Input array.
 * @param talus            Local talus limit.
 * @param iterations       Number of iterations.
 * @param p_deposition_map [out] Reference to the deposition map, provided as an
 *                         output field.
 *
 * **Example**
 * @include ex_thermal_auto_bedrock.cpp
 *
 * **Result**
 * @image html ex_thermal_auto_bedrock.png
 */
void thermal_auto_bedrock(Array       &z,
                          const Array &talus,
                          int          iterations = 10,
                          Array       *p_deposition_map = nullptr);

void thermal_auto_bedrock(Array       &z,
                          const Array *p_mask,
                          const Array &talus,
                          int          iterations = 10,
                          Array *p_deposition_map = nullptr); ///< @overload

void thermal_auto_bedrock(Array &z,
                          float,
                          int    iterations = 10,
                          Array *p_deposition_map = nullptr); ///< @overload

/**
 * @brief Apply iterative thermal flattening erosion on a heightmap.
 *
 * Material is redistributed locally according to the steepest neighbor slope
 * and a per-cell talus threshold. Gentle and steep slopes can be relaxed
 * independently using separate coefficients.
 *
 * The computation is performed on the GPU using OpenCL.
 *
 * @param[in,out] z          Heightmap modified in-place.
 * @param[in]     talus      Per-cell talus threshold.
 * @param[in]     iterations Number of iterations.
 * @param[in]     sigma_inf  Relaxation factor for slopes below talus.
 * @param[in]     sigma_sup  Relaxation factor for slopes above talus.
 */
void thermal_flatten(Array       &z,
                     const Array &talus,
                     int          iterations,
                     float        sigma_inf = 0.5f,
                     float        sigma_sup = 0.f);

void thermal_flatten(Array       &z,
                     const Array *p_mask,
                     const Array &talus,
                     int          iterations,
                     float        sigma_inf = 0.5f,
                     float        sigma_sup = 0.f); ///< @overload

/**
 * @brief Apply thermal weathering erosion.
 *
 * Based on averaging over first neighbors, see \cite Olsen2004. Refer to \cite
 * Musgrave1989 for the original formulation.
 *
 * Thermal erosion refers to the process in which surface sediment weakens due
 * to temperature and detaches, falling down the slopes of the terrain until a
 * resting place is reached, where smooth plateaus tend to form
 * @cite Musgrave1989.
 *
 * @param z                Input array.
 * @param p_mask           Filter mask, expected in [0, 1].
 * @param talus            Talus limit.
 * @param p_bedrock        Lower elevation limit.
 * @param p_deposition_map [out] Reference to the deposition map, provided as an
 *                         output field.
 * @param iterations       Number of iterations.
 */
void thermal_olsen(Array &z, const Array &talus, int iterations);

void thermal_olsen(Array       &z,
                   const Array *p_mask,
                   const Array &talus,
                   int          iterations); ///< @overload

/**
 * @brief Apply thermal weathering erosion to give a scree like effect.
 *
 * @note Only available if OpenCL is enabled.
 *
 * @param z                Input array.
 * @param talus            Talus limit.
 * @param p_deposition_map [out] Reference to the deposition map, provided as an
 *                         output field.
 * @param iterations       Number of iterations.
 *
 * **Example**
 * @include ex_thermal_ridge.cpp
 *
 * **Result**
 * @image html ex_thermal_ridge.png
 */
void thermal_inflate(Array &z, const Array &talus, int iterations = 10);

void thermal_inflate(Array       &z,
                     const Array *p_mask,
                     const Array &talus,
                     int          iterations = 10); ///< @overload

/**
 * @brief Apply thermal erosion using a 'rib' algorithm (taken from Geomorph).
 * @param z          Input heightmap.
 * @param iterations Number of iterations.
 * @param p_bedrock  Lower elevation limit.
 *
 * **Example**
 * @include ex_thermal_rib.cpp
 *
 * **Result**
 * @image html ex_thermal_rib.png
 */
void thermal_rib(Array &z, int iterations);

void thermal_rib(Array &z, const Array *p_mask, int iterations); ///< @overload

/**
 * @brief Apply thermal weathering erosion to give a ridge like effect.
 *
 * @note Based on https://www.fractal-landscapes.co.uk/maths.html
 *
 * @note Only available if OpenCL is enabled.
 *
 * @param z                Input array.
 * @param talus            Talus limit.
 * @param p_deposition_map [out] Reference to the deposition map, provided as an
 *                         output field.
 * @param iterations       Number of iterations.
 *
 * **Example**
 * @include ex_thermal_ridge.cpp
 *
 * **Result**
 * @image html ex_thermal_ridge.png
 */
void thermal_ridge(Array       &z,
                   const Array &talus,
                   int          iterations = 10,
                   Array       *p_deposition_map = nullptr);

void thermal_ridge(Array       &z,
                   const Array *p_mask,
                   const Array &talus,
                   int          iterations = 10,
                   Array       *p_deposition_map = nullptr); ///< @overload

/**
 * @brief Applies the thermal erosion process with a uniform slope threshold.
 *
 * This overload of `thermal_schott` applies thermal erosion with a uniform
 * threshold value. It generates an internal talus map using the specified
 * constant `talus` value and applies the erosion process in the same manner as
 * the other overload.
 *
 * @param z          Reference to the array of elevation values that will be
 *                   modified.
 * @param talus      Constant threshold slope value used for all cells.
 * @param iterations Number of erosion iterations to apply.
 * @param intensity  Intensity factor controlling the amount of change per
 *                   iteration.
 *
 *  **Example**
 * @include ex_thermal_schott.cpp
 *
 * **Result**
 * @image html ex_thermal_schott.png
 */
void thermal_schott(Array       &z,
                    const Array &talus,
                    int          iterations = 10,
                    float        intensity = 0.2f,
                    Array       *p_deposition_map = nullptr);

void thermal_schott(Array       &z,
                    const Array *p_mask,
                    const Array &talus,
                    int          iterations = 10,
                    float        intensity = 0.2f,
                    Array       *p_deposition_map = nullptr); ///< @overload

/**
 * @brief Performs thermal scree erosion on a heightmap.
 *
 * This function applies a thermal erosion process that redistributes material
 * from steeper slopes to flatter areas, simulating talus formation. The process
 * iterates a given number of times to achieve a more stable terrain profile.
 *
 * @param[out] z                The heightmap to be modified in-place by the
 *                              erosion process.
 * @param[in]  talus            The threshold slope angles that determine where
 *                              material is moved.
 * @param[in]  zmax             The maximum allowed elevation for erosion
 *                              effects.
 * @param[in]  iterations       The number of erosion iterations to apply
 *                              (default:
 * 10).
 * @param[out] p_deposition_map Optional pointer to an array that stores the
 *                              deposited material per cell.
 */
void thermal_scree(Array       &z,
                   const Array &talus,
                   const Array &zmax,
                   int          iterations = 10,
                   Array       *p_deposition_map = nullptr);

void thermal_scree(Array       &z,
                   const Array *p_mask,
                   const Array &talus,
                   const Array &zmax,
                   int          iterations = 10,
                   Array       *p_deposition_map = nullptr); ///< @overload

/**
 * @brief Fill valleys using thermal scree deposition and height-based blending.
 *
 * Applies an erosion-based fill to valleys, then blends the result with the
 * original heightmap using a gamma-shaped mask derived from elevation.
 *
 * @param z          Heightmap to modify in place.
 * @param talus      Talus angle map controlling scree deposition.
 * @param iterations Number of erosion iterations.
 * @param gamma      Gamma applied to the height-based mixing mask.
 * @param ratio      Blend ratio controlling valley influence.
 * @param zmin       Minimum height for normalization (auto if zmax <= zmin).
 * @param zmax       Maximum height for normalization (auto if zmax <= zmin).
 *
 *  **Example**
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
 * Drainage basins are computed using D8 flow routing. For each basin, a
 * distance transform to the basin boundary is evaluated and used to lower
 * elevations near watershed divides, forming ridge lines.
 *
 * @param  z            Input elevation field.
 * @param  amplitude    Ridge carving strength.
 * @param  ir           Smoothing.
 * @param  edt_exponent Exponent applied to the distance field to control ridge
 *                      sharpness.
 *
 * @return              Elevation field with watershed ridges emphasized.
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
    float               width = 32.f, // pixels
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
