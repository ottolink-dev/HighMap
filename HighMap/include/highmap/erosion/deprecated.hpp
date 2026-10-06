/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once

#include <cstdint>

#include "highmap/array.hpp"

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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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

// --- Deprecated GPU Hydraulic Functions

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
// this function is going to be removed at some point
[[deprecated("This function is going to be removed at some point.")]]
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
