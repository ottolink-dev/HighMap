/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once

#include "highmap/array.hpp"

namespace hmap::gpu
{

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
             const Array *p_bedrock = nullptr,
             Array       *p_deposition_map = nullptr);

void thermal(Array       &z,
             const Array *p_mask,
             const Array &talus,
             int          iterations = 10,
             const Array *p_bedrock = nullptr,
             Array       *p_deposition_map = nullptr); ///< @overload

void thermal(Array       &z,
             float        talus,
             int          iterations = 10,
             const Array *p_bedrock = nullptr,
             Array       *p_deposition_map = nullptr); ///< @overload

/**
 * @brief Apply thermal weathering erosion with automatic determination of the
 * bedrock.
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
                          Array       *p_deposition_map = nullptr); ///<
// @overload

void thermal_auto_bedrock(Array &z,
                          float,
                          int    iterations = 10,
                          Array *p_deposition_map = nullptr); ///< @overload

/**
 * @brief Apply mass-preserving thermal weathering erosion.
 *
 * @param z                Input array.
 * @param talus            Talus limit.
 * @param iterations       Number of iterations.
 * @param rate             Erosion rate factor.
 * @param p_bedrock        Lower elevation limit.
 * @param p_deposition_map [out] Reference to the deposition map, provided as an
 *                         output field.
 */
void thermal_conserve(Array       &z,
                      const Array &talus,
                      int          iterations = 10,
                      float        rate = 0.5f,
                      const Array *p_bedrock = nullptr,
                      Array       *p_deposition_map = nullptr);

void thermal_conserve(Array       &z,
                      const Array *p_mask,
                      const Array &talus,
                      int          iterations = 10,
                      float        rate = 0.5f,
                      const Array *p_bedrock = nullptr,
                      Array       *p_deposition_map = nullptr); ///< @overload

void thermal_conserve(Array       &z,
                      float        talus,
                      int          iterations = 10,
                      float        rate = 0.5f,
                      const Array *p_bedrock = nullptr,
                      Array       *p_deposition_map = nullptr); ///< @overload

/**
 * @brief Apply iterative thermal flattening erosion on a heightmap.
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
 * @brief Apply thermal weathering erosion (Olsen model).
 *
 * @param z          Input array.
 * @param talus      Talus limit.
 * @param iterations Number of iterations.
 */
void thermal_olsen(Array &z, const Array &talus, int iterations);

void thermal_olsen(Array       &z,
                   const Array *p_mask,
                   const Array &talus,
                   int          iterations); ///< @overload

/**
 * @brief Apply thermal weathering erosion to give a scree like effect.
 *
 * @param z          Input array.
 * @param talus      Talus limit.
 * @param iterations Number of iterations.
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
 * @brief Apply thermal erosion using a 'rib' algorithm.
 *
 * @param z          Input heightmap.
 * @param iterations Number of iterations.
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
 * @param z                Input array.
 * @param talus            Talus limit.
 * @param iterations       Number of iterations.
 * @param p_deposition_map Optional deposition map output.
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
 * @brief Applies the thermal erosion process with a uniform slope threshold
 * (Schott model).
 *
 * @param z                Reference to elevation array.
 * @param talus            Constant threshold slope value or talus array.
 * @param iterations       Number of erosion iterations to apply.
 * @param intensity        Intensity factor controlling amount of change.
 * @param p_deposition_map Optional deposition map output.
 *
 * **Example**
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
 * @param[out] z                The heightmap modified in-place.
 * @param[in]  talus            Threshold slope angles.
 * @param[in]  zmax             Maximum allowed elevation.
 * @param[in]  iterations       Number of iterations.
 * @param[out] p_deposition_map Optional deposition map output.
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

} // namespace hmap::gpu
