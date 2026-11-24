/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file curvature.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file defining a collection of functions for curvature analysis.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @brief Computes a signed level-set curvature.
 *
 * @param  array        Input array.
 * @param  prefilter_ir Optional smoothing radius for curvature evaluation.
 * @return              Signed curvature.
 *
 * See unit tests: @ref test_level_set_curvature.cpp
 */
Array level_set_curvature(const Array &array, int prefilter_ir);

/**
 * \brief Types of curvature measures.
 *
 * Defines the different curvature descriptors that can be computed on a scalar
 * field (e.g. heightmap). They include classical differential curvatures as
 * well as terrain-analysis specific measures.
 */
// clang-format off
enum CurvatureType : int
{
	CT_MIN,          //!< Minimum curvature (principal)
	CT_MAX,          //!< Maximum curvature (principal)
	CT_MEAN,         //!< Mean curvature
	CT_GAUSSIAN,     //!< Gaussian curvature
	CT_PROFILE,      //!< Profile curvature (along slope direction)
	CT_CONTOUR,      //!< Contour curvature (across slope)
	CT_TANGENTIAL,   //!< Tangential curvature
	CT_ACCUMULATION, //!< Accumulation curvature (H² − K²)
	CT_SHAPE_INDEX,  //!< Shape index (normalized curvature classification)
	CT_UNSPHERICITY, //!< Unsphericity (√(H² − K))
	CT_RING,         //!< Ring curvature (planform torsion-like measure)
	CT_ROTOR,        //!< Rotor curvature (rotational component)
};
// clang-format on

/**
 * \brief Compute curvature measures on a scalar field.
 *
 * Computes various curvature descriptors from the input array \p z. If \p ir ==
 * 0, local derivatives (3×3 stencil) are used. Otherwise, a quadric surface is
 * fitted over a (2·ir+1) window.
 *
 * \param z Input scalar field (e.g. heightmap)
 * \param ir Neighborhood radius (0 for local, >0 for quadric fit)
 * \param curvature_type Type of curvature to compute.
 *
 * \return Array of same size containing the selected curvature
 *
 * **Example**
 * @include ex_curvature_quadric.cpp
 *
 * **Result**
 * @image html ex_curvature_quadric.png
 */
Array curvature_quadric(const Array &z, int ir, CurvatureType curvature_type);

} // namespace hmap

namespace hmap::gpu
{

/*! @brief See hmap::level_set_curvature */
Array level_set_curvature(const Array &array, int prefilter_ir);

/*! @brief See hmap::curvature_quadric */
Array curvature_quadric(const Array  &z,
                        int           ir,
                        CurvatureType curvature_type,
                        bool          approx_algo = false);

} // namespace hmap::gpu