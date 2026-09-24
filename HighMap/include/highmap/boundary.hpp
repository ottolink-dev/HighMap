/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file boundary.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for boundary condition functions and utilities.
 *
 * This header file contains declarations for functions and utilities that
 * manage and enforce boundary conditions in heightmap arrays. The functions
 * include methods for filling, extrapolating, and smoothing values at the
 * boundaries, as well as making arrays periodic and handling symmetry.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once
#include "highmap/array.hpp"
#include "highmap/math/distance_functions.hpp"
#include "highmap/math/profiles.hpp"

namespace hmap
{

/**
 * @enum PeriodicityType
 * @brief Describes the periodicity mode applied along map axes.
 */
// clang-format off
enum PeriodicityType : int
{
	PERIODICITY_X, ///< Periodic only along the X axis.
	PERIODICITY_Y, ///< Periodic only along the Y axis.
	PERIODICITY_XY ///< Periodic along both the X and Y axes.
};
// clang-format on

/**
 * @enum DomainBoundary
 * @brief Describes which domain boundary.
 */
// clang-format off
enum DomainBoundary : int
{
	BOUNDARY_LEFT, ///< i = 0
	BOUNDARY_RIGHT, ///< i = nx - 1
	BOUNDARY_TOP, ///< j = ny - 1
	BOUNDARY_BOTTOM ///< j = 0
};
// clang-format on

/**
 * @brief Performs linear extrapolation of values at the borders of an array
 * (e.g., `i = 0`, `j = 0`, etc.) based on the inner values of the array.
 *
 * This function modifies the borders of the input array by applying a linear
 * extrapolation method, which uses the values inside the array to estimate and
 * fill in the border values. The extrapolation is influenced by a relaxation
 * coefficient (`sigma`) and a buffer depth (`nbuffer`), which determines how
 * many layers of the border are extrapolated.
 *
 * @param array   Reference to the input array whose borders need extrapolation.
 * @param nbuffer Optional parameter specifying the buffer depth, i.e., the
 *                number of layers at the border to extrapolate. Default is 1.
 * @param sigma   Optional relaxation coefficient that adjusts the influence of
 *                inner values on the extrapolated border values. Default is
 *                0.0.
 *
 * @see           fill_borders()
 */
void extrapolate_borders(Array &array, int nbuffer = 1, float sigma = 0.f);

/**
 * @brief Applies a falloff effect to the input array based on distance.
 *
 * This function modifies the elements of the input `array` by applying a
 * falloff effect, which decreases the values based on their distance from a
 * central point or area. The strength of the falloff, the type of distance
 * function used, and optional noise can be specified.
 *
 * @param array    Reference to the array that will be modified by the falloff
 *                 effect.
 * @param strength The strength of the falloff effect. A higher value results in
 *                 a stronger falloff. Default is 1.0f.
 * @param dist_fct The distance function to be used for calculating the falloff.
 *                 Options include Euclidian and others, with
 * `DistanceFunction::EUCLIDIAN` as the default.
 * @param p_noise  Optional pointer to an array that provides noise to be added
 *                 to the falloff effect. If nullptr (default), no noise is
 *                 added.
 * @param bbox     A 4D vector representing the bounding box within which the
 *                 falloff effect is applied. The default is {0.f, 1.f, 0.f,
 *                 1.f}.
 *
 * **Example**
 * @include ex_falloff.cpp
 *
 * **Result**
 * @image html ex_falloff.png
 */
void falloff(Array           &array,
             float            strength = 1.f,
             DistanceFunction dist_fct = DistanceFunction::EUCLIDIAN,
             const Array     *p_noise = nullptr,
             glm::vec4        bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Fills the border values of an array (e.g., `i = 0`, `j = 0`, etc.)
 * based on the values of the first neighboring cells.
 *
 * This function modifies the border values of the input array by copying values
 * from their immediate neighbors. The operation ensures that border values are
 * consistent with their adjacent cells, typically used to prepare the array for
 * further processing.
 *
 * @param array Reference to the input array whose borders need to be filled.
 *
 * @see         extrapolate_borders()
 */
void fill_borders(Array &array);

void fill_borders(Array &array, int nbuffer); ///< @overload

/**
 * @brief Creates and returns a new array with additional buffer zones at the
 * boundaries, where the buffer values are filled either by symmetry or by
 * zero-padding.
 *
 * This function takes an input array and generates a new array with extra
 * layers (buffers) added to its boundaries. The size of these buffer zones is
 * specified by the `buffers` parameter, and the values within these buffers can
 * be filled either by reflecting the values at the boundaries (symmetry) or by
 * padding with zeros.
 *
 * @param  array        Reference to the input array.
 * @param  buffers      A vector specifying the buffer sizes for the east, west,
 *                      south, and north boundaries.
 * @param  zero_padding Optional boolean flag to use zero-padding instead of
 *                      symmetry for filling the buffer values. Default is
 *                      `false`.
 * @return              Array A new array with buffers added at the boundaries.
 */
Array generate_buffered_array(const Array &array,
                              glm::ivec4   buffers,
                              bool         zero_padding = false);

/**
 * @brief Adjusts the input array to be periodic in both directions by
 * transitioning smoothly at the boundaries.
 *
 * This function modifies the input array so that the values at the boundaries
 * transition smoothly, making the array periodic in both the horizontal and
 * vertical directions. The width of the transition zone at the boundaries is
 * controlled by the `nbuffer` parameter.
 *
 * @param array   Reference to the input array to be made periodic.
 * @param nbuffer The width of the transition zone at the boundaries.
 *
 * **Example**
 * @include ex_make_periodic.cpp
 *
 * **Result**
 * @image html ex_make_periodic0.png
 * @image html ex_make_periodic1.png
 */
void make_periodic(
    Array                 &array,
    int                    nbuffer,
    const PeriodicityType &periodicity_type = PeriodicityType::PERIODICITY_XY);

/**
 * @brief Creates a periodic array in both directions using a stitching
 * operation that minimizes errors at the boundaries.
 *
 * This function generates a new array that is periodic in both directions by
 * applying a stitching operation at the boundaries. The stitching process aims
 * to minimize discrepancies, creating a seamless transition between the edges
 * of the array. The `overlap` parameter determines the extent of the overlap
 * during the stitching, based on the half-size of the domain. If
 * `overlap` is set to 1, the transition spans the entire domain.
 *
 * @param  array   Reference to the input array to be made periodic.
 * @param  overlap A float value representing the overlap based on the domain's
 *                 half-size. An overlap of 1 means the transition spans the
 *                 whole domain.
 * @return         Array A new array that is periodic in both directions with
 *                 minimized boundary errors.
 *
 * **Example**
 * @include ex_make_periodic_stitching.cpp
 *
 * **Result**
 * @image html ex_make_periodic_stitching0.png
 * @image html ex_make_periodic_stitching1.png
 */
Array make_periodic_stitching(const Array &array, float overlap);

/**
 * @brief Creates a tiled, periodic array by applying a transition with overlap
 * in both directions.
 *
 * This function generates a new array that is periodic and tiled according to
 * the specified tiling dimensions. The periodicity is achieved by applying a
 * transition at the boundaries, where the extent of overlap is determined by
 * the `overlap` parameter. The `tiling` parameter specifies the number of tiles
 * in the horizontal and vertical directions.
 *
 * @param  array   Reference to the input array to be tiled and made periodic.
 * @param  overlap A float value representing the overlap based on the domain's
 *                 half-size. If `overlap` is 1, the transition spans the entire
 *                 domain on both sides.
 * @param  tiling  A 2D vector specifying the number of tiles in the horizontal
 *                 and vertical directions.
 * @return         Array A new tiled array that is periodic in both directions.
 *
 * **Example**
 * @include ex_make_periodic_tiling.cpp
 *
 * **Result**
 * @image html ex_make_periodic_tiling.png
 */
Array make_periodic_tiling(const Array &array,
                           float        overlap,
                           glm::ivec2   tiling);

/**
 * @brief Picks a boundary cell using weighted random selection.
 *
 * Selects a cell on the specified boundary, where each candidate is assigned a
 * weight. Higher weights increase the probability of selection.
 *
 * Optional weighting factors:
 * - favor_boundary_center: favors cells near the boundary center.
 * - favor_lower_elevation: favors lower values in @p z (normalized).
 * - favor_sinks: favors local minima (cells <= neighbors).
 *
 * Weights are combined multiplicatively. Selection is deterministic given @p
 * seed.
 *
 * @param  z                     Input 2D scalar field.
 * @param  boundary              Boundary to sample from.
 * @param  seed                  Random seed.
 * @param  favor_boundary_center Enable center bias.
 * @param  favor_lower_elevation Enable low-value bias.
 * @param  favor_sinks           Enable sink bias.
 *
 * @return                       Selected boundary cell coordinates.
 *
 * **Example**
 * @include ex_pick_boundary_cell.cpp
 *
 * **Result**
 * @image html ex_pick_boundary_cell.png
 */
glm::ivec2 pick_boundary_cell(const Array   &z,
                              DomainBoundary boundary,
                              std::uint32_t  seed,
                              bool           favor_boundary_center = true,
                              bool           favor_lower_elevation = true,
                              bool           favor_sinks = true);

/**
 * @brief Enforces specific values at the boundaries of the array.
 *
 * This function sets the values at the borders of the input array to specified
 * values for each side (east, west, south, and north). The size of the buffer
 * zone at each border can also be defined, allowing precise control over how
 * much of the boundary is modified.
 *
 * @param array         Reference to the input array whose borders are to be
 *                      modified.
 * @param border_values A vector specifying the values to set at the east, west,
 *                      south, and north borders.
 * @param buffer_sizes  A vector specifying the size of the buffer zones at the
 *                      east, west, south, and north borders.
 *
 * **Example**
 * @include ex_set_borders.cpp
 *
 * **Result**
 * @image html ex_set_borders.png
 */
void set_borders(Array     &array,
                 glm::vec4  border_values,
                 glm::ivec4 buffer_sizes);

/**
 * @brief Enforces specific values at the domain boundaries based on distance.
 *
 * This overloaded function sets values near the boundaries of the unit square
 * domain [0, 1] x [0, 1] based on continuous distances. The array corresponds
 * to the specified bounding box `bbox`.
 *
 * @param array         Reference to the input array whose borders are to be
 *                      modified.
 * @param border_values A vector specifying the values to set at the west, east,
 *                      south, and north borders.
 * @param buffer_sizes  A vector specifying the buffer distance at the west,
 *                      east, south, and north borders.
 * @param bbox          Domain bounding box {xmin, xmax, ymin, ymax}. Default is
 *                      {0.f, 1.f, 0.f, 1.f}.
 */
void set_borders(Array    &array,
                 glm::vec4 border_values,
                 glm::vec4 buffer_sizes,
                 glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Enforces a uniform value at all boundaries of the array.
 *
 * This overloaded function sets the same value at all borders of the input
 * array. The size of the buffer zone at each border can be defined with a
 * single value that applies uniformly to all sides.
 *
 * @param array         Reference to the input array whose borders are to be
 *                      modified.
 * @param border_values The value to set at all borders.
 * @param buffer_sizes  The size of the buffer zone to apply uniformly at all
 *                      borders.
 */
void set_borders(Array &array, float border_values, int buffer_sizes);

/**
 * @brief Enforces a uniform value at all domain boundaries based on distance.
 *
 * This overloaded function sets the same value at all borders of the unit
 * square domain [0, 1] x [0, 1] based on a continuous buffer distance. The
 * array corresponds to the specified bounding box `bbox`.
 *
 * @param array         Reference to the input array whose borders are to be
 *                      modified.
 * @param border_values The value to set at all borders.
 * @param buffer_sizes  The buffer distance to apply uniformly at all borders.
 * @param bbox          Domain bounding box {xmin, xmax, ymin, ymax}. Default is
 *                      {0.f, 1.f, 0.f, 1.f}.
 */
void set_borders(Array    &array,
                 float     border_values,
                 float     buffer_sizes,
                 glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Fills the values at the domain borders using symmetry over a specified
 * buffer depth.
 *
 * This function fills the borders of the input array by reflecting the values
 * inside the array symmetrically. The depth of the buffer at each border (east,
 * west, south, north) is specified by the `buffer_sizes` parameter.
 *
 * @param array        Reference to the input array whose borders are to be
 *                     filled using symmetry.
 * @param buffer_sizes A vector specifying the buffer sizes at the east, west,
 *                     south, and north borders.
 */
void sym_borders(Array &array, glm::ivec4 buffer_sizes);

/**
 * @brief Fills the border values (e.g., `i = 0`, `j = 0`, etc.) of the array
 * with zeros.
 *
 * This function sets all the values at the borders of the input array to zero.
 * It can be used to zero out boundary values, preparing the array for specific
 * computations.
 *
 * @param array Reference to the input array whose borders are to be zeroed.
 *
 * @see         fill_borders()
 */
void zeroed_borders(Array &array);

/**
 * @brief Smoothly attenuates array values toward zero near the domain
 * boundaries.
 *
 * Applies a radial profile centered at `center` with the specified `radius`.
 * The effect can be shaped using `radial_profile`, `profile_param`, and an
 * optional noise modulation array.
 *
 * @param array          Input/output array.
 * @param radial_profile Radial attenuation profile.
 * @param profile_param  Profile-specific parameter.
 * @param amount         Effect strength in [0, 1].
 * @param distance       Distance function.
 * @param dfa            Distance axes used for distance evaluation.
 * @param center         Profile center in normalized coordinates.
 * @param radius         Profile radius.
 * @param p_noise_r      Optional modulation array.
 * @param bbox           Domain bounding box.
 *
 * @include ex_zeroed_edges.cpp
 * @image html ex_zeroed_edges.png
 */
void zeroed_edges(Array        &array,
                  RadialProfile radial_profile = RadialProfile::RP_SMOOTHSTEP,
                  float         profile_param = 0.f,
                  float         amount = 1.f,
                  DistanceFunction     distance = DistanceFunction::EUCLIDIAN,
                  DistanceFunctionAxis dfa = DistanceFunctionAxis::DFA_XY,
                  glm::vec2            center = {0.5f, 0.5f},
                  float                radius = 0.5f,
                  const Array         *p_noise_r = nullptr,
                  glm::vec4            bbox = {0.f, 1.f, 0.f, 1.f});

} // namespace hmap
