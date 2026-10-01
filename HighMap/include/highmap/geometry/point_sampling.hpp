/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file point_sampling.hpp
 * @copyright Copyright (c) 2025
 */
#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "point_sampler.hpp"
#include "point_sampler/point.hpp"

#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/geometry/inverse_sampler_2d.hpp"
#include "highmap/geometry/point.hpp"

namespace hmap
{

/**
 * @brief Enumeration of point sampling methods.
 */
enum PointSamplingMethod : int
{
  RND_RANDOM,     ///< Purely random uniform sampling.
  RND_HALTON,     ///< Low-discrepancy Halton sequence sampling.
  RND_HAMMERSLEY, ///< Low-discrepancy Hammersley sequence sampling.
  RND_LHS,        ///< Latin Hypercube Sampling.
};

/**
 * @brief String-to-enum mapping for PointSamplingMethod.
 *
 * This map allows converting human-readable sampling method names into their
 * corresponding PointSamplingMethod values.
 */
static std::map<std::string, int> point_sampling_method_as_string = {
    {"Random", PointSamplingMethod::RND_RANDOM},
    {"Halton sequence", PointSamplingMethod::RND_HALTON},
    {"Hammersley sequence", PointSamplingMethod::RND_HAMMERSLEY},
    {"Latin Hypercube Sampling", PointSamplingMethod::RND_LHS},
};

/**
 * @brief Converts a 2D bounding box into coordinate ranges.
 *
 * This function takes a bounding box defined by four values (xmin, xmax, ymin,
 * ymax) stored in a glm::vec4 and returns an array of two ranges:
 *   - ranges[0] = {xmin, xmax}
 *   - ranges[1] = {ymin, ymax}
 *
 * @param  bbox Bounding box in the format {xmin, xmax, ymin, ymax}.
 * @return      An array of two pairs:
 *         - First pair: x-axis range (min, max)
 *         - Second pair: y-axis range (min, max)
 */
std::array<std::pair<float, float>, 2> bbox_to_ranges2d(const glm::vec4 &bbox);

/**
 * @brief Expand grid by translating and copying the values of the current
 * bounding box to the 8 first neighboring bounding boxes.
 * @param x     `x` coordinates.
 * @param y     `y` coordinates.
 * @param value values.
 * @param bbox  Bounding box.
 */
void expand_points_domain(std::vector<float> &x,
                          std::vector<float> &y,
                          std::vector<float> &value,
                          glm::vec4           bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Expand the grid by adding points on the boundaries of the bounding
 * box.
 * @param x            `x` coordinates.
 * @param y            `y` coordinates.
 * @param value        values.
 * @param bbox         Bounding box.
 * @param corner_value Value at the boundary points.
 */
void expand_points_at_domain_boundaries(std::vector<float> &x,
                                        std::vector<float> &y,
                                        std::vector<float> &value,
                                        glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f},
                                        float     boundary_value = 0.f);

/**
 * @brief Expand the grid by adding four points at the corner of the bounding
 * box.
 * @param x            `x` coordinates.
 * @param y            `y` coordinates.
 * @param value        values.
 * @param bbox         Bounding box.
 * @param corner_value Value at the four corner points.
 */
void expand_points_domain_corners(std::vector<float> &x,
                                  std::vector<float> &y,
                                  std::vector<float> &value,
                                  glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f},
                                  float     corner_value = 0.f);

/**
 * @brief Create a continuous 2D function from a sampled array.
 *
 * This adaptor takes a discrete 2D array (indexed by i, j) defined over a
 * bounding box and returns a continuous function \f$f(x, y)\f$ that
 * interpolates values at arbitrary 2D positions inside the bounding box.
 *
 * The function maps a point \f$p = (x, y)\f$ in world coordinates to array
 * indices \f$(i, j)\f$ and interpolates the array values.
 *
 * @param  array Input array containing sampled values on a regular grid.
 * @param  bbox  Bounding box of the array in world coordinates, given as
 *             `(xmin, xmax, ymin, ymax)`.
 * @return       A function `f(point)` that returns the interpolated value at a
 *               2D position.
 *
 * @par Example
 * @code {.cpp}
 * Array heightmap = [...];
 * glm::vec4 bbox = {0.0f, 100.0f, 0.0f, 100.0f};
 * auto f = make_pointwise_function_from_array(heightmap, bbox);
 *
 * ps::Point<float, 2> p{25.5f, 42.3f};
 * float value = f(p); // interpolated value at (25.5, 42.3)
 * @endcode
 */
std::function<float(const ps::Point<float, 2> &)>
make_pointwise_function_from_array(const Array &array, const glm::vec4 &bbox);

/**
 * @brief Generates random 2D points within a bounding box using a sampling
 * method.
 *
 * @param  count  Number of points to generate.
 * @param  seed   Random number generator seed.
 * @param  method Sampling method to use. Defaults to
 *                PointSamplingMethod::RND_RANDOM.
 * @param  bbox   Bounding box in which to generate the points (a,b,c,d = xmin,
 *                xmax, ymin, ymax). Defaults to the unit square {0.f, 1.f, 0.f,
 *                1.f}.
 * @return        A pair of float vectors {x_coords, y_coords}.
 */
std::array<std::vector<float>, 2> random_points(
    size_t                     count,
    std::uint32_t              seed,
    const PointSamplingMethod &method = PointSamplingMethod::RND_RANDOM,
    const glm::vec4           &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates random 2D points within a bounding box based on a density
 * map.
 *
 * Points are sampled with probability proportional to the density values.
 *
 * @param  count   Number of points to generate.
 * @param  density 2D array representing spatial density values.
 * @param  seed    Random number generator seed.
 * @param  bbox    Bounding box in which to generate the points (a,b,c,d = xmin,
 *                 xmax, ymin, ymax). Defaults to the unit square {0.f, 1.f,
 *                 0.f, 1.f}.
 * @return         A pair of float vectors {x_coords, y_coords}.
 */
std::array<std::vector<float>, 2> random_points_density(
    size_t           count,
    const Array     &density,
    std::uint32_t    seed,
    const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates random 2D points within a bounding box using inverse
 * transform sampling on a 2D density array.
 *
 * Points are sampled according to continuous row-marginal and conditional
 * column cumulative distributions interpolated from the 2D array.
 *
 * @param  count Number of points to generate.
 * @param  array 2D array representing spatial density / probability field.
 * @param  seed  Random number generator seed.
 * @param  bbox  Bounding box in which to generate the points (xmin, xmax, ymin,
 *               ymax). Defaults to the unit square {0.f, 1.f, 0.f, 1.f}.
 * @return       A pair of float vectors {x_coords, y_coords}.
 */
std::array<std::vector<float>, 2> random_points_inverse_sampling(
    size_t           count,
    const Array     &array,
    std::uint32_t    seed,
    const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates random 2D points with a minimum separation distance.
 *
 * Produces a blue-noise-like distribution where no two points are closer than
 * min_dist.
 *
 * @param  min_dist Minimum allowed distance between points.
 * @param  seed     Random number generator seed.
 * @param  bbox     Bounding box in which to generate the points (a,b,c,d =
 *                  xmin, xmax, ymin, ymax). Defaults to the unit square {0.f,
 * 1.f, 0.f, 1.f}.
 * @return          A pair of float vectors {x_coords, y_coords}.
 */
std::array<std::vector<float>, 2> random_points_distance(
    float            min_dist,
    std::uint32_t    seed,
    const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates random 2D points with distance constraints and a density
 * map.
 *
 * Points are separated by distances between @p min_dist and @p max_dist, and
 * are sampled proportionally to the given density map.
 *
 * @param  min_dist Minimum allowed distance between points.
 * @param  max_dist Maximum allowed distance between points.
 * @param  density  2D array representing spatial density values.
 * @param  seed     Random number generator seed.
 * @param  bbox     Bounding box in which to generate the points (a,b,c,d =
 *                  xmin, xmax, ymin, ymax). Defaults to the unit square {0.f,
 * 1.f, 0.f, 1.f}.
 * @return          A pair of float vectors {x_coords, y_coords}.
 */
std::array<std::vector<float>, 2> random_points_distance(
    float            min_dist,
    float            max_dist,
    const Array     &density,
    std::uint32_t    seed,
    const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates random 2D points with distances drawn from a power-law
 * distribution.
 *
 * Distances between points follow a power-law distribution between
 * @p dist_min and @p dist_max, with exponent @p alpha.
 *
 * @param  dist_min Minimum possible distance between points.
 * @param  dist_max Maximum possible distance between points.
 * @param  alpha    Power-law exponent (larger alpha favors shorter distances).
 * @param  seed     Random number generator seed.
 * @param  bbox     Bounding box in which to generate the points (a,b,c,d =
 *                  xmin, xmax, ymin, ymax). Defaults to the unit square {0.f,
 * 1.f, 0.f, 1.f}.
 * @return          A pair of float vectors {x_coords, y_coords}.
 */
std::array<std::vector<float>, 2> random_points_distance_power_law(
    float            dist_min,
    float            dist_max,
    float            alpha,
    std::uint32_t    seed,
    const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates random 2D points with distances drawn from a Weibull
 * distribution.
 *
 * Distances between points follow a Weibull distribution parameterized by
 * @p lambda (scale) and @p k (shape).
 *
 * @param  dist_min Minimum possible distance between points.
 * @param  lambda   Weibull distribution scale parameter.
 * @param  k        Weibull distribution shape parameter.
 * @param  seed     Random number generator seed.
 * @param  bbox     Bounding box in which to generate the points (a,b,c,d =
 *                  xmin, xmax, ymin, ymax). Defaults to the unit square {0.f,
 * 1.f, 0.f, 1.f}.
 * @return          A pair of float vectors {x_coords, y_coords}.
 */
std::array<std::vector<float>, 2> random_points_distance_weibull(
    float            dist_min,
    float            lambda,
    float            k,
    std::uint32_t    seed,
    const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates jittered grid-based 2D points.
 *
 * Points are placed on a regular grid, offset by jitter_amount, and optionally
 * staggered according to stagger_ratio.
 *
 * @param  count         Number of points to generate.
 * @param  jitter_amount Maximum jitter to apply in each axis (x, y).
 * @param  stagger_ratio Ratio of staggering between consecutive rows or
 *                       columns.
 * @param  seed          Random number generator seed.
 * @param  bbox          Bounding box in which to generate the points (a,b,c,d =
 * xmin, xmax, ymin, ymax). Defaults to the unit square {0.f, 1.f, 0.f, 1.f}.
 * @return               A pair of float vectors {x_coords, y_coords}.
 */
std::array<std::vector<float>, 2> random_points_jittered(
    size_t           count,
    const glm::vec2 &jitter_amount,
    const glm::vec2 &stagger_ratio,
    std::uint32_t    seed,
    const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Remove grid points that are outside a given bounding box.
 * @param x     `x` coordinates.
 * @param y     `y` coordinates.
 * @param value Values at points.
 * @param bbox  Bounding box.
 */
void remove_points_outside_bbox(std::vector<float> &x,
                                std::vector<float> &y,
                                std::vector<float> &value,
                                glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Rescale coordinate (x, y) so that they fit in a unit-square box based
 * on a given initial bounding box.
 * @param x[in, out] `x` coordinates (output).
 * @param y[in, out] `y` coordinates (output).
 * @param bbox       Initial bounding box.
 */
void rescale_points_to_unit_square(std::vector<float> &x,
                                   std::vector<float> &y,
                                   glm::vec4           bbox);

} // namespace hmap