/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file cloud.hpp
 * @brief Definition of the `Cloud` class for manipulating sets of 2D points.
 *
 * This file contains the definition of the `Cloud` class, which is used to
 * manage and manipulate unordered sets of points in a 2D space. The `Cloud`
 * class provides functionalities for creating clouds of points, including
 * random generation, value assignment, and data export. It also supports
 * operations such as computing bounding boxes, distances, and interpolating
 * values from arrays.
 *
 * The class includes methods for adding and removing points, randomizing
 * positions and values, remapping values, and exporting data to CSV files.
 * Additional functionalities include calculating Signed Distance Functions
 * (SDF) and projecting point clouds onto arrays.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once
#include <cmath>
#include <optional>

#include "highmap/array.hpp"
#include "highmap/geometry/point.hpp"
#include "highmap/geometry/point_sampling.hpp"
#include "highmap/interpolate/interpolate2d.hpp"

namespace hmap
{

class Graph;

/**
 * @class Cloud
 * @brief Represents a collection of unordered points in 2D space.
 *
 * The `Cloud` class provides functionality to manage and manipulate an
 * unordered set of points in a 2D space. It supports various operations such as
 * adding points, calculating the centroid, merging with other point clouds, and
 * more. This class is useful for applications involving point cloud processing,
 * geometric computations, and spatial analysis.
 *
 * See unit tests: @ref test_cloud.cpp
 */
class Cloud
{
public:
  std::vector<Point> points = {}; ///< Points of the cloud.

  // ==========================================================================
  //  Constructors
  // ==========================================================================

  /**
   * @brief Default constructor for the Cloud class.
   *
   * Initializes an empty cloud with no points.
   */
  Cloud(){};

  virtual ~Cloud() = default;

  /**
   * @brief Constructs a new Cloud object with random positions and values.
   *
   * @param npoints Number of points to generate.
   * @param seed    Random seed used to generate the points. Using the same seed
   *                will produce the same set of points.
   * @param bbox    Bounding box within which the points will be generated. The
   *                bounding box is defined as {xmin, xmax, ymin, ymax}.
   */
  Cloud(int npoints, std::uint32_t seed, glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f});

  /**
   * @brief Constructs a new Cloud object based on a list of existing points.
   *
   * @param points A vector of `Point` objects representing the cloud's points.
   */
  Cloud(const std::vector<Point> &points) : points(points){};

  /**
   * @brief Constructs a new Cloud object from lists of `x` and `y` coordinates.
   *
   * @param x             A vector of `x` coordinates for the points.
   * @param y             A vector of `y` coordinates for the points.
   * @param default_value The default value assigned to each point. Defaults to
   * 0 if not specified.
   */
  Cloud(const std::vector<float> &x,
        const std::vector<float> &y,
        float                     default_value = 0.f);

  /**
   * @brief Constructs a new Cloud object from lists of `x` and `y` coordinates
   * with assigned values.
   *
   * @param x A vector of `x` coordinates for the points.
   * @param y A vector of `y` coordinates for the points.
   * @param v A vector of values associated with each point.
   */
  Cloud(const std::vector<float> &x,
        const std::vector<float> &y,
        const std::vector<float> &v);

  /**
   * @brief Construct a point cloud from grid indices mapped to a bounding box.
   *
   * Each index is normalized by the grid shape and remapped to the given
   * bounding box. The resulting points are stored as 3D positions (z = 1).
   *
   * @param indices Input grid indices.
   * @param shape   Grid dimensions.
   * @param bbox    Bounding box (xmin, xmax, ymin, ymax).
   */
  Cloud(const std::vector<glm::ivec2> &indices,
        const glm::ivec2              &shape,
        const glm::vec4               &bbox = {0.f, 1.f, 0.f, 1.f});

  /**
   * @brief  Constructs a new Cloud object from lists of xyz data as glm::vec3.
   * */
  Cloud(const std::vector<glm::vec3> &xyv);

  /**
   * @brief Add a new point to the cloud.
   *
   * @param p The point to be added to the cloud.
   */
  void add_point(const Point &p);

  /**
   * @brief Remove a point from the cloud.
   * @param point_idx Index of the point to be removed.
   */
  void remove_point(int point_idx);

  // ==========================================================================
  //  Accessors
  // ==========================================================================

  /**
   * @brief Get the bounding box of the cloud.
   *
   * @return glm::vec4 The bounding box of the cloud in the format `[xmin, xmax,
   *         ymin, ymax]`.
   */
  glm::vec4 get_bbox() const;

  /**
   * @brief Calculates the centroid of a set of points.
   *
   * This function computes the center (or centroid) of a collection of points
   * by averaging their coordinates. It sums up all the point coordinates and
   * then divides the result by the total number of points to obtain the average
   * position. The centroid represents the geometric center of the point cloud.
   *
   * @return Point The computed center, represented as a `Point` object, which
   *         contains the average (x, y) coordinates of the points.
   */
  Point get_center() const;

  /**
   * @brief Computes the indices of the points that form the convex hull of a
   * set of points.
   *
   * @return std::vector<int> A vector containing the indices of the points that
   *         make up the convex hull, listed in order.
   *
   * **Example**
   * @include ex_cloud_get_convex_hull.cpp
   *
   * **Result**
   * @image html ex_cloud_get_convex_hull.png
   */
  std::vector<int> get_convex_hull() const;

  /**
   * @brief Get the values assigned to the points in the cloud.
   * @return std::vector<float> A vector containing the values of all points in
   *         the cloud.
   */
  std::vector<float> get_values() const;

  /**
   * @brief Get the maximum value among the points in the cloud.
   * @return float The maximum value among the points.
   */
  float get_values_max() const;

  /**
   * @brief Get the minimum value among the points in the cloud.
   * @return float The minimum value among the points.
   */
  float get_values_min() const;

  /**
   * @brief Get the `x` coordinates of the points in the cloud.
   * @return std::vector<float> A vector containing the `x` coordinates of the
   *         points.
   */
  std::vector<float> get_x() const;

  /**
   * @brief Get the concatenated `x` and `y` coordinates of the points in the
   * cloud.
   *
   * This method returns a vector containing the `x` and `y` coordinates of all
   * points in the cloud, arranged in the form `[x0, y0, x1, y1, ...]`.
   *
   * @return std::vector<float> A vector containing the concatenated `x` and `y`
   * coordinates of the points.
   */
  std::vector<float> get_xy() const;

  /**
   * @brief Get the `y` coordinates of the points in the cloud.
   * @return std::vector<float> A vector containing the `y` coordinates of the
   *         points.
   */
  std::vector<float> get_y() const;

  /**
   * \brief Find the index of the nearest point in the cloud.
   *
   * Computes the point whose (x,y) position is closest to \p xy using squared
   * Euclidean distance.
   *
   * \note This implementation is not optimized and is intended for occasional
   * single queries. It performs a linear search and is not suitable for
   * large-scale or repeated nearest-neighbor queries.
   *
   * \param xy Query position
   * \return Index of the nearest point
   */
  size_t nearest_point(const glm::vec2 &xy) const;

  /**
   * @brief Set points of the using x, y coordinates.
   *
   * @param x A vector containing the `x` coordinates of the points.
   * @param y A vector containing the `y` coordinates of the points.
   */
  void set_points(const std::vector<float> &x, const std::vector<float> &y);

  /**
   * @brief Set new values for the cloud points.
   * @param new_values A vector of new values to assign to the points.
   */
  void set_values(const std::vector<float> &new_values);

  /**
   * @brief Set a single value for all cloud points.
   *
   * This method assigns the same value to all points in the cloud.
   *
   * @param new_value The value to assign to all points.
   */
  void set_values(float new_value);

  /**
   * @brief Set the values of the cloud points using values from an underlying
   * array.
   *
   * This method assigns values to the cloud points by interpolating values from
   * an input array. The positions of the cloud points are mapped to the array
   * using the specified bounding box.
   *
   * @param array The input array from which to derive the values.
   * @param bbox  The bounding box that defines the mapping from the cloud
   * points' coordinates to the array's coordinates.
   */
  void set_values_from_array(const Array     &array,
                             const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

  /**
   * @brief Sets point values based on their distance to the bounding box
   * border.
   *
   * For each point in the cloud, this method computes the shortest distance to
   * the edges of the given bounding box and stores it as the point's value.
   *
   * @param bbox Bounding box in the format {xmin, xmax, ymin, ymax}.
   */
  void set_values_from_border_distance(
      const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

  /**
   * @brief Set the values of the cloud points based on the distance to the
   * convex hull of the cloud.
   *
   * This method assigns values to the cloud points based on their distance to
   * the convex hull of the cloud. The convex hull must be initialized using
   * `to_graph_delaunay()` before this method is called to ensure the correct
   * indices of the convex hull points are available.
   */
  void set_values_from_chull_distance();

  /**
   * @brief Sets point values based on the distance to their nearest neighbor.
   */
  void set_values_from_min_distance();

  /**
   * @brief Check whether the cloud has no points.
   * @return true if the cloud contains no points, false otherwise.
   */
  bool empty() const;

  /**
   * @brief Get the number of points in the cloud.
   * @return size_t The number of points in the cloud.
   */
  size_t size() const;

  // ==========================================================================
  //  Basic Ops
  // ==========================================================================

  /**
   * @brief Clear all data from the cloud.
   */
  void clear();

  /**
   * @brief Print information about the cloud's points.
   */
  void print();

  /**
   * @brief Randomize the positions and values of the cloud points.
   *
   * @param seed Random seed number for generating positions and values.
   * @param bbox Bounding box within which the points will be randomized.
   */
  void randomize(std::uint32_t seed, glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f});

  /**
   * @brief Remap the values of the cloud points to a target range.
   *
   * @param vmin The lower bound of the target range.
   * @param vmax The upper bound of the target range.
   */
  void remap_values(float vmin, float vmax);

  /**
   * @brief Randomly perturbs the positions and values of all points in the
   * cloud.
   *
   * @param dx   Scale factor for the random displacement along the X-axis.
   * @param dy   Scale factor for the random displacement along the Y-axis.
   * @param seed Seed value for the pseudo-random number generator, ensuring
   *             reproducibility.
   * @param dv   Scale factor for the random displacement applied to the point's
   *             value component `v`.
   */
  void shuffle(float dx, float dy, std::uint32_t seed, float dv = 0.f);

  /**
   * @brief Snap points to the bounding box edges and corners.
   *
   * Points within a tolerance distance from the bounding box edges are
   * projected onto the closest edge. Afterwards, the closest point to each
   * corner is snapped exactly to that corner.
   *
   * @param bbox            Bounding box defined as (xmin, xmax, ymin, ymax).
   * @param tolerance_ratio Ratio used to compute the snapping distance relative
   *                        to the average point spacing.
   */
  void snap_points_to_bounding_box(const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f},
                                   float            tolerance_ratio = 1.f);

  // ==========================================================================
  //  Conversion / IO
  // ==========================================================================

  /**
   * @brief Loads point data from a CSV file into the Cloud object.
   *
   * This function reads a CSV file where each line contains either 2D (X, Y) or
   * 3D (X, Y, Z) point data. The function automatically detects the
   * dimensionality of the points based on the number of values per line. The
   * loaded points are stored in the `points` member of the Cloud object.
   *
   * @param  fname The path to the CSV file to be read.
   * @return       true if the file was successfully read and the points were
   *               loaded, false otherwise.
   *
   * @note The CSV file must be well-formed, with each line containing either 2
   * or 3 comma-separated values. Lines with an unexpected number of values will
   * cause the function to return false.
   *
   * @warning If the file cannot be opened or contains invalid data (e.g.,
   * non-numeric values), the function will log an error and return false.
   */
  bool from_csv(const std::string &fname);

  /**
   * @brief Project the cloud points onto an array.
   *
   * This method projects the cloud points' values onto a given array, mapping
   * their positions within the specified bounding box. The resulting array will
   * be populated with the values of the points at their corresponding
   * positions.
   *
   * @param array The input array where the cloud points' values will be
   *              projected.
   * @param bbox  The bounding box that defines the mapping from the cloud
   * points' coordinates to the array's coordinates.
   */
  void to_array(Array &array, glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f}) const;

  /*! @brief See hmap::to_array */
  Array to_array(glm::ivec2 shape, glm::vec4 bbox) const;

  /**
   * @brief Interpolate the values of an array using the cloud points.
   *
   * This method populates an array with interpolated values based on the
   * positions and values of the cloud points. The interpolation method can be
   * specified, and optional noise arrays can be used for domain warping.
   *
   * @param array                The output array that will be populated with
   *                             interpolated values.
   * @param bbox                 The bounding box that defines the cloud's
   *                             coordinate system.
   * @param interpolation_method The method used for interpolation (e.g.,
   * nearest neighbor, bilinear).
   * @param p_noise_x            Optional reference to a noise array applied to
   *                             the x-coordinates for domain warping (not in
   *                             pixels).
   * @param p_noise_y            Optional reference to a noise array applied to
   *                             the y-coordinates for domain warping (not in
   *                             pixels).
   * @param bbox_array           The bounding box of the destination array.
   *
   * **Example**
   * @include ex_cloud_to_array_interp.cpp
   *
   * **Result**
   * @image html ex_cloud_to_array_interp.png
   */
  void to_array_interp(Array                &array,
                       glm::vec4             bbox = {0.f, 1.f, 0.f, 1.f},
                       InterpolationMethod2D interpolation_method =
                           InterpolationMethod2D::ITP2D_DELAUNAY,
                       Array    *p_noise_x = nullptr,
                       Array    *p_noise_y = nullptr,
                       glm::vec4 bbox_array = {0.f, 1.f, 0.f, 1.f}) const;

  /**
   * @brief Export the cloud data to a CSV file.
   * @param fname The name of the output CSV file.
   */
  void to_csv(const std::string &fname) const;

  /**
   * @brief Convert the cloud to a graph using Delaunay triangulation.
   * @return Graph The resulting graph from Delaunay triangulation.
   */
  Graph to_graph_delaunay();

  /**
   * @brief Saves the current data as a PNG image file.
   *
   * @param fname The file name for the output PNG image. This should include
   * the file extension (e.g., "output.png").
   * @param cmap  An integer specifying the colormap to be used for rendering
   * the data. This index refers to a predefined colormap.
   * @param bbox  A `glm::vec4` specifying the bounding box of the data to be
   *              included in the image. It is given as {xmin, xmax, ymin,
   * ymax}. The default is {0.f, 1.f, 0.f, 1.f}.
   * @param depth An integer specifying the bit depth of the image. It should be
   * a value defined by OpenCV (e.g., `CV_8U` for 8-bit unsigned). The default
   * is `CV_8U`.
   * @param shape A `glm::ivec2` specifying the dimensions of the output image.
   * It is given as {width, height}. The default is {512, 512}.
   */
  void to_png(const std::string &fname,
              int                cmap,
              glm::vec4          bbox = {0.f, 1.f, 0.f, 1.f},
              int                depth = CV_8U,
              glm::ivec2         shape = {512, 512});

  /**
   * @brief Convert path points to a vector of 3D positions.
   * @return Vector of points as (x, y, v).
   */
  std::vector<glm::vec3> to_vec3() const;
};

// ==========================================================================
//  Functions
// ==========================================================================

/**
 * @brief Compute a distance field from a point cloud.
 *
 * @param  cloud      Input point cloud.
 * @param  shape      Output array dimensions.
 * @param  bbox_array Bounding box of the output domain.
 * @param  p_noise_x  Optional x-direction noise (domain warp).
 * @param  p_noise_y  Optional y-direction noise (domain warp).
 * @return            Array of distances to the nearest point.
 */
Array cloud_sdf_to_array(const Cloud &cloud,
                         glm::ivec2   shape,
                         glm::vec4    bbox_array = {0.f, 1.f, 0.f, 1.f},
                         const Array *p_noise_x = nullptr,
                         const Array *p_noise_y = nullptr);

/**
 * @brief Checks whether the point cloud contains duplicate points.
 *
 * @param  cloud   Input point cloud.
 * @param  eps     Distance tolerance for considering two points identical.
 * @param  xy_only If true, compares only the X and Y coordinates.
 * @return         True if duplicate points are found, false otherwise.
 */
bool has_duplicates(const Cloud &cloud, float eps = 1e-9f, bool xy_only = true);

/**
 * @brief Interpolate values from an array at the points' `(x, y)` locations.
 *
 * This method computes interpolated values for each point in the cloud based on
 * its `(x, y)` coordinates and an underlying array, using bilinear
 * interpolation within the specified bounding box.
 *
 * @param  array The input array from which to interpolate values.
 * @param  bbox  The bounding box of the array.
 * @return       std::vector<float> A vector containing the interpolated values
 *               for each point.
 */
std::vector<float> interpolate_values_from_array(const Cloud &cloud,
                                                 const Array &array,
                                                 glm::vec4    bbox);

/**
 * @brief Merges two point clouds into one.
 *
 * @param  cloud1 The first point cloud to be merged. This cloud will be the
 *                base cloud to which the points from the second cloud are
 *                added.
 * @param  cloud2 The second point cloud whose points will be appended to the
 *                first cloud.
 * @return        Cloud The resulting point cloud that includes all points from
 *                both `cloud1` and `cloud2`.
 */
Cloud merge_cloud(const Cloud &cloud1, const Cloud &cloud2);

/**
 * @brief Merges multiple point clouds into a single cloud.
 *
 * @param clouds A vector of Cloud objects to be merged.
 */
Cloud merge_clouds(const std::vector<Cloud> &clouds);

/**
 * @brief Generates a random cloud of points within a bounding box.
 *
 * @param  count  Number of points to generate.
 * @param  seed   Random number generator seed.
 * @param  method Sampling method to use. Defaults to
 *                PointSamplingMethod::RND_RANDOM.
 * @param  bbox   Bounding box in which to generate the points (a,b,c,d = xmin,
 *                xmax, ymin, ymax). Defaults to the unit square {0.f, 1.f, 0.f,
 *                1.f}.
 * @return        A Cloud containing the generated points.
 *
 * **Example**
 * @include ex_point_sampling.cpp
 *
 * **Result**
 * @image html ex_point_sampling0.png
 * @image html ex_point_sampling1.png
 * @image html ex_point_sampling2.png
 * @image html ex_point_sampling3.png
 */
Cloud random_cloud(
    size_t                     count,
    std::uint32_t              seed,
    const PointSamplingMethod &method = PointSamplingMethod::RND_LHS,
    const glm::vec4           &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a random cloud of points based on a spatial density map.
 *
 * The probability of placing a point at a location is proportional to the
 * density value at that location.
 *
 * @param  count   Number of points to generate.
 * @param  density 2D array representing spatial density values.
 * @param  seed    Random number generator seed.
 * @param  bbox    Bounding box in which to generate the points (a,b,c,d = xmin,
 *                 xmax, ymin, ymax). Defaults to the unit square {0.f, 1.f,
 *                 0.f, 1.f}.
 * @return         A Cloud containing the generated points.
 *
 * **Example**
 * @include ex_point_sampling.cpp
 *
 * **Result**
 * @image html ex_point_sampling0.png
 * @image html ex_point_sampling1.png
 * @image html ex_point_sampling2.png
 * @image html ex_point_sampling3.png
 */
Cloud random_cloud_density(size_t           count,
                           const Array     &density,
                           std::uint32_t    seed,
                           const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a random cloud of points separated by at least a given
 * minimum distance.
 *
 * Points are distributed randomly but maintain a minimum separation, producing
 * a blue-noise-like distribution.
 *
 * @param  min_dist Minimum allowed distance between any two points.
 * @param  seed     Random number generator seed.
 * @param  bbox     Bounding box in which to generate the points (a,b,c,d =
 *                  xmin, xmax, ymin, ymax). Defaults to the unit square {0.f,
 * 1.f, 0.f, 1.f}.
 * @return          A Cloud containing the generated points.
 *
 * **Example**
 * @include ex_point_sampling.cpp
 *
 * **Result**
 * @image html ex_point_sampling0.png
 * @image html ex_point_sampling1.png
 * @image html ex_point_sampling2.png
 * @image html ex_point_sampling3.png
 */
Cloud random_cloud_distance(float            min_dist,
                            std::uint32_t    seed,
                            const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a random cloud of points separated by a distance range,
 * influenced by a density map.
 *
 * Points maintain a separation between @p min_dist and @p max_dist, and are
 * distributed according to the provided density map.
 *
 * @param  min_dist Minimum allowed distance between points.
 * @param  max_dist Maximum allowed distance between points.
 * @param  density  2D array representing spatial density values.
 * @param  seed     Random number generator seed.
 * @param  bbox     Bounding box in which to generate the points (a,b,c,d =
 *                  xmin, xmax, ymin, ymax). Defaults to the unit square {0.f,
 * 1.f, 0.f, 1.f}.
 * @return          A Cloud containing the generated points.
 *
 * **Example**
 * @include ex_point_sampling.cpp
 *
 * **Result**
 * @image html ex_point_sampling0.png
 * @image html ex_point_sampling1.png
 * @image html ex_point_sampling2.png
 * @image html ex_point_sampling3.png
 */
Cloud random_cloud_distance(float            min_dist,
                            float            max_dist,
                            const Array     &density,
                            std::uint32_t    seed,
                            const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a random cloud of points with distances drawn from a
 * power-law distribution.
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
 * @return          A Cloud containing the generated points.
 *
 * **Example**
 * @include ex_point_sampling.cpp
 *
 * **Result**
 * @image html ex_point_sampling0.png
 * @image html ex_point_sampling1.png
 * @image html ex_point_sampling2.png
 * @image html ex_point_sampling3.png
 */
Cloud random_cloud_distance_power_law(
    float            dist_min,
    float            dist_max,
    float            alpha,
    std::uint32_t    seed,
    const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a random cloud of points with distances drawn from a Weibull
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
 * @return          A Cloud containing the generated points.
 *
 * **Example**
 * @include ex_point_sampling.cpp
 *
 * **Result**
 * @image html ex_point_sampling0.png
 * @image html ex_point_sampling1.png
 * @image html ex_point_sampling2.png
 * @image html ex_point_sampling3.png
 */
Cloud random_cloud_distance_weibull(
    float            dist_min,
    float            lambda,
    float            k,
    std::uint32_t    seed,
    const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a jittered grid cloud of points.
 *
 * Points are placed on a grid and then offset by a jitter amount, optionally
 * staggered for a more irregular pattern.
 *
 * @param  count         Number of points to generate.
 * @param  jitter_amount Maximum jitter to apply in each axis (x, y).
 * @param  stagger_ratio Ratio of staggering between consecutive rows or
 *                       columns.
 * @param  seed          Random number generator seed.
 * @param  bbox          Bounding box in which to generate the points (a,b,c,d =
 * xmin, xmax, ymin, ymax). Defaults to the unit square {0.f, 1.f, 0.f, 1.f}.
 * @return               A Cloud containing the generated points.
 *
 * **Example**
 * @include ex_point_sampling.cpp
 *
 * **Result**
 * @image html ex_point_sampling0.png
 * @image html ex_point_sampling1.png
 * @image html ex_point_sampling2.png
 * @image html ex_point_sampling3.png
 */
Cloud random_cloud_jittered(size_t           count,
                            const glm::vec2 &jitter_amount,
                            const glm::vec2 &stagger_ratio,
                            std::uint32_t    seed,
                            const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Filter a point cloud using rejection sampling based on a density mask.
 *
 * Each point in the input cloud is retained or discarded according to a
 * spatially varying probability derived from the provided 2D density mask. The
 * mask is mapped to the bounding box and evaluated as a continuous function
 * over the domain.
 *
 * @param density_mask 2D array defining the spatial density field over the
 *                     bounding box. Values should lie in \f$[0,1]\f$.
 * @param seed         Random seed used for reproducible rejection sampling.
 * @param bbox         Bounding box `(xmin, ymin, xmax, ymax)` defining the
 *                     spatial extent of the density mask in world coordinates.
 */
void rejection_filter_density(Cloud           &cloud,
                              const Array     &density_mask,
                              std::uint32_t    seed,
                              const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Scales the point coordinates in a cloud relative to a center point.
 *
 * Scales the (x, y) coordinates of all points in the cloud relative to a
 * specified center point while preserving each point's value (v). A scale
 * factor smaller than 1 shrinks the point distribution toward the center
 * (creating a margin from domain boundaries), while a scale factor greater than
 * 1 expands it.
 *
 * @param  cloud  Input point cloud.
 * @param  scale  Scale factor(s) for x and y axes.
 * @param  center Center of scaling (defaults to (0.5, 0.5)).
 * @return        Cloud  A new point cloud with scaled point coordinates.
 */
Cloud scale(const Cloud &cloud,
            glm::vec2    scale,
            glm::vec2    center = {0.5f, 0.5f});

Cloud scale(const Cloud &cloud, float scale, glm::vec2 center = {0.5f, 0.5f});

} // namespace hmap
