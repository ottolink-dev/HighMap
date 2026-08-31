/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file path.hpp
 * @brief Path class for manipulating and analyzing paths in 2D space.
 *
 * This file defines the `Path` class, which extends the `Cloud` class to
 * provide functionalities specific to paths. The `Path` class includes methods
 * for creating, manipulating, and analyzing paths, including smoothing the path
 * with various curves, resampling, and performing operations like meandering
 * and fractalizing.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once
#include "highmap/boundary.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/interpolate/interpolate1d.hpp"

namespace hmap
{

/**
 * @class Path
 * @brief Represents an ordered set of points in 2D, forming a polyline (open or
 * closed).
 *
 * The `Path` class extends the `Cloud` class to represent a sequence of points
 * in 2D space that are connected in a specific order, forming a polyline. The
 * polyline can be either open or closed, depending on whether the first and
 * last points are connected. This class provides methods for various geometric
 * operations, including path length calculation, point insertion, and more.
 *
 * **Example**
 * @include ex_path.cpp
 *
 * **Result**
 * @image html ex_path.png
 */
class Path : public Cloud
{
public:
  enum class EdgeDivisionMode : int
  {
    EDM_PER_EDGE,
    EDM_FULL_ARC,
  };

  // ==========================================================================
  //  Constructors
  // ==========================================================================

  /**
   * @brief Construct a new Path object with default properties. Initializes an
   * empty path with the `closed` property set to `false`.
   */
  Path() = default;

  /**
   * @brief Construct a new Path object with random positions and values.
   * Initializes a path with a specified number of points, random values, and
   * the option to be open or closed.
   * @param npoints Number of points to generate.
   * @param seed    Random seed number for generating random values.
   * @param bbox    Bounding box for random point generation.
   */
  Path(int npoints, std::uint32_t seed, glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f})
      : Cloud(npoints, seed, bbox){};

  /**
   * @brief Construct a new Path object based on a list of points. Initializes a
   * path with the specified points and an option to be open or closed.
   * @param points List of points defining the path.
   */
  Path(std::vector<Point> points) : Cloud(points){};

  /**
   * @brief Construct a new Path object based on `x` and `y` coordinates.
   * Initializes a path with the specified `x` and `y` coordinates and an option
   * to be open or closed.
   * @param x List of `x` coordinates for the points.
   * @param y List of `y` coordinates for the points.
   */
  Path(std::vector<float> x, std::vector<float> y) : Cloud(x, y){};

  /**
   * @brief Construct a new Path object based on `x`, `y` coordinates, and
   * values. Initializes a path with the specified `x` and `y` coordinates,
   * associated values, and an option to be open or closed.
   * @param x List of `x` coordinates for the points.
   * @param y List of `y` coordinates for the points.
   * @param v List of values associated with the points.
   */
  Path(std::vector<float> x, std::vector<float> y, std::vector<float> v)
      : Cloud(x, y, v){};

  /**
   * @brief Construct a point cloud from grid indices mapped to a bounding box.
   *
   * Each index is normalized by the grid shape and remapped to the given
   * bounding box. The resulting points are stored as 3D positions (z = 1).
   *
   * @param indices Input grid indices.
   * @param shape   Grid dimensions.
   * @param bbox    Bounding box (xmin, xmax, ymin, ymax). Default is {0.f, 1.f,
   *                0.f, 1.f}.
   */
  Path(const std::vector<glm::ivec2> &indices,
       const glm::ivec2              &shape,
       const glm::vec4               &bbox = {0.f, 1.f, 0.f, 1.f})
      : Cloud(indices, shape, bbox){};

  /**
   * @brief  Constructs a new Path object from lists of xyz data as glm::vec3.
   * */
  Path(const std::vector<glm::vec3> &xyv) : Cloud(xyv){};

  // ==========================================================================
  //  Accessors
  // ==========================================================================

  /**
   * @brief Get the arc length of the path.
   *
   * The arc length is the cumulative distance along the path, normalized to the
   * range [0, 1]. This represents the distance traveled from the start to each
   * point along the path as a fraction of the total path length.
   *
   * @return std::vector<float> Vector of arc length values, where each entry
   *         corresponds to a point on the path and represents the normalized
   *         distance from the start of the path to that point.
   */
  std::vector<float> get_arc_length() const;

  /**
   * @brief Get the cumulative distance of the path.
   *
   * This method computes the cumulative distance along the path, which is the
   * total distance traveled up to each point on the path. It accumulates the
   * distances from the start of the path to each point.
   *
   * @return std::vector<float> Vector of cumulative distance values, where each
   *         entry represents the distance from the start of the path to the
   *         respective point.
   */
  std::vector<float> get_cumulative_distance() const;

  /**
   * @brief Computes the signed curvature at each point of the path.
   * @param  normalized If true, scales values to [-1, 1] using the maximum
   *                    absolute curvature.
   * @return            Vector of curvature values.
   */
  std::vector<float> get_curvature(bool normalized = false) const;

  /**
   * @brief Compute midpoints of path edges.
   *
   * Returns the midpoint of each segment between consecutive points.
   *
   * @return Vector of edge center points.
   */
  std::vector<Point> get_edge_centers() const;

  /**
   * @brief Computes unit normals at each point of the path.
   * @return Vector of normal vectors.
   */
  std::vector<glm::vec2> get_normals() const;

  /**
   * @brief Computes unit tangents at each point of the path.
   * @return Vector of tangent vectors.
   */
  std::vector<glm::vec2> get_tangents() const;

  /**
   * @brief Get the values assigned to the points on the path.
   *
   * This method retrieves the values assigned to each point on the path. These
   * values can represent any attribute associated with the points, such as
   * color, intensity, or other metrics.
   *
   * @return std::vector<float> Vector of values assigned to the points on the
   *         path.
   */
  std::vector<float> get_values() const;

  /**
   * @brief Get the `x` coordinates of the points on the path.
   *
   * This method retrieves the `x` coordinates of all points in the path. Each
   * value in the returned vector corresponds to the `x` coordinate of a point
   * along the path.
   *
   * @return std::vector<float> Vector of `x` values of the points on the path.
   */
  std::vector<float> get_x() const;

  /**
   * @brief Get the coordinates of the points as a single vector.
   *
   * This method returns the coordinates of the points in the path as a single
   * vector, where the coordinates are interleaved: `(x0, y0, x1, y1, ...)`.
   * This format is useful for operations that require a flat representation of
   * the coordinates.
   *
   * @return std::vector<float> Vector of interleaved `x` and `y` coordinates of
   *         the points on the path.
   */
  std::vector<float> get_xy() const;

  /**
   * @brief Get the `y` coordinates of the points on the path.
   *
   * This method retrieves the `y` coordinates of all points in the path. Each
   * value in the returned vector corresponds to the `y` coordinate of a point
   * along the path.
   *
   * @return std::vector<float> Vector of `y` values of the points on the path.
   */
  std::vector<float> get_y() const;

  /**
   * @brief Check whether the path is closed.
   * @return True if the path is closed, false otherwise.
   */
  bool is_closed() const;

  /**
   * @brief Set whether the path is closed.
   * @param new_value True to close the path, false to keep it open.
   */
  void set_closed(bool new_value);

  // ==========================================================================
  //  Basic Ops
  // ==========================================================================

  /**
   * @brief Clear the path data.
   *
   * This method removes all points and associated data from the path,
   * effectively resetting it to an empty state.
   */
  void clear();

  /**
   * @brief Enforces monotonicity on the values of the points in the path.
   *
   * This method adjusts the `v` values of the points in the path to ensure that
   * they are either monotonically decreasing or increasing, based on the input
   * parameter.
   *
   * @param decreasing If true, enforces a monotonically decreasing order for
   * the values. If false, enforces a monotonically increasing order for the
   * values.
   *
   * @note This method modifies the path in place.
   */
  void enforce_monotonic_values(bool decreasing = true);

  /**
   * @brief Reorder points using a nearest neighbor search.
   *
   * This method reorders the points in the path to minimize the total distance
   * by performing a nearest neighbor search starting from the specified
   * `start_index`. This approach is useful for optimizing the path or improving
   * its order.
   *
   * @param start_index Index of the starting point for the nearest neighbor
   *                    search. Default is 0.
   */
  void reorder_nns(int start_index = 0);

  /**
   * @brief Reverse the order of points in the path.
   *
   * This method reverses the sequence of points in the path, which can be
   * useful for various applications, such as changing the direction of the path
   * traversal.
   */
  void reverse();

  // ==========================================================================
  //  Sampling
  // ==========================================================================

  /**
   * @brief Divide the path by adding a point between each pair of consecutive
   * points.
   *
   * This method adds new points in the middle of each segment of the path to
   * create a denser set of points along the path. This is useful for increasing
   * the resolution of the path.
   */
  void divide();

  /**
   * @brief Resample the path so that there is at least one point per pixel
   * for a given grid shape and domain bounding box.
   *
   * @param shape      Grid dimensions (width, height).
   * @param bbox       Bounding box of the domain (xmin, xmax, ymin, ymax).
   * @param itp_method Interpolation method used to resample the path.
   */
  void resample_by_grid_resolution(
      glm::ivec2            shape,
      glm::vec4             bbox = {0.f, 1.f, 0.f, 1.f},
      InterpolationMethod1D itp_method = InterpolationMethod1D::LINEAR);

  /**
   * @brief Resample the path to achieve an approximately constant distance
   * between points.
   *
   * This method adjusts the points in the path to ensure that the distance
   * between each consecutive point is approximately equal to the specified
   * `delta`. This is useful for creating a path with evenly spaced points.
   *
   * @param delta      Target distance between consecutive points.
   * @param itp_method Interpolation method used to resample the path.
   */
  void resample_by_spacing(
      float                 delta,
      InterpolationMethod1D itp_method = InterpolationMethod1D::LINEAR);

  /**
   * @brief Resamples the path using cubic interpolation along arc length.
   *
   * The path is reparameterized by arc length and interpolated using cubic
   * splines for x, y, and value components.
   *
   * @param npoints    Number of points to resample.
   * @param itp_method Interpolation method to use (default is CUBIC).
   *
   * **Example**
   * @include ex_path_resample.cpp
   *
   * **Result**
   * @image html ex_path_resample.png
   */
  void resample_interp(
      int                   npoints,
      InterpolationMethod1D itp_method = InterpolationMethod1D::CUBIC);

  /**
   * @brief Resample the path to achieve fairly uniform distance between
   * consecutive points.
   *
   * This method adjusts the path so that the distance between each consecutive
   * point is as uniform as possible. It redistributes the points to ensure more
   * even spacing along the path.
   *
   * @param itp_method Interpolation method used for resampling.
   */
  void resample_uniform(
      InterpolationMethod1D itp_method = InterpolationMethod1D::LINEAR);

  /**
   * @brief Sample the path using normalized arc-length parameterization.
   *
   * Interpolates position and value at parameter @p t in [0, 1]. Optionally
   * uses a precomputed cumulative arc-length vector and returns the local unit
   * tangent at the sampled point.
   *
   * @param  t         Normalized arc-length parameter.
   * @param  p_arc     Optional cumulative arc-length vector to avoid
   *                   recomputation.
   * @param  p_tangent Optional output pointer receiving the unit tangent.
   * @return           Sampled position (x, y) and interpolated value (z).
   */
  glm::vec3 sample_at(float                     t,
                      const std::vector<float> *p_arc = nullptr,
                      glm::vec2                *p_tangent = nullptr) const;

  /**
   * @brief Subsample the path by keeping only every n-th point.
   *
   * This method reduces the number of points in the path by retaining only
   * every 'step'-th point, effectively subsampling the path. This can be useful
   * for simplifying the path or reducing data size.
   *
   * @param step The interval of points to keep. For example, a step of 2 will
   *             keep every second point.
   */
  void subsample(int step);

  // ==========================================================================
  //  Conversion / IO
  // ==========================================================================

  /**
   * @brief Project path points to an array.
   *
   * This method projects the points of the path onto a 2D array, filling the
   * array based on the path's points. Optionally, the contour of the path can
   * be filled using flood fill if the `filled` parameter is set to true.
   *
   * **Example**
   * @include ex_path_to_array.cpp
   *
   * @param array  The array to which the path points will be projected.
   * @param bbox   Bounding box defining the domain of the array.
   * @param filled Boolean flag indicating whether to perform flood filling of
   * the path's contour.
   */
  void to_array(Array    &array,
                glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f},
                bool      filled = false) const;

  /*! @brief See hmap::to_array */
  Array to_array(glm::ivec2 shape,
                 glm::vec4  bbox = {0.f, 1.f, 0.f, 1.f},
                 bool       filled = false) const;

  /*! @brief See hmap::to_array */
  void to_array_mask(Array    &array,
                     glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f},
                     bool      filled = false) const;

  /**
   * @brief Export path as PNG image file.
   *
   * This method generates a PNG image representing the path and saves it to the
   * specified file. The resolution of the image can be adjusted using the
   * `shape` parameter.
   *
   * **Example**
   * @include ex_path_to_png.cpp
   *
   * @param fname The filename for the output PNG image.
   * @param shape Resolution of the image, specified as width and height.
   * Default is {512, 512}.
   */
  void to_png(std::string fname, glm::ivec2 shape = {512, 512});

private:
  enum class PathClosure : int
  {
    PT_OPEN,
    PT_CLOSE,
  } path_closure = PathClosure::PT_OPEN;
};

// ==========================================================================
//  Functions
// ==========================================================================

/**
 * @brief Smooth the path using Bezier curves.
 *
 * This method applies Bezier curve smoothing to the path. The
 * `curvature_ratio` controls the amount of curvature applied, with typical
 * values in the range of [-1, 1], where positive values generally result in
 * more pronounced curvature. The `edge_divisions` parameter determines the
 * number of subdivisions per edge to create a smoother curve.
 *
 * **Example**
 * @include ex_path_bezier.cpp
 *
 * **Result**
 * @image html ex_path_bezier.png
 *
 * @param curvature_ratio Amount of curvature, usually in the range [-1, 1],
 *                        with positive values resulting in more curvature.
 * @param edge_divisions  Number of subdivisions per edge to achieve smooth
 *                        curves.
 * @param edm             The mode for dividing edges. Default is
 *                        Path::EdgeDivisionMode::EDM_PER_EDGE.
 *
 * See unit tests: @ref test_splines.cpp
 */
Path bezier(const Path            &path,
            float                  curvature_ratio = 0.3f,
            int                    edge_divisions = 10,
            Path::EdgeDivisionMode edm = Path::EdgeDivisionMode::EDM_PER_EDGE);

/**
 * @brief Smooth the path using Bezier curves (alternative method).
 *
 * This alternative method applies Bezier curve smoothing to the path. The
 * `curvature_ratio` parameter affects the curvature amount, similar to the
 * `bezier` method. The `edge_divisions` parameter specifies how finely each
 * edge is divided to create a smoother curve.
 *
 * **Example**
 * @include ex_path_bezier_round.cpp
 *
 * **Result**
 * @image html ex_path_bezier_round.png
 *
 * @param curvature_ratio Amount of curvature, typically within [-1, 1], with
 *                        positive values for increased curvature.
 * @param edge_divisions  Number of edge subdivisions for smoothness.
 * @param edm             The mode for dividing edges. Default is
 *                        Path::EdgeDivisionMode::EDM_PER_EDGE.
 *
 * See unit tests: @ref test_splines.cpp
 */
Path bezier_round(
    const Path            &path,
    float                  curvature_ratio = 0.3f,
    int                    edge_divisions = 10,
    Path::EdgeDivisionMode edm = Path::EdgeDivisionMode::EDM_PER_EDGE);

/**
 * @brief Smooth the path using B-Spline curves.
 *
 * This method smooths the path using B-Spline curves, which provide a smooth
 * curve that passes through the control points. The `edge_divisions`
 * parameter defines the number of subdivisions per edge for achieving smooth
 * curves.
 *
 * **Important**: This function does not correctly handle closed polylines
 * (circular contours). If the path is closed, the smoothing may not correctly
 * close the loop, potentially leaving a gap between the start and end points.
 *
 * **Example**
 * @include ex_path_bspline.cpp
 *
 * **Result**
 * @image html ex_path_bspline.png
 *
 * @param edge_divisions Number of subdivisions per edge to achieve a smooth
 *                       B-Spline curve.
 *
 * @warning This function does not correctly handle closed polylines.
 * @param edm            Mode for dividing edges. Default is EDM_PER_EDGE.
 *
 * See unit tests: @ref test_splines.cpp
 */
Path bspline(const Path            &path,
             int                    edge_divisions = 10,
             Path::EdgeDivisionMode edm = Path::EdgeDivisionMode::EDM_PER_EDGE);

/**
 * @brief Smooth the path using Catmull-Rom curves.
 *
 * This method applies Catmull-Rom curve smoothing to the path. Catmull-Rom
 * splines are interpolating splines that pass through each control point. The
 * `edge_divisions` parameter determines the number of subdivisions per edge for
 * smoothing.
 *
 * @param edge_divisions Number of edge subdivisions to create a smooth
 *                       Catmull-Rom curve.
 * @param edm            Mode for dividing edges. Default is EDM_PER_EDGE.
 *
 * **Example**
 * @include ex_path_catmullrom.cpp
 *
 * **Result**
 * @image html ex_path_catmullrom.png
 *
 * See unit tests: @ref test_splines.cpp
 */
Path catmullrom(
    const Path            &path,
    int                    edge_divisions = 10,
    Path::EdgeDivisionMode edm = Path::EdgeDivisionMode::EDM_PER_EDGE);

/**
 * @brief Smooth the path using De Casteljau curves.
 *
 * This function smooths a path by applying De Casteljau's algorithm to generate
 * intermediate points along the path, effectively creating a Bézier curve that
 * approximates the original path. The path is divided into segments, and the De
 * Casteljau algorithm is applied to each segment, resulting in a smooth curve.
 *
 * The parameter `edge_divisions` controls the number of divisions
 * (sub-segments) created along each segment of the path. A higher number of
 * divisions will result in a smoother curve, but will also increase the
 * computational cost.
 *
 * @param edge_divisions The number of divisions for each edge of the path.
 *                       Default is 10, which provides a balanced level of
 *                       smoothing.
 * @param edm            Mode for dividing edges (default is EDM_PER_EDGE).
 *
 * **Example**
 * @include ex_path_decasteljau.cpp
 *
 * **Result**
 * @image html ex_path_decasteljau.png
 *
 * See unit tests: @ref test_splines.cpp
 */
Path decasteljau(
    const Path            &path,
    int                    edge_divisions = 10,
    Path::EdgeDivisionMode edm = Path::EdgeDivisionMode::EDM_PER_EDGE);

/**
 * @brief Simplifies the current path using the Visvalingam-Whyatt algorithm.
 *
 * This method reduces the number of points in the path to the specified target,
 * `n_points_target`, while preserving the overall shape. It calculates the area
 * of triangles formed by consecutive points and removes points corresponding to
 * the smallest areas iteratively.
 *
 * @param path            The input path to be decimated.
 * @param n_points_target The desired number of points to retain in the path. If
 *                        the current number of points is less than
 *                        `n_points_target` or the path contains fewer than 3
 * points, the method returns without modifying the path.
 *
 * **Example**
 * @include ex_path_decimate.cpp
 *
 * **Result**
 * @image html ex_path_decimate.png
 *
 * See unit tests: @ref test_decimate_vw.cpp
 */
Path decimate_vw(const Path &path, int n_points_target = 3);

/**
 * @brief Applies fractalization to the path by adding points and randomly
 * displacing their positions.
 *
 * This method enhances the complexity of a path by iteratively adding new
 * points between existing ones and displacing them using Gaussian noise. The
 * process can simulate natural phenomena like terrain generation or random walk
 * paths. The number of iterations determines the level of detail added to the
 * path.
 *
 * @param path          The input path to be fractalized.
 * @param iterations    Number of iterations to apply the fractalization
 *                      process.
 * @param seed          Seed value for random number generation, ensuring
 *                      reproducibility.
 * @param sigma         Standard deviation of the Gaussian displacement,
 *                      relative to the distance between points.
 * @param orientation   Determines the displacement direction: `0` for random,
 * `1` for inflation, `-1` for deflation.
 * @param persistence   Factor that adjusts the noise intensity across
 *                      iterations.
 * @param control_field Optional pointer to an array that locally modifies the
 *                      displacement amplitude.
 * @param bbox          Bounding box that defines the valid area for the control
 *                      field's influence.
 * @param bounded       If true, clamps each newly displaced midpoint to the
 *                      bounding box defined by the initial edge it belongs to.
 *                      Defaults to false.
 *
 * **Example**
 * @include ex_path_fractalize.cpp
 *
 * **Result**
 * @image html ex_path_fractalize.png
 */
Path fractalize(const Path   &path,
                int           iterations,
                std::uint32_t seed,
                float         sigma = 0.2f,
                int           orientation = 0,
                float         persistence = 1.f,
                Array        *p_control_field = nullptr,
                glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f},
                bool          bounded = false);

/**
 * @brief Applies uniform fractalization along the path by first resampling the
 * path uniformly according to a target edge spacing (or number of samples) so
 * that all segments have equal base length, and then applying relative midpoint
 * displacement (scaled by sigma relative to segment length).
 *
 * @param path          The input path.
 * @param iterations    Number of recursive midpoint displacement levels.
 * @param seed          Seed for random number generation.
 * @param sigma         Standard deviation of the Gaussian displacement relative
 *                      to the resampled segment length.
 * @param spacing       Target uniform distance between consecutive points prior
 *                      to midpoint displacement (<= 0 to auto-compute from
 * total length and path size).
 * @param orientation   `0` for random lateral, `1` for positive normal, `-1`
 *                      for negative normal.
 * @param persistence   Roughness decay factor applied to sigma each level
 *                      (e.g., 1.0 or 0.5).
 * @param p_ctrl_array  Optional control array to modulate local amplitude.
 * @param bbox          Domain bounding box.
 * @param bounded       If true, clamps displaced midpoints to initial edge
 *                      bounding boxes.
 * @return              Fractalized path with uniform detail density.
 */
Path fractalize_uniform(const Path   &path,
                        int           iterations,
                        std::uint32_t seed,
                        float         sigma = 0.2f,
                        float         spacing = 0.f,
                        int           orientation = 0,
                        float         persistence = 1.f,
                        Array        *p_ctrl_array = nullptr,
                        glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f},
                        bool          bounded = false);

/**
 * @brief Inflate (offset) a path along its normals using curvature.
 *
 * Each point is displaced along its normal direction by an amount proportional
 * to local curvature and the given radius.
 *
 * @param  path     Input path
 * @param  radius   Inflation radius (scales displacement)
 * @param  resample If true, resample output to preserve spatial resolution
 *
 * @return          Inflated path
 *
 * **Example**
 * @include ex_path_inflate.cpp
 *
 * **Result**
 * @image html ex_path_inflate.png
 */
Path inflate(const Path &path, float strength, bool resample = true);

/**
 * @brief Add "meanders" to the path.
 *
 * This method introduces meandering effects to the path by adding random
 * deviations. The amplitude of the meanders is controlled by the `ratio`
 * parameter, while the `noise_ratio` controls the amount of randomness. The
 * `seed` parameter is used to initialize the random number generator, ensuring
 * reproducibility. The `iterations` parameter defines how many times the
 * meandering process is applied, and `edge_divisions` controls how finely each
 * edge is subdivided during the meandering.
 *
 * @param path           The input path to be meanderized.
 * @param ratio          Amplitude ratio of the meanders. Typically a positive
 *                       value.
 * @param noise_ratio    Ratio of randomness introduced during meandering.
 *                       Default is 0.1.
 * @param seed           Seed for random number generation. Default is 1.
 * @param iterations     Number of iterations to apply meandering. Default is 1.
 * @param edge_divisions Number of sub-divisions of each edge. Default is 10.
 * @param edm            The mode for dividing edges during the meandering
 *                       process. Default is EDM_PER_EDGE.
 *
 * **Example**
 * @include ex_path_meanderize.cpp
 *
 * **Result**
 * @image html ex_path_meanderize.png
 */
Path meanderize(
    const Path            &path,
    float                  ratio,
    float                  noise_ratio = 0.1f,
    std::uint32_t          seed = 1,
    int                    iterations = 1,
    int                    edge_divisions = 10,
    Path::EdgeDivisionMode edm = Path::EdgeDivisionMode::EDM_PER_EDGE);

/**
 * @brief Compute a distance field from a point path.
 *
 * For each grid cell, computes the distance to the nearest point in the path,
 * optionally applying domain warping using noise fields.
 *
 * @param  path       Input point path.
 * @param  shape      Output array dimensions.
 * @param  bbox_array Bounding box of the output domain.
 * @param  p_noise_x  Optional x-direction noise (domain warp).
 * @param  p_noise_y  Optional y-direction noise (domain warp).
 * @return            Array of distances to the nearest point.
 */
Array path_sdf_to_array(const Path  &path,
                        glm::ivec2   shape,
                        glm::vec4    bbox_array = {0.f, 1.f, 0.f, 1.f},
                        const Array *p_noise_x = nullptr,
                        const Array *p_noise_y = nullptr);

/**
 * @brief Removes geometric loops in a 2D path caused by self-intersections.
 *
 * This function iteratively detects self-intersections in the path and removes
 * the segments forming loops, optionally inserting the exact intersection
 * points.
 *
 * @param  path The input path as a vector of 2D points.
 * @return      Path The simplified path with self-intersecting loops removed.
 *
 * **Example**
 * @include ex_path_remove_geometry_loops.cpp
 *
 * **Result**
 * @image html ex_path_remove_geometry_loops.png
 */
Path remove_geometric_loops(const Path &path);

/**
 * @brief Applies a smoothing operation to the path points using a moving
 * average filter.
 *
 * This method smooths the path points based on a specified number of
 * neighboring points, an averaging intensity, and an inertia factor. The
 * smoothing involves calculating the average of neighboring points within a
 * range defined by `navg`, and then applying an intensity-based weighted
 * average to blend the original and smoothed values. Additionally, an inertia
 * effect can be applied to gradually adjust point positions based on previous
 * points.
 *
 * @param path                The input path to be smoothed.
 * @param navg                Number of neighboring points to consider on each
 *                            side of the current point during the smoothing
 *                            process. Higher values result in broader
 *                            smoothing.
 * @param averaging_intensity The weight given to the averaged points. A value
 *                            of 1.0 applies full intensity, resulting in a
 *                            complete averaging of the neighboring points.
 *                            Lower values retain more of the original point's
 *                            position.
 * @param inertia             The factor by which each point is influenced by
 *                            its previous point after the initial smoothing
 *                            pass. A value of 0 has no inertia effect, while a
 *                            higher value blends the current point's position
 *                            with that of the preceding point, creating a
 *                            trailing effect.
 *
 * **Example**
 * @include ex_path_smooth.cpp
 *
 * **Result**
 * @image html ex_path_smooth.png
 */
Path smooth(const Path &path,
            int         navg = 1,
            float       averaging_intensity = 1.f,
            float       inertia = 0.f);

/**
 * @brief Generates a continuous, non-self-intersecting squiggle curve by
 * applying the squiggle recursive triangle subdivision algorithm to each edge
 * of the input path.
 *
 * Implements the squiggle curve construction from Prusinkiewicz et al.,
 * "Generating Mountainous Terrain". Constrained within domain triangles for
 * each edge, ensuring continuity and absence of self-intersections.
 *
 * @param  path         The input path to apply the squiggle algorithm to.
 * @param  iterations   Number of recursive subdivision levels per edge (e.g., 2
 *                      to 6).
 * @param  seed         Random seed number for deterministic edge selection.
 * @param  height_ratio Height of the bounding triangle relative to edge length
 *                      (default 0.5).
 * @param  orientation  Displacement direction: 0 to alternate sides, 1 for
 *                      positive normal, -1 for negative normal.
 * @param  p_weights    Optional pointer to an Array defining spatial weights to
 *                      bias path progression.
 * @param  p_mask       Optional pointer to a binary/mask Array (1 = allowed, 0
 *                      = forbidden) preventing the path from entering forbidden
 * areas.
 * @param  bbox         Bounding box for spatial sampling.
 * @return              Generated continuous Path.
 *
 * **Example**
 * @include ex_path_squiggle.cpp
 *
 * **Result**
 * @image html ex_path_squiggle.png
 *
 * Reference: Prusinkiewicz et al., Generating Mountainous Terrain (1993).
 */
Path squiggle(const Path   &path,
              int           iterations = 4,
              std::uint32_t seed = 0,
              float         height_ratio = 0.5f,
              int           orientation = 0,
              const Array  *p_weights = nullptr,
              const Array  *p_mask = nullptr,
              glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates squiggle paths with secondary tributary branches along each
 * edge.
 *
 * @param  path               The input path.
 * @param  iterations         Number of recursive subdivision levels.
 * @param  seed               Random seed number.
 * @param  branch_probability Probability of initiating a tributary branch at
 *                            subdivision steps.
 * @param  max_branches       Maximum number of secondary branches.
 * @param  height_ratio       Height of the bounding triangle relative to edge
 *                            length.
 * @param  orientation        Displacement direction: 0 to alternate sides, 1
 *                            for positive normal, -1 for negative normal.
 * @param  p_weights          Optional spatial weights Array.
 * @param  p_mask             Optional mask Array.
 * @param  bbox               Bounding box for array sampling.
 * @return                    Vector of Paths where the first Path is the main
 *                            trunk and subsequent Paths are tributary branches.
 */
std::vector<Path> squiggle_branches(const Path   &path,
                                    int           iterations = 4,
                                    std::uint32_t seed = 0,
                                    float         branch_probability = 0.3f,
                                    int           max_branches = 8,
                                    float         height_ratio = 0.5f,
                                    int           orientation = 0,
                                    const Array  *p_weights = nullptr,
                                    const Array  *p_mask = nullptr,
                                    glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

// ==========================================================================
//  Verification Functions
// ==========================================================================

/**
 * @brief Asserts that the start and end points of two paths are within a
 * specified tolerance.
 *
 * @param path1   The first path to compare.
 * @param path2   The second path to compare.
 * @param tol     Tolerance for comparing the start and end points. Default is
 *                1e-6f.
 * @param verbose If true, prints detailed information during comparison.
 *                Default is false.
 */
bool assert_start_end_points(const Path &path1,
                             const Path &path2,
                             float       tol = 1e-6f,
                             bool        verbose = false);

/**
 * @brief Calculate the chamfer distance between two paths.
 *
 * @param a The first path.
 * @param b The second path.
 */
float chamfer_distance(const Path &a, const Path &b);

/**
 * @brief Check if a path contains duplicate points within a given tolerance.
 *
 * @param  path The path to check for duplicates.
 * @param  tol  Tolerance for considering two points as the same. Default is
 *              1e-6.
 * @return      true If there are duplicate points in the path.
 * @return      false If there are no duplicate points in the path.
 */
bool has_duplicates(const Path &path, float tol = 1e-6f);

} // namespace hmap
