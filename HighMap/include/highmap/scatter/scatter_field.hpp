/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file scatter_field.hpp
 * @copyright Copyright (c) 2026 Otto Link.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/geometry/cloud.hpp"
#include "highmap/scatter/scatter_item.hpp"

namespace hmap
{

/**
 * @enum ScatterShape
 * @brief Defines the geometric shape used when converting a ScatterField to a
 * heightmap.
 */
enum ScatterShape : int
{
  SCATTER_SHAPE_DISK,        ///< Circular dome / hemisphere profile.
  SCATTER_SHAPE_POLYGON,     ///< Irregular polygon profile.
  SCATTER_SHAPE_CONE,        ///< Linear cone profile.
  SCATTER_SHAPE_PYRAMID,     ///< Polygonal pyramid profile.
  SCATTER_SHAPE_SMOOTH_DOME, ///< Smooth cosine dome profile.
};

/**
 * @class ScatterField
 * @brief Container and processor for spatial ScatterItem collections.
 */
class ScatterField
{
public:
  // ==========================================================================
  //  Constructors
  // ==========================================================================

  /**
   * @brief Default constructor initializing an empty field.
   */
  ScatterField() = default;

  /**
   * @brief Constructs a ScatterField from a vector of items.
   * @param items Vector of ScatterItem instances.
   */
  ScatterField(const std::vector<ScatterItem> &items);

  /**
   * @brief Move-constructs a ScatterField from a vector of items.
   * @param items Rvalue vector of ScatterItem instances.
   */
  ScatterField(std::vector<ScatterItem> &&items) noexcept;

  /**
   * @brief Constructs a ScatterField from a Cloud of 2D points.
   *
   * If a point's value `v` is non-zero, it is used as the item's radius;
   * otherwise @p default_radius is assigned.
   *
   * @param cloud          Input point cloud.
   * @param class_id       Class identifier assigned to all imported items.
   * @param default_radius Default radius assigned if point value is 0.
   */
  ScatterField(const Cloud &cloud,
               uint32_t     class_id = 0,
               float        default_radius = HMAP_DEFAULT_SCATTER_RADIUS);

  virtual ~ScatterField() = default;

  // ==========================================================================
  //  Container Interface
  // ==========================================================================

  ScatterItem       &at(size_t index);
  const ScatterItem &at(size_t index) const;

  ScatterItem       &back();
  const ScatterItem &back() const;

  auto begin() noexcept
  {
    return items.begin();
  }
  auto begin() const noexcept
  {
    return items.begin();
  }
  auto cbegin() const noexcept
  {
    return items.cbegin();
  }

  size_t capacity() const noexcept
  {
    return items.capacity();
  }
  void clear() noexcept
  {
    items.clear();
  }

  ScatterItem *data() noexcept
  {
    return items.data();
  }
  const ScatterItem *data() const noexcept
  {
    return items.data();
  }

  template <typename... Args> ScatterItem &emplace_back(Args &&...args)
  {
    return items.emplace_back(std::forward<Args>(args)...);
  }

  bool empty() const noexcept
  {
    return items.empty();
  }

  auto end() noexcept
  {
    return items.end();
  }
  auto end() const noexcept
  {
    return items.end();
  }
  auto cend() const noexcept
  {
    return items.cend();
  }

  ScatterItem       &front();
  const ScatterItem &front() const;

  ScatterItem &operator[](size_t index)
  {
    return items[index];
  }
  const ScatterItem &operator[](size_t index) const
  {
    return items[index];
  }

  void push_back(const ScatterItem &item)
  {
    items.push_back(item);
  }
  void push_back(ScatterItem &&item)
  {
    items.push_back(std::move(item));
  }

  void reserve(size_t new_cap)
  {
    items.reserve(new_cap);
  }
  size_t size() const noexcept
  {
    return items.size();
  }

  // ==========================================================================
  //  Operations
  // ==========================================================================

  /**
   * @brief Densifies the scatter field by adding Voronoi vertices
   * (circumcenters of Delaunay triangles).
   *
   * @param default_radius Default radius assigned to newly added items.
   */
  void densify(float default_radius = HMAP_DEFAULT_SCATTER_RADIUS);

  /**
   * @brief Returns a new ScatterField containing only items matching the given
   * class_id.
   * @param  class_id Class identifier to filter by.
   * @return          ScatterField Filtered field.
   */
  ScatterField filter_by_class(uint32_t class_id) const;

  /**
   * @brief Computes the 2D bounding box enclosing all item positions.
   * @return glm::vec4 Bounding box as {xmin, xmax, ymin, ymax}.
   */
  glm::vec4 get_bbox() const;

  /**
   * @brief Returns a sorted list of unique class identifiers present in the
   * field.
   * @return std::vector<uint32_t> Unique class IDs.
   */
  std::vector<uint32_t> get_class_ids() const;

  /**
   * @brief Extracts the characteristic radius of all items.
   * @return std::vector<float> Vector of item radii.
   */
  std::vector<float> get_radius() const;

  /**
   * @brief Extracts the x-coordinate of all items.
   * @return std::vector<float> Vector of x-coordinates.
   */
  std::vector<float> get_x() const;

  /**
   * @brief Extracts the y-coordinate of all items.
   * @return std::vector<float> Vector of y-coordinates.
   */
  std::vector<float> get_y() const;

  /**
   * @brief Extracts the z-coordinate (elevation) of all items.
   * @return std::vector<float> Vector of z-coordinates.
   */
  std::vector<float> get_z() const;

  /**
   * @brief Slightly perturbs the 2D (x, y) coordinates of scatter items.
   *
   * @param dx   Maximum displacement along the x-axis.
   * @param dy   Maximum displacement along the y-axis (defaults to 0.0).
   * @param seed Random seed for reproducibility.
   */
  void perturb_positions(float dx, float dy = 0.0f, uint32_t seed = 0);

  /**
   * @brief Resolves overlapping items by discarding smaller items in collision.
   */
  void prune_collisions();

  /**
   * @brief Prunes items in this field that overlap with (fall within the
   * influence radius of) items in another field.
   *
   * The other field remains constant and unmodified. For each item in this
   * field, if its distance to any item in @p other is strictly less than the
   * sum of their radii, it is pruned.
   *
   * @param other Const reference to the reference ScatterField (e.g. rocks).
   *
   * @overload
   */
  void prune_collisions(const ScatterField &other);

  /**
   * @brief Prunes items in the field using density-based acceptance sampling
   * modulated to reach an approximate target retention ratio.
   *
   * @param density_mask 2D array defining the spatial density field.
   * @param target_ratio Approximate target fraction of items to retain in [0,
   *                     1].
   * @param seed         Random seed for reproducible pruning.
   * @param bbox         Bounding box defining the domain of the density mask.
   */
  void prune_density(const Array     &density_mask,
                     float            target_ratio = 0.8f,
                     uint32_t         seed = 0,
                     const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

  /**
   * @brief Relaxes/regularizes item positions using k-nearest neighbor
   * repulsion.
   *
   * Uses PointSampler's relaxation algorithm to reduce clustering and obtain a
   * more uniform or blue-noise-like spatial distribution.
   *
   * @param k_neighbors Number of nearest neighbors to consider for repulsion.
   * @param step_size   Step size per relaxation iteration.
   * @param iterations  Number of relaxation iterations.
   */
  void regularize_positions(size_t k_neighbors = 8,
                            float  step_size = 0.1f,
                            size_t iterations = 10);

  /**
   * @brief Reinforces spatial clustering of classes by iteratively assigning
   * each item the dominant class among its nearest neighbors.
   *
   * @param iterations   Number of smoothing/reinforcement iterations.
   * @param k_neighbors  Number of spatial nearest neighbors to query.
   * @param include_self If true, considers the item's own current class in the
   *                     majority vote.
   */
  void reinforce_class_clusters(size_t iterations = 2,
                                size_t k_neighbors = 4,
                                bool   include_self = true);

  /**
   * @brief Resolves overlapping items by iteratively moving colliding items
   * apart using Delaunay triangulation neighborhood and repulsive relaxation.
   *
   * Preserves item count and overall distribution while moving items until
   * overlap is resolved or max iterations are reached.
   *
   * @param iterations            Maximum number of relaxation iterations.
   * @param tolerance             Collision tolerance ratio in [0, 1) where 0
   *                              means no overlap allowed.
   * @param step_size             Relaxation step size / movement rate.
   * @param triangulation_substep Frequency (in iterations) to recompute the
   *                              Delaunay triangulation.
   */
  void resolve_collisions(size_t iterations = 20,
                          float  tolerance = 0.f,
                          float  step_size = 2.,
                          size_t triangulation_substep = 5);

  /**
   * @brief Resolves collisions with another static ScatterField by pushing
   * items in this field away from colliding items in the other field, using
   * joint Delaunay triangulation neighborhood and repulsive relaxation.
   *
   * The other field remains constant and unmodified.
   *
   * @param other                 Const reference to the reference ScatterField.
   * @param iterations            Maximum number of relaxation iterations.
   * @param tolerance             Collision tolerance ratio in [0, 1).
   * @param step_size             Relaxation step size.
   * @param triangulation_substep Frequency (in iterations) to recompute the
   *                              Delaunay triangulation.
   *
   * @overload
   */
  void resolve_collisions(const ScatterField &other,
                          size_t              iterations = 20,
                          float               tolerance = 0.f,
                          float               step_size = 2.f,
                          size_t              triangulation_substep = 5);

  /**
   * @brief Sets the elevation (z-coordinate) of all items by sampling a terrain
   * heightmap.
   *
   * @param elevation Terrain heightmap array.
   * @param bbox      Bounding box defining the domain of the elevation array.
   */
  void set_elevation_from_terrain(const Array     &elevation,
                                  const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

  /**
   * @brief Randomly shuffles class identifiers between neighboring items.
   *
   * @param ratio       Fraction of items to attempt class shuffling on in [0,
   *                    1].
   * @param k_neighbors Number of nearest spatial neighbors to consider.
   * @param seed        Random seed for reproducibility.
   */
  void shuffle_classes(float    ratio = 0.1f,
                       size_t   k_neighbors = 4,
                       uint32_t seed = 0);

  /**
   * @brief Converts the field into a Cloud of 2D points (x, y) with value set
   * to radius.
   * @return Cloud Point cloud representation.
   */
  Cloud to_cloud() const;

  /**
   * @brief Exports the local density map of the field as a 2D Array using 2D
   * Kernel Density Estimation (Gaussian splatting / KDE).
   *
   * @param  shape            Dimensions of the output density map {nx, ny}.
   * @param  sigma            Gaussian kernel standard deviation (bandwidth) in
   *                          world/domain coordinates. If <= 0, automatically
   *                          derived.
   * @param  class_id         Optional class identifier to compute density for a
   *                          specific class only.
   * @param  weighted_by_area If true, weights density by footprint area ($\pi
   *                          r^2$).
   * @param  bbox             Bounding box defining the domain {xmin, xmax,
   *                          ymin, ymax}.
   * @return                  Array 2D continuous density map array.
   */
  Array to_density_map(glm::ivec2              shape,
                       float                   sigma = 0.0f,
                       std::optional<uint32_t> class_id = std::nullopt,
                       bool                    weighted_by_area = false,
                       glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f}) const;

  /**
   * @brief Converts the scatter field into a heightmap by representing each
   * scattered object as a local elevation increase on GPU.
   *
   * @param  shape               Dimensions of the output heightmap {nx, ny}.
   * @param  shape_type          Object geometry shape (default
   *                             SCATTER_SHAPE_DISK).
   * @param  height_radius_ratio Ratio of object maximum height relative to its
   *                             radius (default 1.0f).
   * @param  class_id            Optional class identifier to filter items.
   * @param  p_rock_map          Optional output array for the rock map with
   *                             normalized amplitude in [0, 1].
   * @param  seed                Random seed for procedural shape variation
   *                             (e.g. polygon).
   * @param  bbox                Domain bounding box {xmin, xmax, ymin, ymax}.
   * @return                     Array               2D heightmap array
   *                             containing accumulated elevations.
   */
  Array to_heightmap(glm::ivec2              shape,
                     ScatterShape            shape_type = SCATTER_SHAPE_DISK,
                     float                   height_radius_ratio = 1.0f,
                     std::optional<uint32_t> class_id = std::nullopt,
                     Array                  *p_rock_map = nullptr,
                     uint32_t                seed = 0,
                     glm::vec4               bbox = {0.f, 1.f, 0.f, 1.f}) const;

  /**
   * @brief Exports the items to a CSV file (x, y, z, class_id, radius).
   * @param fname Output file path.
   */
  void to_csv(const std::string &fname) const;

  /**
   * @brief Renders the scatter field to an 8-bit RGB image buffer.
   *
   * @param  shape      Image dimensions {width, height}.
   * @param  background Optional background terrain array.
   * @param  bbox       Bounding box {xmin, xmax, ymin, ymax} of the domain.
   * @param  flip_y     If true, row 0 is top (screen/image space). If false,
   *                    row 0 is bottom (Cartesian domain space).
   * @return            std::vector<uint8_t> Interleaved RGB 8-bit image data
   *                    (size = width * height * 3).
   */
  std::vector<uint8_t> to_img_8bit(glm::ivec2   shape,
                                   const Array &background = {},
                                   glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f},
                                   bool         flip_y = false) const;

  /**
   * @brief Exports a visual representation of the scatter field as a PNG image.
   *
   * @param fname      Output PNG file path.
   * @param shape      Image dimensions {width, height}.
   * @param background Optional background terrain array.
   * @param bbox       Bounding box {xmin, xmax, ymin, ymax} of the domain.
   */
  void to_png(const std::string &fname,
              glm::ivec2         shape,
              const Array       &background = {},
              glm::vec4          bbox = {0.f, 1.f, 0.f, 1.f}) const;

  /**
   * @brief Returns a multi-line formatted summary string of the scatter field.
   * @return std::string Pretty-printed summary.
   */
  virtual std::string to_string() const;

protected:
  std::vector<ScatterItem> items = {}; ///< List of scatter item instances.
};

// ============================================================================
//  Functions (Alphabetically Sorted)
// ============================================================================

/**
 * @brief Merges multiple ScatterField instances into a single instance.
 *
 * @param  fields         List of ScatterField instances.
 * @param  merge_by_class If true, preserves existing class IDs. If false,
 *                        assigns a new unique class ID for each input field
 *                        corresponding to its index in @p fields.
 * @return                ScatterField   Merged scatter field.
 */
ScatterField merge_scatter_fields(const std::vector<ScatterField> &fields,
                                  bool merge_by_class = true);

/**
 * @brief Merges multiple ScatterField instances with associated bounding boxes
 * into a single instance, resolving overlaps by keeping items from earlier
 * fields in the input order.
 *
 * @param  fields         List of ScatterField instances.
 * @param  bboxs          Bounding box for each input field {xmin, xmax, ymin,
 *                        ymax}.
 * @param  merge_by_class If true, preserves existing class IDs. If false,
 *                        assigns a new unique class ID for each input field
 *                        corresponding to its index in @p fields.
 * @return                ScatterField   Merged scatter field.
 */
ScatterField merge_scatter_fields(const std::vector<ScatterField> &fields,
                                  const std::vector<glm::vec4>    &bboxs,
                                  bool merge_by_class = true);

/**
 * @brief Converts a ScatterField into a heightmap by representing each
 * scattered object as a local elevation increase on GPU.
 *
 * @param  field               Input ScatterField.
 * @param  shape               Dimensions of the output heightmap {nx, ny}.
 * @param  shape_type          Object geometry shape.
 * @param  height_radius_ratio Ratio of object maximum height relative to its
 *                             radius.
 * @param  class_id            Optional class identifier to filter items.
 * @param  p_rock_map          Optional output array for the rock map with
 *                             normalized amplitude in [0, 1].
 * @param  seed                Random seed for procedural shape variation.
 * @param  bbox                Domain bounding box {xmin, xmax, ymin, ymax}.
 * @return                     Array               2D heightmap array.
 */
Array scatter_field_to_heightmap(
    const ScatterField     &field,
    glm::ivec2              shape,
    ScatterShape            shape_type = SCATTER_SHAPE_DISK,
    float                   height_radius_ratio = 1.0f,
    std::optional<uint32_t> class_id = std::nullopt,
    Array                  *p_rock_map = nullptr,
    uint32_t                seed = 0,
    glm::vec4               bbox = {0.f, 1.f, 0.f, 1.f});

} // namespace hmap
