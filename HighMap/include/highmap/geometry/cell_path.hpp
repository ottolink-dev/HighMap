/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file cell_path.hpp
 * @brief Path class for manipulating and analyzing paths through array/grid
 * indices (i, j).
 * @copyright Copyright (c) 2026 Otto Link
 */
#pragma once
#include <vector>

#include <glm/glm.hpp>

#include "highmap/functions.hpp"

namespace hmap
{

/**
 * @brief Ordered path of grid cell indices.
 *
 * Stores a sequence of 2D integer coordinates (i, j) representing a path
 * through a grid or array.
 */
class CellPath
{
public:
  // ==========================================================================
  //  Constructors
  // ==========================================================================

  /**
   * @brief Construct a new CellPath object.
   */
  CellPath() = default;

  /**
   * @brief Construct a path by copying cell indices.
   * @param indices Ordered cell indices.
   */
  CellPath(const std::vector<glm::ivec2> &indices);

  /**
   * @brief Construct a path by moving cell indices.
   * @param indices Ordered cell indices.
   */
  CellPath(std::vector<glm::ivec2> &&indices) noexcept;

  // ==========================================================================
  //  Container interface
  // ==========================================================================

  /**
   * @brief Return an iterator to the beginning.
   */
  std::vector<glm::ivec2>::iterator begin() noexcept;

  /**
   * @brief Return a const iterator to the beginning.
   */
  std::vector<glm::ivec2>::const_iterator begin() const noexcept;

  /**
   * @brief Return an iterator to the end.
   */
  std::vector<glm::ivec2>::iterator end() noexcept;

  /**
   * @brief Return a const iterator to the end.
   */
  std::vector<glm::ivec2>::const_iterator end() const noexcept;

  /**
   * @brief Access the first element.
   */
  glm::ivec2 &front();

  /**
   * @brief Access the first element (const).
   */
  const glm::ivec2 &front() const;

  /**
   * @brief Access the last element.
   */
  glm::ivec2 &back();

  /**
   * @brief Access the last element (const).
   */
  const glm::ivec2 &back() const;

  /**
   * @brief Access element by index.
   * @param index Element index.
   */
  glm::ivec2 &operator[](size_t index);

  /**
   * @brief Access element by index (const).
   * @param index Element index.
   */
  const glm::ivec2 &operator[](size_t index) const;

  /**
   * @brief Clear all cell indices.
   */
  void clear() noexcept;

  /**
   * @brief Check whether the cell path is empty.
   * @return True if empty, false otherwise.
   */
  bool empty() const noexcept;

  /**
   * @brief Add a cell index to the end of the path.
   * @param idx Grid cell index to append.
   */
  void push_back(const glm::ivec2 &idx);

  /**
   * @brief Get the number of cells in the path.
   * @return Number of cell indices.
   */
  size_t size() const noexcept;

  /**
   * @brief Get the cell indices.
   * @return Mutable reference to the indices.
   */
  std::vector<glm::ivec2> &get_indices();

  /**
   * @brief Get the cell indices.
   * @return Const reference to the indices.
   */
  const std::vector<glm::ivec2> &get_indices() const;

private:
  std::vector<glm::ivec2> indices;
};

/**
 * @brief Rasterizes a segment between two integer points using Bresenham's
 * algorithm.
 *
 * Ensures full grid connectivity by filling all intermediate integer cells
 * between points a and b.
 *
 * @param out Output container receiving the rasterized points.
 * @param a   Start point.
 * @param b   End point.
 */
void add_line_bresenham(CellPath &out, glm::ivec2 a, glm::ivec2 b);

/**
 * @brief Applies transverse noise displacement along an integer path.
 *
 * The function computes a global direction from first to last point, derives a
 * perpendicular axis, and displaces each point using procedural noise.
 *
 * @param path        Input/output integer path to deform.
 * @param seed        Random seed for noise generation.
 * @param kw          Noise spatial frequency (scale).
 * @param amp         Maximum displacement amplitude.
 * @param shape       Optional domain scaling / bounds for noise sampling.
 * @param noise_type  Type of procedural noise to use.
 * @param octaves     Number of noise octaves for fractal noise.
 * @param weight      Blend factor between octaves.
 * @param persistence Amplitude decay per octave.
 * @param lacunarity  Frequency increase per octave.
 *
 * **Example**
 * @include ex_add_noise_ipath.cpp
 *
 * **Result**
 * @image html ex_add_noise_ipath.png
 */
void add_noise(CellPath         &path,
               std::uint32_t     seed,
               float             kw,
               float             amp,
               const glm::ivec2 &shape,
               NoiseType         noise_type = NoiseType::PERLIN,
               int               octaves = 8,
               float             weight = 0.7f,
               float             persistence = 0.5f,
               float             lacunarity = 2.f);

/**
 * @brief Ensures that an integer path is grid-adjacent and gap-free.
 *
 * Reconstructs the path so that every consecutive step moves to a neighboring
 * grid cell (8-connected), filling missing cells if needed.
 *
 * @param path Input/output path to repair and densify.
 */
void enforce_path_adjacency(CellPath &path);

/**
 * @brief Checks if a grid path has no gaps (8-connectivity).
 *
 * Verifies that each consecutive cell is adjacent using Moore neighborhood (max
 * |dx| <= 1 and max |dy| <= 1).
 *
 * @param  path Input path of grid cells.
 * @return      true if the path is fully 8-connected, false otherwise.
 */
bool is_path_adjacent(const CellPath &path);

} // namespace hmap