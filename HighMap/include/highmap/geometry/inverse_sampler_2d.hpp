/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file inverse_sampler_2d.hpp
 * @copyright Copyright (c) 2025
 */
#pragma once
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <random>
#include <utility>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @class InverseSampler2D
 * @brief Sequential 2D inverse transform point sampler based on a 2D density
 * array.
 *
 * Supports generating points one-by-one and dynamically updating density in
 * local bounding box or index sub-regions [imin, imax] x [jmin, jmax].
 */
class InverseSampler2D
{
public:
  // ==========================================================================
  //  Constructors
  // ==========================================================================

  InverseSampler2D() = default;

  /**
   * @brief Constructs an inverse sampler for a 2D density array.
   *
   * @param density 2D density array.
   * @param seed    Random number generator seed.
   * @param bbox    Bounding box (xmin, xmax, ymin, ymax). Defaults to unit
   *                square.
   */
  InverseSampler2D(const Array     &density,
                   std::uint32_t    seed = 0,
                   const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

  // ==========================================================================
  //  Sampling
  // ==========================================================================

  /**
   * @brief Samples a single point sequentially using the internal RNG.
   *
   * @return 2D point coordinates within the bounding box.
   */
  glm::vec2 sample();

  /**
   * @brief Samples a single point with user-provided uniform random numbers.
   *
   * @param  u  Uniform random number in [0, 1) for row selection.
   * @param  u2 Uniform random number in [0, 1) for column selection.
   * @return    2D point coordinates within the bounding box.
   */
  glm::vec2 sample(float u, float u2);

  /**
   * @brief Samples multiple points sequentially.
   *
   * @param  count Number of points to sample.
   * @return       A pair of coordinate vectors {x_coords, y_coords}.
   */
  std::array<std::vector<float>, 2> sample(size_t count);

  // ==========================================================================
  //  Updates & Modifiers
  // ==========================================================================

  /**
   * @brief Updates the density in the sub-region [imin, imax] x [jmin, jmax]
   * from another array.
   *
   * @param updated_array Array containing updated density values.
   * @param range_i       Column index range {imin, imax} (inclusive).
   * @param range_j       Row index range {jmin, jmax} (inclusive).
   */
  void update_region(const Array      &updated_array,
                     const glm::ivec2 &range_i,
                     const glm::ivec2 &range_j);

  /**
   * @brief Updates the density in the sub-region [imin, imax] x [jmin, jmax]
   * using a modifier function.
   *
   * @param updater Callback f(i, j, current_val) returning the new density
   * value.
   * @param range_i Column index range {imin, imax} (inclusive).
   * @param range_j Row index range {jmin, jmax} (inclusive).
   */
  void update_region(
      const std::function<float(int i, int j, float current)> &updater,
      const glm::ivec2                                        &range_i,
      const glm::ivec2                                        &range_j);

  /**
   * @brief Reinitializes the sampler with a new density array.
   *
   * @param density 2D density array.
   * @param seed    Random number generator seed.
   * @param bbox    Bounding box (xmin, xmax, ymin, ymax).
   */
  void reset(const Array     &density,
             std::uint32_t    seed = 0,
             const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

  // ==========================================================================
  //  Getters & Setters
  // ==========================================================================

  /**
   * @brief Gets the bounding box.
   */
  const glm::vec4 &get_bbox() const;

  /**
   * @brief Gets the current density array.
   */
  const Array &get_density() const;

  /**
   * @brief Gets the total integrated weight.
   */
  float get_total_weight() const;

  /**
   * @brief Sets the bounding box.
   */
  void set_bbox(const glm::vec4 &bbox);

  /**
   * @brief Sets the RNG seed.
   */
  void set_seed(std::uint32_t seed);

private:
  Array                                 density_;
  glm::vec4                             bbox_ = {0.f, 1.f, 0.f, 1.f};
  std::vector<float>                    row_sums_;
  std::vector<float>                    r_cdf_;
  std::mt19937                          rng_{0};
  std::uniform_real_distribution<float> dist_{0.f, 1.f};

  void rebuild_cdf();
  void update_row_cdf_from(int jmin);
};

} // namespace hmap
