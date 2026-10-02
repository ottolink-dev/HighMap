/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file forest_growth.hpp
 * @copyright Copyright (c) 2025 Otto Link.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/flora/forest.hpp"

namespace hmap
{

/**
 * @struct InteractionMatrix
 * @brief Generic S x S matrix representing pair-wise interaction properties
 * between species.
 */
struct InteractionMatrix
{
  size_t             size = 0;    ///< Number of species (matrix dimension).
  std::vector<float> values = {}; ///< Row-major matrix elements.

  // ==========================================================================
  //  Constructors
  // ==========================================================================

  InteractionMatrix() = default;

  /**
   * @brief Constructs an InteractionMatrix of size num_species x
   * num_species.
   * @param num_species Number of species.
   * @param default_val Initial value for all elements.
   */
  explicit InteractionMatrix(size_t num_species, float default_val = 0.0f);

  // ==========================================================================
  //  Accessors & Mutators
  // ==========================================================================

  /**
   * @brief Fills the entire matrix with a constant value.
   * @param value Value to fill.
   */
  void fill(float value);

  /**
   * @brief Gets the value at entry (s1, s2).
   * @param  s1 Row species index.
   * @param  s2 Column species index.
   * @return    float Matrix value.
   */
  float get(size_t s1, size_t s2) const;

  /**
   * @brief Sets the value at entry (s1, s2).
   * @param s1    Row species index.
   * @param s2    Column species index.
   * @param value Matrix value.
   */
  void set(size_t s1, size_t s2, float value);

  /**
   * @brief Sets symmetric values at (s1, s2) and (s2, s1).
   * @param s1    First species index.
   * @param s2    Second species index.
   * @param value Matrix value.
   */
  void set_symmetric(size_t s1, size_t s2, float value);

  // ==========================================================================
  //  Factory Methods
  // ==========================================================================

  /**
   * @brief Creates a matrix with specified diagonal values and constant
   * off-diagonal values.
   * @param  diag         Vector of diagonal elements.
   * @param  off_diag_val Off-diagonal value.
   * @return              InteractionMatrix Resulting matrix.
   */
  static InteractionMatrix diagonal(const std::vector<float> &diag,
                                    float off_diag_val = 0.0f);

  /**
   * @brief Creates an interaction matrix computed from individual species
   * radii.
   *
   * Each entry (s1, s2) is set to multiplier * (radii[s1] + radii[s2]).
   *
   * @param  radii      Vector of species radii.
   * @param  multiplier Scaling factor applied to the sum of radii
   *                    (defaults to 2.0f).
   * @return            InteractionMatrix Resulting matrix.
   */
  static InteractionMatrix from_radii(const std::vector<float> &radii,
                                      float multiplier = 2.0f);

  /**
   * @brief Creates a matrix initialized with random values perturbed
   * around 1.0.
   *
   * Each entry is sampled uniformly in [1.0 - random_offset, 1.0 +
   * random_offset].
   *
   * @param  num_species   Number of species.
   * @param  seed          Random seed.
   * @param  random_offset Maximum random variation relative to 1.0.
   * @param  symmetric     Whether to enforce symmetric values (A_ij =
   *                       A_ji).
   * @return               InteractionMatrix Resulting matrix.
   */
  static InteractionMatrix random(size_t   num_species,
                                  uint32_t seed,
                                  float    random_offset = 0.1f,
                                  bool     symmetric = true);

  /**
   * @brief Creates a matrix initialized to a uniform value.
   * @param  num_species Number of species.
   * @param  val         Uniform value.
   * @return             InteractionMatrix Resulting matrix.
   */
  static InteractionMatrix uniform(size_t num_species, float val);
};

// ============================================================================
//  Growth & Thinning Functions (Alphabetically Sorted)
// ============================================================================

/**
 * @brief Performs multi-species Strauss soft-core thinning on an input forest.
 *
 * @param  forest              Input candidate forest.
 * @param  repulsion_distances Interaction distance matrix between species.
 * @param  repulsion_strengths Repulsion strength matrix between species in [0,
 *                             1].
 * @param  target_count        Optional target number of trees (0 means keep all
 *                             accepted).
 * @param  seed                Random number generator seed.
 * @param  bbox                Bounding box {xmin, xmax, ymin, ymax}.
 * @return                     Forest Thinned forest.
 */
Forest thin_forest_soft_core(const Forest            &forest,
                             const InteractionMatrix &repulsion_distances,
                             const InteractionMatrix &repulsion_strengths,
                             size_t                   target_count = 0,
                             uint32_t                 seed = 0,
                             const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

} // namespace hmap
