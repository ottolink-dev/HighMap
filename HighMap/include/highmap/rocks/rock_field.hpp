/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file rock_field.hpp
 * @copyright Copyright (c) 2026 Otto Link.
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/rocks/rock.hpp"
#include "highmap/rocks/rock_distribution.hpp"
#include "highmap/rocks/rock_simulation.hpp"
#include "highmap/scatter/scatter_field.hpp"

namespace hmap
{

/**
 * @class RockField
 * @brief Container class representing a field of rocks/boulders (scree, talus,
 * or boulder field).
 */
class RockField : public ScatterField
{
public:
  // ==========================================================================
  //  Constructors
  // ==========================================================================

  /**
   * @brief Default constructor initializing an empty rock field.
   */
  RockField() = default;

  /**
   * @brief Constructs a RockField from a vector of rocks.
   * @param rocks Vector of Rock instances.
   */
  RockField(const std::vector<Rock> &rocks);

  /**
   * @brief Move-constructs a RockField from a vector of rocks.
   * @param rocks Rvalue vector of Rock instances.
   */
  RockField(std::vector<Rock> &&rocks) noexcept;

  /**
   * @brief Constructs a RockField from a base ScatterField.
   * @param field Base scatter field.
   */
  RockField(const ScatterField &field);

  /**
   * @brief Move-constructs a RockField from a base ScatterField.
   * @param field Rvalue base scatter field.
   */
  RockField(ScatterField &&field) noexcept;

  /**
   * @brief Constructs a RockField from a Cloud of 2D points.
   *
   * @param cloud          Input point cloud.
   * @param class_id       Class identifier assigned to all imported rocks.
   * @param default_radius Default radius assigned if point value is 0.
   */
  RockField(const Cloud &cloud,
            uint32_t     class_id = 0,
            float        default_radius = HMAP_DEFAULT_ROCK_RADIUS);

  // ==========================================================================
  //  Rock-Specific Operations
  // ==========================================================================

  /**
   * @brief Applies a heavy-tailed power-law (Pareto) radius distribution to all
   * rocks.
   *
   * Radii are sampled according to the probability density $P(R \ge r) \propto
   * r^{-\alpha}$, bounded within $[r_{\min}, r_{\max}]$.
   *
   * @param alpha Power-law exponent (e.g. 2.0).
   * @param r_min Minimum rock radius.
   * @param r_max Maximum rock radius.
   * @param seed  Random seed for sampling.
   */
  void apply_power_law_distribution(float    alpha = 2.0f,
                                    float    r_min = 1e-3f,
                                    float    r_max = 5e-2f,
                                    uint32_t seed = 0);

  /**
   * @brief Applies a power-law radius distribution according to a
   * RockDistribution definition.
   *
   * @param dist RockDistribution parameters.
   * @param seed Random seed for sampling.
   */
  void apply_power_law_distribution(const RockDistribution &dist,
                                    uint32_t                seed = 0);

  /**
   * @brief Simulates gravitational slope sorting (scree/talus dynamics).
   *
   * On natural scree and talus slopes, larger boulders carry more momentum and
   * roll further downhill toward gentler base slopes, while smaller stones
   * remain caught on steeper sections.
   *
   * @param slope            2D slope gradient array (e.g., from
   *                         `gradient_magnitude` or `slope`).
   * @param sorting_strength Strength factor in [0, 1] governing the correlation
   *                         between slope and rock size.
   * @param bbox             Domain bounding box {xmin, xmax, ymin, ymax}.
   */
  void apply_slope_sorting(const Array     &slope,
                           float            sorting_strength = 1.0f,
                           const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

  /**
   * @brief Packs smaller interstitial pebbles/stones in the voids between
   * existing larger rocks without displacing the existing rocks.
   *
   * @param count    Target number of small rocks to pack.
   * @param density  Optional 2D spatial probability density array for candidate
   *                 sampling.
   * @param r_min    Minimum interstitial rock radius.
   * @param r_max    Maximum interstitial rock radius.
   * @param class_id Class identifier for the packed stones.
   * @param seed     Random seed for candidate sampling.
   * @param bbox     Domain bounding box {xmin, xmax, ymin, ymax}.
   */
  void pack_interstitial_rocks(size_t           count,
                               const Array     &density = {},
                               float            r_min = 5e-4f,
                               float            r_max = 2e-3f,
                               uint32_t         class_id = 0,
                               uint32_t         seed = 0,
                               const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

  /**
   * @brief Simulates physical movement (sliding, rolling, bouncing, and
   * inter-rock collision avoidance) for this RockField across terrain on GPU.
   *
   * @param elevation      Terrain elevation heightmap.
   * @param p_friction_map Optional soil friction map.
   * @param options        Simulation options.
   */
  void simulate_physics(const Array                 &elevation,
                        const Array                 *p_friction_map = nullptr,
                        const RockSimulationOptions &options = {});

  /**
   * @brief Returns a multi-line formatted summary string of the rock field.
   * @return std::string Pretty-printed summary.
   */
  std::string to_string() const override;
};

} // namespace hmap
