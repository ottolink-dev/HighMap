/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file poi.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for terrain-aware points of interest (POI) extraction and
 * procedural point generation.
 *
 * @copyright Copyright (c) 2026 Otto Link
 */
#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/geometry/cloud.hpp"

namespace hmap
{

/**
 * @enum PoiType
 * @brief Classification types for terrain-aware points of interest.
 */
enum class PoiType : int
{
  CUSTOM = 0,           ///< Custom or unclassified score field
  APEX,                 ///< Local summits / peaks
  SINK,                 ///< Depressions / pits / sinkholes
  SADDLE,               ///< Col / pass / saddle points
  RIDGE,                ///< Elongated high terrain
  VALLEY,               ///< Elongated low terrain
  CLIFF,                ///< High slope / drop-offs
  PLATEAU,              ///< Elevated flat areas
  FLAT_LOWLAND,         ///< Low flat areas / plains
  VIEWPOINT,            ///< High elevation + high prominence / openness
  SETTLEMENT_CANDIDATE, ///< Flat enough, accessible, good relief
};

/**
 * @struct PointOfInterest
 * @brief Represents a terrain-aware point of interest with physical
 * descriptors.
 */
struct PointOfInterest
{
  glm::vec2 position = {0.f, 0.f}; ///< Coordinates in domain or normalized units.
  float   score = 0.f;             ///< Suitability or feature score in [0, 1].
  float   elevation = 0.f;         ///< Terrain elevation at point.
  float   slope = 0.f;             ///< Slope magnitude / talus.
  float   aspect = 0.f;            ///< Aspect angle in radians [-π, π].
  PoiType type = PoiType::CUSTOM;  ///< Point type classification.
};

/**
 * @struct PoiGenerationOptions
 * @brief Configuration parameters for procedural point-of-interest generation.
 */
struct PoiGenerationOptions
{
  int scale_radius = 8;                  ///< Neighborhood radius for local metrics.
  float prominence = 0.05f;              ///< Minimum elevation prominence over neighborhood.
  float threshold = 0.5f;                ///< Score / suitability acceptance threshold in [0, 1].
  float min_distance = 0.05f;            ///< Minimum Euclidean separation distance between points.
  size_t max_points = 50;                ///< Maximum number of points to generate.
  float border_margin = 0.02f;           ///< Exclusion margin near bounding box edges.
  std::uint32_t seed = 1337;             ///< PRNG seed for deterministic generation and jitter.
  float jitter = 0.0f;                   ///< Stochastic displacement ratio in [0, 1].
  const Array *p_mask = nullptr;         ///< Optional exclusion / region-of-interest mask.
  glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f}; ///< Domain bounding box {xmin, xmax, ymin, ymax}.
};

// ============================================================================
//  POI Functions (Alphabetically Sorted)
// ============================================================================

/**
 * @brief Computes a normalized suitability/score field in [0, 1] for a given
 * terrain feature type.
 *
 * @param  heightmap Input terrain heightmap array.
 * @param  type      Target point of interest type.
 * @param  options   Generation options controlling neighborhood scale and
 *                   parameters.
 * @return           Array     Suitability field in [0, 1].
 */
Array compute_poi_suitability_field(const Array                &heightmap,
                                    PoiType                     type,
                                    const PoiGenerationOptions &options);

/**
 * @brief Generates a Cloud of points of interest from a heightmap.
 *
 * @param  heightmap Input terrain heightmap array.
 * @param  type      Target point of interest type.
 * @param  options   Generation options.
 * @return           Cloud     Point cloud carrying the score as the point
 *                   value.
 */
Cloud generate_points_of_interest_cloud(const Array                &heightmap,
                                        PoiType                     type,
                                        const PoiGenerationOptions &options);

/**
 * @brief Generates a list of PointOfInterest objects from a heightmap.
 *
 * @param  heightmap Input terrain heightmap array.
 * @param  type      Target point of interest type.
 * @param  options   Generation options.
 * @return           std::vector<PointOfInterest> List of extracted points of
 *                   interest.
 */
std::vector<PointOfInterest> generate_points_of_interest(
    const Array                &heightmap,
    PoiType                     type,
    const PoiGenerationOptions &options);

/**
 * @brief Generates 2D positions of points of interest from a heightmap.
 *
 * @param  heightmap Input terrain heightmap array.
 * @param  type      Target point of interest type.
 * @param  options   Generation options.
 * @return           std::vector<glm::vec2> List of 2D positions.
 */
std::vector<glm::vec2> generate_points_of_interest_positions(
    const Array                &heightmap,
    PoiType                     type,
    const PoiGenerationOptions &options);

/**
 * @brief Generates a list of PointOfInterest objects from a precomputed score
 * field.
 *
 * @param  score_field Precomputed suitability / score field in [0, 1].
 * @param  heightmap   Input terrain heightmap used to sample elevation, slope,
 *                     and aspect.
 * @param  options     Generation options.
 * @param  type        Optional classification type to assign to the points.
 * @return             std::vector<PointOfInterest> List of extracted points of
 *                     interest.
 */
std::vector<PointOfInterest> generate_points_from_score_field(
    const Array                &score_field,
    const Array                &heightmap,
    const PoiGenerationOptions &options,
    PoiType                     type = PoiType::CUSTOM);

} // namespace hmap
