/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

#pragma once
#include "highmap/boundary.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/math/profiles.hpp"

namespace hmap
{

class Graph;

// clang-format off
enum ElevationLongitudinalProfile : int
{
	ELP_UNCHANGED,
	ELP_FLAT,       ///< 
	ELP_DECREASING, ///< 
	ELP_INCREASING, ///< 
};
// clang-format on

/**
 * @brief Carves a continuous river network into a heightmap along a directed
 * flow Graph.
 *
 * A triangular riverbed with a given slope (talus) is generated from the
 * channel bottom using expand_talus, and merged with the original heightmap
 * using the network signed distance field (SDF), scaling channel width and
 * depth with flow accumulation.
 *
 * @param[in,out] z                        Target heightmap to carve in-place.
 * @param[in]     graph                    Input flow network Graph.
 * @param[in]     talus                    Slope of the triangular riverbed
 *                                         valley walls.
 * @param[in]     outlet_width             Maximum channel influence width at
 *                                         downstream outlet(s).
 * @param[in]     headwater_width          Minimum channel influence width at
 *                                         headwaters (default: 0).
 * @param[in]     width_exponent           Scaling exponent for channel width
 *                                         with flow accumulation (default:
 *                                         0.5).
 * @param[in]     depth                    Incision depth below graph elevation
 *                                         for river bottom (default: 0.05).
 * @param[in]     depth_exponent           Scaling exponent for depth scaling
 *                                         with flow accumulation (default:
 *                                         0.3).
 * @param[in]     seed                     Random seed for expand_talus
 *                                         (default: 0).
 * @param[in]     talus_noise_ratio        Roughness / perturbation ratio for
 *                                         valley expansion (default: 0).
 * @param[in]     radial_profile           Blending profile for SDF transition
 *                                         (default: RP_SMOOTHSTEP_UPPER).
 * @param[in]     radial_profile_parameter Shaping parameter for the radial
 *                                         profile (default: 2.0).
 * @param[out]    p_river_mask             Optional output mask in [0, 1] of
 *                                         carved river channel.
 * @param[in]     bbox                     Domain bounding box.
 */
void carve_river(
    Array        &z,
    const Graph  &graph,
    float         talus,
    float         outlet_width,
    float         headwater_width = 0.f,
    float         width_exponent = 0.5f,
    float         depth = 0.05f,
    float         depth_exponent = 0.3f,
    std::uint32_t seed = 0,
    float         talus_noise_ratio = 0.f,
    RadialProfile radial_profile = RadialProfile::RP_SMOOTHSTEP_UPPER,
    float         radial_profile_parameter = 2.f,
    Array        *p_river_mask = nullptr,
    glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Dig a path on a heightmap.
 *
 * This function modifies a heightmap array by "digging" a path into it based on
 * the provided path. The path is represented by a `Path` object, and the
 * function adjusts the height values in the `z` array to create the appearance
 * of a dug path. The width, border decay, and flattening radius parameters
 * control the characteristics of the dig. Optionally, the path can be forced to
 * have a monotonically decreasing elevation.
 *
 * **Example**
 * @include ex_dig_path.cpp
 *
 * **Result**
 * @image html ex_dig_path.png
 *
 * @param z                 Input array representing the heightmap to be
 *                          modified.
 * @param path              The path to be dug into the heightmap, with
 *                          coordinates with respect to a unit-square. The path
 *                          will be processed to create the dig effect.
 * @param width             Radius of the path width in pixels. This determines
 *                          how wide the dug path will be.
 * @param decay             Radius of the path border decay in pixels. This
 *                          controls how quickly the effect of the path fades
 *                          out towards the edges.
 * @param flattening_radius Radius used to flatten the elevation of the path,
 *                          creating a smooth transition. This is measured in
 *                          pixels.
 * @param force_downhill    If `true`, the path's elevation will be forced to
 *                          decrease monotonically, creating a downhill effect.
 * @param bbox              Bounding box specifying the region of the heightmap
 *                          to consider for the digging operation. It defines
 *                          the area where the path is applied.
 * @param depth             Optional depth parameter to specify the maximum
 *                          depth of the dig. If not provided, the default depth
 *                          of 0.f is used.
 */
void dig_path(Array    &z,
              Path     &path,
              int       width = 1,
              int       decay = 2,
              int       flattening_radius = 16,
              bool      force_downhill = false,
              glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f},
              float     depth = 0.f);

/**
 * @brief Blends a flatbed carve into an existing heightmap.
 *
 * Carves and blends a flatbed shape along a path into the provided array, using
 * a smooth falloff mask.
 *
 * @param z                        Heightmap to modify in place.
 * @param path                     Path defining the bed centerline.
 * @param bottom_extent            Radius of the flatbed region.
 * @param vmin                     Base height value.
 * @param depth                    Bed depth.
 * @param falloff_distance         Falloff width outside the bed.
 * @param outer_slope              Linear slope beyond the bed.
 * @param preserve_bedshape        Preserve relative bed shape under radial
 *                                 noise.
 * @param radial_profile           Radial profile type.
 * @param radial_profile_parameter Shape parameter for the profile.
 * @param p_falloff_mask           Optional output falloff mask.
 * @param p_noise_r                Optional radial noise.
 * @param bbox                     World-space bounding box.
 *
 * **Example**
 * @include ex_flatbed_carve.cpp
 *
 * **Result**
 * @image html ex_flatbed_carve.png
 */
void flatbed_carve(Array        &z,
                   const Path   &path,
                   float         bottom_extent = 32.f, // pixels
                   float         vmin = 0.1f,
                   float         depth = 0.05f,
                   float         falloff_distance = 128.f,
                   float         outer_slope = 0.1f,
                   bool          preserve_bedshape = true,
                   RadialProfile radial_profile = RadialProfile::RP_GAIN,
                   float         radial_profile_parameter = 2.f,
                   Array        *p_falloff_mask = nullptr,
                   const Array  *p_noise_r = nullptr,
                   glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Carves a trench along a given path into a heightmap.
 *
 * This function modifies the input heightmap by carving a trench defined by a
 * polyline path. The trench shape is controlled both longitudinally (along the
 * path) and radially (across the path). Optional noise and scaling can be
 * applied to introduce variations in width and depth.
 *
 * @param[in,out] z                          Heightmap to modify (2D array of
 *                                           elevations).
 * @param[in]     path                       Polyline defining the trench
 *                                           centerline.
 * @param[in]     width                      Base width of the trench.
 * @param[in]     enable_width_depth_scaling If true, modulates trench width
 *                                           based on local elevation
 *                                           difference.
 * @param[in]     radial_profile             Type of radial profile used to
 *                                           shape the trench cross-section.
 * @param[in]     radial_profile_parameter   Parameter controlling the radial
 *                                           profile shape.
 * @param[in]     longitudinal_profile       Controls how elevation evolves
 *                                           along the path (e.g., flat,
 *                                           increasing, decreasing).
 * @param[in]     elevation_shift            Global elevation offset applied
 *                                           along the path.
 * @param[in]     shift_ramp_start_ratio     Ratio of the path length used to
 *                                           smoothly ramp the elevation shift
 *                                           at the beginning.
 * @param[in]     shift_ramp_end_ratio       Ratio of the path length used to
 *                                           smoothly ramp the elevation shift
 *                                           at the end.
 * @param[in]     min_slope                  Minimum slope enforced between
 *                                           consecutive path points when
 *                                           monotonicity is required.
 * @param[in]     k_neighbors                Number of nearest neighbors used
 *                                           for interpolation.
 * @param[in]     p_noise_r                  Optional pointer to a noise array
 *                                           used to locally perturb trench
 *                                           width (can be nullptr).
 * @param[out]    p_bending_mask             Optional output mask representing
 *                                           blending weights (can be nullptr).
 * @param[in]     tile_size                  Tile size in pixels for spatial BVH
 *                                           partitioning (default: 32). A
 *                                           recommended rule of thumb is to
 *                                           choose tile_size on the order of
 *                                           the trench width in pixels:
 *                                           tile_size ~ clamp(width / dx, 16,
 *                                           128).
 * @param[in]     bbox                       Bounding box of the domain as
 *                                           (xmin, xmax, ymin, ymax).
 *
 * @warning The function assumes that @p elevation_shift is non-zero when
 * width-depth scaling is enabled, to avoid division by zero.
 *
 * **Example**
 * @include ex_trench.cpp
 *
 * **Result**
 * @image html ex_trench.png
 */
void trench(Array        &z,
            const Path   &path,
            float         width,
            bool          enable_width_depth_scaling = true,
            bool          enable_width_distance_scaling = true,
            bool          enable_width_curvature_scaling = false,
            float         curvature_radius_min = 1.f,
            float         curv_width_ratio_min = 0.5f,
            float         curv_width_ratio_max = 2.f,
            RadialProfile radial_profile = RadialProfile::RP_SMOOTHSTEP_UPPER,
            float         radial_profile_parameter = 2.f,
            ElevationLongitudinalProfile longitudinal_profile =
                ElevationLongitudinalProfile::ELP_DECREASING,
            float        elevation_shift = -0.05f,
            float        shift_ramp_start_ratio = 0.f,
            float        shift_ramp_end_ratio = 0.f,
            float        min_slope = 0.001f,
            size_t       k_neighbors = 8,
            const Array *p_noise_r = nullptr,
            Array       *p_bending_mask = nullptr,
            int          tile_size = 32,
            glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

void trench(Array                   &z,
            const std::vector<Path> &paths,
            float                    width,
            bool                     enable_width_depth_scaling = true,
            bool                     enable_width_distance_scaling = true,
            bool                     enable_width_curvature_scaling = false,
            float                    curvature_radius_min = 1.f,
            float                    curv_width_ratio_min = 0.5f,
            float                    curv_width_ratio_max = 2.f,
            RadialProfile radial_profile = RadialProfile::RP_SMOOTHSTEP_UPPER,
            float         radial_profile_parameter = 2.f,
            ElevationLongitudinalProfile longitudinal_profile =
                ElevationLongitudinalProfile::ELP_DECREASING,
            float        elevation_shift = -0.05f,
            float        shift_ramp_start_ratio = 0.f,
            float        shift_ramp_end_ratio = 0.f,
            float        min_slope = 0.001f,
            size_t       k_neighbors = 8,
            const Array *p_noise_r = nullptr,
            Array       *p_bending_mask = nullptr,
            int          tile_size = 32,
            glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Carves river channels for an entire flow network represented by a
 * Graph with per-node relative widths.
 *
 * Each node in the graph provides a relative width ratio in [0, 1] (scaled by
 * @p width to obtain the node's channel width). Along each edge between node
 * `u` and node `v`, width is linearly interpolated as a function of the
 * projection parameter `t`:
 *
 * `w(t) = width * ((1 - t) * r_u + t * r_v)`
 *
 * @param[in,out] z                          Target heightmap to carve in-place.
 * @param[in]     graph                      Input flow network Graph.
 * @param[in]     relative_widths            Relative width ratio in [0, 1] per
 *                                           graph node.
 * @param[in]     width                      Maximum base width.
 * @param[in]     enable_width_depth_scaling Modulate width according to depth.
 * @param[in]     radial_profile             Transverse cross-section profile.
 * @param[in]     radial_profile_parameter   Shaping parameter for radial
 *                                           profile.
 * @param[in]     longitudinal_profile       Downhill elevation monotonicity
 *                                           mode.
 * @param[in]     elevation_shift            Depth of trench incision.
 * @param[in]     min_slope                  Minimum downstream slope.
 * @param[in]     p_noise_r                  Optional noise array for radial
 *                                           perturbation.
 * @param[out]    p_bending_mask             Optional blending mask.
 * @param[in]     tile_size                  Spatial tile partition size.
 * @param[in]     bbox                       Domain bounding box.
 */
void trench(Array                    &z,
            const Graph              &graph,
            const std::vector<float> &relative_widths,
            float                     width,
            bool                      enable_width_depth_scaling = true,
            RadialProfile radial_profile = RadialProfile::RP_SMOOTHSTEP_UPPER,
            float         radial_profile_parameter = 2.f,
            ElevationLongitudinalProfile longitudinal_profile =
                ElevationLongitudinalProfile::ELP_DECREASING,
            float        elevation_shift = -0.05f,
            float        min_slope = 0.001f,
            const Array *p_noise_r = nullptr,
            Array       *p_bending_mask = nullptr,
            int          tile_size = 32,
            glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Carves river channels for an entire flow network represented by a
 * Graph, computing relative widths from distance to outlet.
 *
 * @param[in,out] z                          Target heightmap to carve in-place.
 * @param[in]     graph                      Input flow network Graph.
 * @param[in]     width                      Target river width at the
 *                                           outlet(s).
 * @param[in]     width_exponent             Scaling exponent for width with
 *                                           distance to outlet (default: 1.0).
 * @param[in]     enable_width_depth_scaling Modulate width according to depth.
 * @param[in]     radial_profile             Transverse cross-section profile.
 * @param[in]     radial_profile_parameter   Shaping parameter for radial
 *                                           profile.
 * @param[in]     longitudinal_profile       Downhill elevation monotonicity
 *                                           mode.
 * @param[in]     elevation_shift            Depth of trench incision.
 * @param[in]     min_slope                  Minimum downstream slope.
 * @param[in]     p_noise_r                  Optional noise array for radial
 *                                           perturbation.
 * @param[out]    p_bending_mask             Optional blending mask.
 * @param[in]     tile_size                  Spatial tile partition size.
 * @param[in]     bbox                       Domain bounding box.
 */
void trench(Array        &z,
            const Graph  &graph,
            float         width,
            float         width_exponent = 1.0f,
            bool          enable_width_depth_scaling = true,
            RadialProfile radial_profile = RadialProfile::RP_SMOOTHSTEP_UPPER,
            float         radial_profile_parameter = 2.f,
            ElevationLongitudinalProfile longitudinal_profile =
                ElevationLongitudinalProfile::ELP_DECREASING,
            float        elevation_shift = -0.05f,
            float        min_slope = 0.001f,
            const Array *p_noise_r = nullptr,
            Array       *p_bending_mask = nullptr,
            int          tile_size = 32,
            glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

} // namespace hmap