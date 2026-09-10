/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file geo.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2023
 */
#pragma once

#include "highmap/array.hpp"
#include "highmap/functions.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/math/profiles.hpp"

namespace hmap
{

/**
 * @brief Return a caldera-shaped heightmap.
 *
 * @param  shape         Array shape.
 * @param  radius        Crater radius at the ridge.
 * @param  sigma_inner   Inner half-width.
 * @param  sigma_outer   Outer half-width.
 * @param  z_bottom      Bottom elevation (ridge is at elevation `1`).
 * @param  p_noise       Displacement noise.
 * @param  noise_r_amp   Radial noise amplitude.
 * @param  noise_ratio_z Vertical noise relative scale (in [0, 1]).
 * @param  center        Primitive reference center.
 * @param  bbox          Domain bounding box.
 * @return               Array Resulting array.
 *
 * **Example**
 * @include ex_caldera.cpp
 *
 * **Result**
 * @image html ex_caldera.png
 */
Array caldera(glm::ivec2   shape,
              float        radius,
              float        sigma_inner,
              float        sigma_outer,
              float        z_bottom,
              const Array *p_noise,
              float        noise_r_amp,
              float        noise_ratio_z,
              glm::vec2    center = {0.5f, 0.5f},
              glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f});

Array caldera(glm::ivec2 shape,
              float      radius,
              float      sigma_inner,
              float      sigma_outer,
              float      z_bottom,
              glm::vec2  center = {0.5f, 0.5f},
              glm::vec4  bbox = {0.f, 1.f, 0.f, 1.f}); ///< @overload

/**
 * @brief Generate a procedural impact crater heightmap.
 *
 * Creates a crater-shaped terrain with configurable inner cavity, outer lip,
 * asymmetry, terraces, and optional central peak. Optional masks for the crater
 * region and inner crater region can also be generated.
 *
 * @param  shape               Output array shape.
 * @param  radius              Crater radius.
 * @param  center              Crater center position.
 * @param  angle               Main asymmetry orientation angle in degrees.
 * @param  inner_depth         Depth of the crater cavity.
 * @param  inner_exp           Exponent controlling the inner profile shape.
 * @param  lip_height          Height of the crater rim.
 * @param  lip_extent          Extent of the outer rim.
 * @param  lip_exp             Exponent controlling rim falloff.
 * @param  asym_ratio          Asymmetry factor applied along the main
 *                             direction.
 * @param  central_peak_height Height of the central peak.
 * @param  central_peak_extent Radius of the central peak region.
 * @param  n_terraces          Number of inner terraces.
 * @param  terrace_extent      Terrace spacing factor.
 * @param  terrace_exp         Terrace profile exponent.
 * @param  terrace_persistence Terrace amplitude decay factor.
 * @param  p_noise_r           Optional radial noise distortion.
 * @param  bbox                Bounding box coordinates.
 * @param  p_crater_mask       Optional output crater mask.
 * @param  p_inner_crater_mask Optional output inner crater mask.
 * @return                     Generated crater heightmap.
 *
 * **Example**
 * @include ex_crater.cpp
 *
 * **Result**
 * @image html ex_crater.png
 */
Array crater(glm::ivec2   shape,
             float        radius,
             glm::vec2    center = {0.5f, 0.5f},
             float        angle = 30.f, // degs,
             float        inner_depth = 0.1f,
             float        inner_exp = 3.f,
             float        lip_height = 0.02f,
             float        lip_extent = 0.25f,
             float        lip_exp = 2.f,
             float        asym_ratio = 0.8f,
             float        central_peak_height = 0.01f,
             float        central_peak_extent = 0.4f,
             int          n_terraces = 0,
             float        terrace_extent = 0.1f,
             float        terrace_exp = 2.f,
             float        terrace_persistence = 0.5f,
             const Array *p_noise_r = nullptr,
             glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f},
             Array       *p_crater_mask = nullptr,
             Array       *p_inner_crater_mask = nullptr);

/**
 * @brief Generates a 2D island mask by perturbing a radial boundary with noise.
 *
 * Creates a binary mask of a single connected island using radial distance from
 * a center point and adding an FBM-based angular displacement.
 *
 * @param  shape        Output array size (width, height).
 * @param  radius       Base island radius.
 * @param  seed         Noise seed.
 * @param  displacement Amplitude of boundary perturbation.
 * @param  noise_type   Type of underlying noise function.
 * @param  kw           Angular frequency for sampling the noise.
 * @param  octaves      Number of FBM octaves.
 * @param  weight       Base amplitude of the FBM.
 * @param  persistence  Amplitude falloff per octave.
 * @param  lacunarity   Frequency multiplier per octave.
 * @param  center       Center of the island in normalized coordinates.
 * @param  bbox         Bounding box for coordinate remapping.
 *
 * @return              Binary mask where 1.f indicates land and 0.f indicates
 *                      water.
 *
 * **Example**
 * @include ex_island.cpp
 *
 * **Result**
 * @image html ex_island.png
 */
Array island_land_mask(glm::ivec2       shape,
                       float            radius,
                       std::uint32_t    seed,
                       float            displacement = 0.2f,
                       NoiseType        noise_type = NoiseType::SIMPLEX2S,
                       float            kw = 4.f,
                       int              octaves = 8,
                       float            weight = 0.f,
                       float            persistence = 0.5f,
                       float            lacunarity = 2.f,
                       const glm::vec2 &center = {0.5f, 0.5f},
                       const glm::vec4 &bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a land mask of an island chain scattered along a path.
 *
 * Spawns `island_count` islands along the input path. Each island is an
 * fbm-perturbed radial blob (same construction as `island_land_mask`). Island
 * positions follow the path arc length, optionally jittered along and across
 * the path by `scatter`. Island sizes follow a signed falloff:
 * `size_falloff = 0` gives uniform sizes, `+1` shrinks the islands toward the
 * path end (hotspot-chain look, largest island at the head), `-1` shrinks them
 * toward the start; `size_jitter` adds per-island random variation.
 *
 * @param  shape         Output array shape.
 * @param  path          Input path (at least two points), in normalized domain
 *                       coordinates.
 * @param  seed          Random seed number.
 * @param  island_count  Number of islands spawned along the path.
 * @param  island_radius Base island radius in normalized domain units.
 * @param  size_falloff  Signed size falloff along the path, in [-1, 1].
 * @param  size_jitter   Per-island random radius variation ratio.
 * @param  scatter       Random island displacement along and across the path.
 * @param  displacement  Amplitude of the fbm boundary perturbation.
 * @param  noise_type    Noise type used for the boundary perturbation.
 * @param  kw            Noise wavenumber.
 * @param  octaves       Number of fbm octaves.
 * @param  weight        Octave weighting.
 * @param  persistence   Octave persistence.
 * @param  lacunarity    Defines the wavenumber ratio between octaves.
 * @param  bbox          Domain bounding box.
 * @return               Binary land mask (0 or 1).
 *
 * **Example**
 * @include ex_island_chain.cpp
 *
 * **Result**
 * @image html ex_island_chain.png
 */
Array island_chain_land_mask(glm::ivec2    shape,
                             const Path   &path,
                             std::uint32_t seed,
                             int           island_count = 5,
                             float         island_radius = 0.1f,
                             float         size_falloff = 0.5f,
                             float         size_jitter = 0.3f,
                             float         scatter = 0.f,
                             float         displacement = 0.2f,
                             NoiseType     noise_type = NoiseType::SIMPLEX2S,
                             float         kw = 4.f,
                             int           octaves = 8,
                             float         weight = 0.f,
                             float         persistence = 0.5f,
                             float         lacunarity = 2.f,
                             glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generate a rift-shaped primitive.
 *
 * Creates a longitudinal depression centered around a rotated axis. The rift
 * shape is defined by an outer radial profile, an optional bottom profile, and
 * configurable inner/outer slopes. Noise fields can be used to locally perturb
 * the radius and axial offset of the structure.
 *
 * @param  shape                Output array shape.
 * @param  angle                Rift orientation angle in degrees.
 * @param  radius               Rift half-width in normalized domain units.
 * @param  axial_slope          Longitudinal slope applied along the rift axis.
 * @param  depth                Rift depth amplitude.
 * @param  scale_with_depth     If true, scale radial features proportionally to
 *                              the depth value.
 * @param  profile              Radial profile used for the outer rift shape.
 * @param  profile_param        Additional parameter for the outer profile.
 * @param  bottom_extent        Relative extent of the flat or inner bottom
 *                              region inside the rift.
 * @param  bottom_depth         Additional depth contribution applied to the
 *                              bottom region.
 * @param  bottom_profile       Radial profile used for the bottom region.
 * @param  bottom_profile_param Additional parameter for the bottom profile.
 * @param  outer_slope          Slope factor controlling the outer falloff.
 * @param  p_noise_r            Optional noise array used to locally modulate
 *                              the rift radius.
 * @param  p_noise_offset       Optional noise array used to locally offset the
 *                              rift axis position.
 * @param  center               Rift center in normalized coordinates.
 * @param  bbox                 Domain bounding box.
 * @return                      Generated rift array.
 *
 * **Example**
 * @include ex_rift.cpp
 *
 * **Result**
 * @image html ex_rift.png
 */
Array rift(glm::ivec2    shape,
           float         angle,
           float         radius = 0.1f,
           float         axial_slope = 0.2f,
           float         depth = 0.5f,
           bool          scale_with_depth = true,
           RadialProfile profile = RadialProfile::RP_SMOOTHSTEP,
           float         profile_param = 0.f,
           float         bottom_extent = 0.1,
           float         bottom_depth = 0.05f,
           RadialProfile bottom_profile = RadialProfile::RP_SQRT,
           float         bottom_profile_param = 0.f,
           bool          bottom_force_minimum_depth = true,
           float         outer_slope = 0.5f,
           const Array  *p_noise_r = nullptr,
           const Array  *p_noise_offset = nullptr,
           glm::vec2     center = {0.5f, 0.5f},
           glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f},
           Array        *p_rift_mask = nullptr,
           Array        *p_bottom_mask = nullptr);

/**
 * @brief Generate a procedural natural valley head / valley incision
 * deformation heightmap.
 *
 * Creates a valley head incision that smoothly emerges at the origin and
 * progressively deepens and widens downstream following exponential saturation
 * laws and smooth lateral falloff.
 *
 * @param  shape             Output array shape.
 * @param  angle             Valley orientation angle in degrees.
 * @param  max_depth         Maximum incision depth at maturity.
 * @param  min_width         Valley width at the head.
 * @param  max_width         Mature valley width.
 * @param  depth_length      Distance scale over which depth develops.
 * @param  width_length      Distance scale over which width develops.
 * @param  development_power Exponent controlling how quickly the incision
 *                           emerges.
 * @param  profile_power     Cross-section exponent (1 = sharp V, > 1 = rounded
 *                           U-shape).
 * @param  p_noise_offset    Optional cross-section centerline offset noise.
 * @param  p_noise_r         Optional width / radial noise.
 * @param  center            Valley head position.
 * @param  bbox              Bounding box coordinates.
 * @param  p_mask            Optional output mask for incision intensity.
 * @return                   Generated valley head deformation array.
 *
 * **Example**
 * @include ex_valley_head.cpp
 *
 * **Result**
 * @image html ex_valley_head.png
 */
Array valley_head(glm::ivec2   shape,
                  float        angle = 0.f,
                  float        max_depth = 0.2f,
                  float        min_width = 0.02f,
                  float        max_width = 0.15f,
                  float        depth_length = 0.5f,
                  float        width_length = 0.5f,
                  float        development_power = 2.f,
                  float        profile_power = 1.3f,
                  const Array *p_noise_offset = nullptr,
                  const Array *p_noise_r = nullptr,
                  glm::vec2    center = {0.5f, 0.5f},
                  glm::vec4    bbox = {0.f, 1.f, 0.f, 1.f},
                  Array       *p_mask = nullptr);

} // namespace hmap

namespace hmap::gpu
{

/**
 * @brief Generates a synthetic "badlands" terrain heightmap.
 *
 * This function procedurally creates a 2D elevation field resembling badlands
 * (highly eroded terrain with sharp ridges and valleys). It combines fractal
 * noise (FBM) with a Voronoi-based primitive, displaced along a specified
 * orientation.
 *
 * @param  shape          Output array shape (resolution in x and y).
 * @param  kw             Frequency vector controlling the scale of the
 *                        features.
 * @param  seed           Random seed used for noise and Voronoi generation.
 * @param  rugosity       Controls roughness of the fractal noise (higher = more
 *                        irregular).
 * @param  angle          Orientation angle (degrees) of terrain displacements.
 * @param  k_smoothing    Voronoi smoothing parameter (controls ridge
 *                        sharpness).
 * @param  base_noise_amp Amplitude of the base displacement noise.
 * @param  p_noise_x      Optional pointer to external displacement noise field
 *                        (X-axis).
 * @param  p_noise_y      Optional pointer to external displacement noise field
 *                        (Y-axis).
 * @param  bbox           Bounding box of the generation domain in normalized
 *                        coordinates.
 *
 * @return                Array containing the generated badlands heightmap.
 *
 * **Example**
 * @include ex_badlands.cpp
 *
 * **Result**
 * @image html ex_badlands.png
 */
Array badlands(glm::ivec2    shape,
               glm::vec2     kw,
               std::uint32_t seed,
               int           octaves = 8,
               float         rugosity = 0.2f,
               float         angle = 30.f,
               float         k_smoothing = 0.1f,
               float         base_noise_amp = 0.2f,
               const Array  *p_noise_x = nullptr,
               const Array  *p_noise_y = nullptr,
               glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a synthetic procedural terrain resembling basaltic
 * landforms.
 *
 * This function creates a multi-scale procedural field combining large, medium,
 * and small-scale Voronoi-based patterns, noise warping, and optional
 * flattening, simulating the morphology of fractured basalt or volcanic
 * terrains. The terrain is constructed using a combination of Voronoi diagrams
 * (via `voronoi_fbm`) and fractal noise (`noise_fbm`), layered with
 * frequency-domain manipulations and amplitude/gain controls at each scale.
 *
 * The final output is a heightmap represented as an `Array`, normalized and
 * composed of:
 * - Large-scale cellular patterns with smoothed Voronoi edge distances.
 * - Medium and small-scale structures introducing finer surface variation.
 * - Optional rugosity (fine detail) and flattening to simulate erosion or flow
 * effects.
 *
 * @param  shape                   Output resolution (width x height) of the
 *                                 field.
 * @param  kw                      Base wave numbers (frequency) for the terrain
 *                                 features.
 * @param  seed                    Initial seed used for deterministic random
 *                                 generation.
 * @param  warp_kw                 Frequency of the warping noise that displaces
 *                                 Voronoi positions.
 * @param  large_scale_warp_amp    Amplitude of displacement for large-scale
 *                                 Voronoi warping.
 * @param  large_scale_gain        Gain adjustment applied to the large-scale
 *                                 features.
 * @param  large_scale_amp         Final amplitude of the large-scale height
 *                                 contribution.
 * @param  medium_scale_kw_ratio   Scaling factor for the frequency of the
 *                                 medium-scale patterns.
 * @param  medium_scale_warp_amp   Amplitude of warping for the medium-scale
 *                                 displacement.
 * @param  medium_scale_gain       Gain control for medium-scale modulation.
 * @param  medium_scale_amp        Amplitude of the medium-scale heightmap.
 * @param  small_scale_kw_ratio    Frequency ratio for small-scale details.
 * @param  small_scale_amp         Amplitude of small-scale pattern
 *                                 contribution.
 * @param  small_scale_overlay_amp Additional overlay strength for repeating the
 *                                 small-scale pattern.
 * @param  rugosity_kw_ratio       Frequency ratio for high-frequency noise
 *                                 applied as fine roughness.
 * @param  rugosity_amp            Strength of the rugosity (high-frequency
 *                                 modulation).
 * @param  flatten_activate        Enables or disables the final flattening
 *                                 operation.
 * @param  flatten_kw_ratio        Frequency scaling of the flattening noise
 *                                 field.
 * @param  flatten_amp             Amplitude control of the flattening
 *                                 operation.
 * @param  p_noise_x               Optional pointer to a noise field used to
 *                                 displace grid coordinates in X.
 * @param  p_noise_y               Optional pointer to a noise field used to
 *                                 displace grid coordinates in Y.
 * @param  bbox                    The 2D bounding box ({xmin, xmax, ymin,
 *                                 ymax}) over which the terrain is generated.
 *
 * @return                         A procedurally generated `Array` representing
 *                                 the synthetic basalt-like terrain field.
 *
 * @note
 * - This function relies on OpenCL-based kernels via the `gpu::` namespace.
 * - The returned field is normalized in amplitude but may require rescaling to
 * match specific physical units.
 * - Adjusting `seed`, `warp_kw`, and the gain/amplitude values can produce a
 * wide variety of terrain features.
 *
 * **Example**
 * @include ex_basalt_field.cpp
 *
 * **Result**
 * @image html ex_basalt_field.png
 */
Array basalt_field(glm::ivec2    shape,
                   glm::vec2     kw,
                   std::uint32_t seed,
                   float         warp_kw = 4.f,
                   float         large_scale_warp_amp = 0.2f,
                   float         large_scale_gain = 6.f,
                   float         large_scale_amp = 0.2f,
                   float         medium_scale_kw_ratio = 3.f,
                   float         medium_scale_warp_amp = 1.f,
                   float         medium_scale_gain = 7.f,
                   float         medium_scale_amp = 0.08f,
                   float         small_scale_kw_ratio = 10.f,
                   float         small_scale_amp = 0.1f,
                   float         small_scale_overlay_amp = 0.002f,
                   float         rugosity_kw_ratio = 1.f,
                   float         rugosity_amp = 1.f,
                   bool          flatten_activate = true,
                   float         flatten_kw_ratio = 1.f,
                   float         flatten_amp = 0.f,
                   const Array  *p_noise_x = nullptr,
                   const Array  *p_noise_y = nullptr,
                   glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates an island heightmap from a land mask using radial profiles,
 * slope shaping, optional noise, and water attenuation.
 *
 * Applies elevation shaping inside the land mask using distance-based profiles,
 * slope controls, optional radial and slope noise, and water-depth modeling.
 *
 * @param  land_mask              Binary mask defining the island shape.
 * @param  p_noise_r              Optional radial noise field (same size as
 *                                mask).
 * @param  apex_elevation         Maximum elevation at island center.
 * @param  filter_distance        Enable smoothing of the distance transform.
 * @param  filter_ir              Radius of the distance filter.
 * @param  slope_min              Minimum slope.
 * @param  slope_max              Maximum slope.
 * @param  slope_start            Start of slope ramp (0–1).
 * @param  slope_end              End of slope ramp (0–1).
 * @param  slope_noise_intensity  Intensity of slope noise modulation.
 * @param  k_smooth               Smoothing factor for the final profile.
 * @param  radial_noise_intensity Intensity of radial displacement noise.
 * @param  radial_profile_gain    Exponent for radial falloff.
 * @param  water_decay            Transition sharpness to water level.
 * @param  water_depth            Water depth below the shoreline.
 * @param  lee_angle              Angle for lee-side (wind shadow) erosion.
 * @param  lee_amp                Amplitude of lee-side erosion.
 * @param  uplift_amp             Amplitude of uplift (orographic) effect.
 * @param  p_water_depth          Optional output: per-pixel water depth.
 * @param  p_inland_mask          Optional output: mask of inland pixels.
 *
 * @return                        Heightmap array representing the generated
 *                                island.
 *
 * **Example**
 * @include ex_island.cpp
 *
 * **Result**
 * @image html ex_island.png
 */
Array island(const Array &land_mask,
             const Array *p_noise_r = nullptr,
             float        apex_elevation = 1.f,
             bool         filter_distance = true,
             int          filter_ir = 32,
             float        slope_min = 0.1f,
             float        slope_max = 8.f,
             float        slope_start = 0.5f,
             float        slope_end = 1.f,
             float        slope_noise_intensity = 0.5f,
             float        k_smooth = 0.05f,
             float        radial_noise_intensity = 0.3f,
             float        radial_profile_gain = 2.f,
             float        water_decay = 0.05f,
             float        water_depth = 0.3f,
             float        lee_angle = 30.f,
             float        lee_amp = 0.f,
             float        uplift_amp = 0.f,
             Array       *p_water_depth = nullptr,
             Array       *p_inland_mask = nullptr);

/**
 * @brief Generates an island heightmap from a land mask using internally
 * generated FBM noise for radial perturbation and slope modulation.
 *
 * Same as the other overload but creates a noise field procedurally from a
 * seed, including frequency, octaves, orientation, and smoothing controls.
 *
 * @param  land_mask              Binary mask defining the island shape.
 * @param  seed                   Noise seed.
 * @param  noise_amp              Amplitude of generated noise.
 * @param  noise_kw               Noise frequencies (x, y).
 * @param  noise_octaves          Number of FBM octaves.
 * @param  noise_rugosity         Persistence-like roughness control.
 * @param  noise_angle            Rotation of the noise field (radians).
 * @param  noise_k_smoothing      Smoothing applied to the noise field.
 * @param  apex_elevation         Maximum elevation at island center.
 * @param  filter_distance        Enable smoothing of the distance transform.
 * @param  filter_ir              Radius of distance filter.
 * @param  slope_min              Minimum slope.
 * @param  slope_max              Maximum slope.
 * @param  slope_start            Start of slope ramp (0–1).
 * @param  slope_end              End of slope ramp (0–1).
 * @param  slope_noise_intensity  Intensity of slope noise modulation.
 * @param  k_smooth               Smoothing factor for final profile.
 * @param  radial_noise_intensity Intensity of radial displacement noise.
 * @param  radial_profile_gain    Exponent for radial falloff.
 * @param  water_decay            Transition sharpness to water level.
 * @param  water_depth            Water depth below the shoreline.
 * @param  lee_angle              Angle for lee-side erosion.
 * @param  lee_amp                Amplitude for lee-side erosion.
 * @param  uplift_amp             Amplitude of uplift (orographic) effect.
 * @param  p_water_depth          Optional output: per-pixel water depth.
 * @param  p_inland_mask          Optional output: mask of inland pixels.
 *
 * @return                        Heightmap array representing the generated
 *                                island.
 *
 * **Example**
 * @include ex_island.cpp
 *
 * **Result**
 * @image html ex_island.png
 */
Array island(const Array  &land_mask,
             std::uint32_t seed,
             float         noise_amp = 0.07f,
             glm::vec2     noise_kw = {4.f, 4.f},
             int           noise_octaves = 8,
             float         noise_rugosity = 0.7f,
             float         noise_angle = 45.f,
             float         noise_k_smoothing = 0.05f,
             float         apex_elevation = 1.f,
             bool          filter_distance = true,
             int           filter_ir = 32,
             float         slope_min = 0.1f,
             float         slope_max = 8.f,
             float         slope_start = 0.5f,
             float         slope_end = 1.f,
             float         slope_noise_intensity = 0.5f,
             float         k_smooth = 0.05f,
             float         radial_noise_intensity = 0.3f,
             float         radial_profile_gain = 2.f,
             float         water_decay = 0.05f,
             float         water_depth = 0.3f,
             float         lee_angle = 30.f,
             float         lee_amp = 0.f,
             float         uplift_amp = 0.f,
             Array        *p_water_depth = nullptr,
             Array        *p_inland_mask = nullptr);

/**
 * @brief Generates a procedural "mountain cone" heightmap using fractal noise
 * and Voronoi patterns.
 *
 * This function synthesizes a terrain-like array shaped as a cone with
 * noise-based ridges and surface roughness. The generation combines a
 * cone-shaped envelope (`cone_sigmoid`), a fractal base noise
 * (`gpu::noise_fbm`), and a Voronoi ridge structure (`voronoi_fbm`), producing
 * a mountain with controlled smoothness, angular distortion, and noise-based
 * displacements.
 *
 * @param  shape          The output array dimensions (width, height).
 * @param  seed           The random seed for noise generation.
 * @param  scale          Global scaling factor for the mountain size.
 * @param  octaves        Number of fractal noise octaves for both FBM and
 *                        Voronoi.
 * @param  peak_kw        Controls the base frequency (inverse wavelength) of
 *                        the peak noise.
 * @param  rugosity       Controls the amplitude decay across octaves (higher =
 *                        rougher surface).
 * @param  angle          Direction (in degrees) for the displacement noise
 *                        distortion.
 * @param  k_smoothing    Smoothing factor applied in Voronoi blending.
 * @param  gamma          Gamma correction exponent applied at the end.
 * @param  cone_alpha     Controls the cone envelope steepness (higher = sharper
 *                        peak).
 * @param  ridge_amp      Controls the amplitude of ridge enhancement in the
 *                        Voronoi term.
 * @param  base_noise_amp Amplitude multiplier for the displacement noise.
 * @param  center         The 2D position (in normalized coordinates) of the
 *                        mountain cone center.
 * @param  p_noise_x      Optional pointer to an external X displacement noise
 *                        field (can be nullptr).
 * @param  p_noise_y      Optional pointer to an external Y displacement noise
 *                        field (can be nullptr).
 * @param  bbox           The bounding box (min_x, min_y, max_x, max_y) defining
 *                        the spatial extent of the map.
 *
 * @return                An Array representing the final terrain heightmap in
 *                        normalized [0, 1] range.
 *
 * **Example**
 * @include ex_mountain_cone.cpp
 *
 * **Result**
 * @image html ex_mountain_cone.png
 */
Array mountain_cone(glm::ivec2    shape,
                    std::uint32_t seed,
                    float         scale = 1.f,
                    int           octaves = 8,
                    float         peak_kw = 4.f,
                    float         rugosity = 0.f,
                    float         angle = 45.f,
                    float         k_smoothing = 0.f,
                    float         gamma = 0.5f,
                    float         cone_alpha = 1.f,
                    float         ridge_amp = 0.4f,
                    float         base_noise_amp = 0.05f,
                    glm::vec2     center = {0.5f, 0.5f},
                    const Array  *p_noise_x = nullptr,
                    const Array  *p_noise_y = nullptr,
                    glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a synthetic mountain-like inselberg (isolated hill)
 * heightmap.
 *
 * This function procedurally creates a 2D elevation field resembling an
 * inselberg, using a combination of fractal noise (FBM), Gaussian pulses, and
 * Voronoi-based structures. It allows control over scale, shape, orientation,
 * rugosity, and optional geological-like effects such as bulk uplift and
 * deposition smoothing.
 *
 * @param  shape          Output array shape (resolution in x and y).
 * @param  seed           Random seed used for noise and Voronoi generation.
 * @param  scale          Global scaling factor for the inselberg width and
 *                        structure.
 * @param  rugosity       Controls roughness of the fractal noise (higher = more
 *                        irregular).
 * @param  angle          Orientation angle (degrees) for directional
 *                        displacements.
 * @param  gamma          Gamma correction factor applied to the final
 *                        heightmap.
 * @param  round_shape    If true, generates a symmetric round shape;
 *                        if false, adds directional displacement.
 * @param  add_deposition If true, applies a smoothing fill step simulating
 *                        sediment deposition.
 * @param  bulk_amp       Amplitude of bulk uplift applied to the base pulse (0
 *                        = none, >0 = raises and normalizes the feature).
 * @param  base_noise_amp Amplitude of the base displacement noise.
 * @param  k_smoothing    Voronoi smoothing parameter (controls ridge
 *                        sharpness).
 * @param  center         Center of the inselberg in normalized coordinates.
 * @param  p_noise_x      Optional pointer to external displacement noise field
 *                        (X-axis).
 * @param  p_noise_y      Optional pointer to external displacement noise field
 *                        (Y-axis).
 * @param  bbox           Bounding box of the generation domain in normalized
 *                        coordinates.
 *
 * @return                Array containing the generated inselberg heightmap.
 *
 * **Example**
 * @include ex_mountain_inselberg.cpp
 *
 **Result**
 * @image html ex_mountain_inselberg.png
 */
Array mountain_inselberg(glm::ivec2    shape,
                         std::uint32_t seed,
                         float         scale = 1.f,
                         int           octaves = 8,
                         float         rugosity = 0.2f,
                         float         angle = 45.f,
                         float         gamma = 1.1f,
                         bool          round_shape = false,
                         bool          add_deposition = true,
                         float         bulk_amp = 0.2f,
                         float         base_noise_amp = 0.2f,
                         float         k_smoothing = 0.1f,
                         glm::vec2     center = {0.5f, 0.5f},
                         const Array  *p_noise_x = nullptr,
                         const Array  *p_noise_y = nullptr,
                         glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a heightmap representing a radial mountain range.
 *
 * This function creates a heightmap that simulates a mountain range emanating
 * radially from a specified center. The mountain range is influenced by various
 * noise parameters and control attributes.
 *
 * @param  shape              The dimensions of the output heightmap as a 2D
 *                            vector.
 * @param  kw                 The wave numbers (frequency components) as a 2D
 *                            vector.
 * @param  seed               The seed for random noise generation.
 * @param  half_width         The half-width of the radial mountain range,
 *                            controlling its spread. Default is 0.2f.
 * @param  angle_spread_ratio The ratio controlling the angular spread of the
 *                            mountain range. Default is 0.5f.
 * @param  center             The center point of the radial mountain range as
 *                            normalized coordinates within [0, 1]. Default is
 *                           {0.5f, 0.5f}.
 * @param  octaves            The number of octaves for fractal noise
 *                            generation. Default is 8.
 * @param  weight             The initial weight for noise contribution. Default
 *                            is 0.7f.
 * @param  persistence        The amplitude scaling factor for subsequent noise
 *                            octaves. Default is 0.5f.
 * @param  lacunarity         The frequency scaling factor for subsequent noise
 *                            octaves. Default is 2.0f.
 * @param  p_ctrl_param       Optional pointer to an array of control parameters
 *                            influencing the terrain generation.
 * @param  p_noise_x          Optional pointer to a precomputed noise array for
 *                            the X-axis.
 * @param  p_noise_y          Optional pointer to a precomputed noise array for
 *                            the Y-axis.
 * @param  p_angle            Optional pointer to an array to output the angle.
 * @param  bbox               The bounding box of the output heightmap in
 *                            normalized coordinates [xmin, xmax, ymin, ymax].
 *                            Default is {0.0f, 1.0f, 0.0f, 1.0f}.
 *
 * @return                    Array The generated heightmap representing the
 *                            radial mountain range.
 *
 * @note Only available if OpenCL is enabled.
 *
 * **Example**
 * @include ex_mountain_range_radial.cpp
 *
 * **Result**
 * @image html ex_mountain_range_radial.png
 */
Array mountain_range_radial(glm::ivec2    shape,
                            glm::vec2     kw,
                            std::uint32_t seed,
                            float         half_width = 0.2f,
                            float         angle_spread_ratio = 0.5f,
                            float         core_size_ratio = 1.f,
                            glm::vec2     center = {0.5f, 0.5f},
                            int           octaves = 8,
                            float         weight = 0.7f,
                            float         persistence = 0.5f,
                            float         lacunarity = 2.f,
                            const Array  *p_ctrl_param = nullptr,
                            const Array  *p_noise_x = nullptr,
                            const Array  *p_noise_y = nullptr,
                            const Array  *p_angle = nullptr,
                            glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a mountain-like heightmap with a flattened (stump-shaped)
 * peak.
 *
 * This function procedurally creates a terrain resembling a broad mountain with
 * a smooth, truncated summit. The result combines layered fractal noise (FBM),
 * Gaussian shaping, and Voronoi-based ridge formation to produce
 * natural-looking geomorphological features. The shape can be modulated by
 * directional displacement noise, smoothed using deposition simulation, and
 * refined through gamma and smoothing parameters.
 *
 * @param  shape          Dimensions of the output heightmap (width, height).
 * @param  seed           Random seed for noise generation.
 * @param  scale          Global scale factor controlling the size of features.
 * @param  octaves        Number of noise octaves used in FBM generation.
 * @param  peak_kw        Base spatial frequency of the main ridge/peak
 *                        structure.
 * @param  rugosity       Controls the roughness of base noise displacements.
 * @param  angle          Orientation angle (in degrees) of the main ridge
 *                        direction.
 * @param  k_smoothing    Smoothing coefficient applied to the Voronoi FBM
 *                        ridges.
 * @param  gamma          Gamma correction applied to the ridge intensity.
 * @param  add_deposition If true, applies an additional smoothing pass to
 *                        simulate sediment deposition.
 * @param  ridge_amp      Amplitude scaling factor for ridge prominence.
 * @param  base_noise_amp Amplitude scaling factor for base displacement noise.
 * @param  center         Normalized coordinates of the mountain’s center
 *                        (default domain: [0, 1]²).
 * @param  p_noise_x      Optional pointer to horizontal displacement noise
 *                        (nullptr for none).
 * @param  p_noise_y      Optional pointer to vertical displacement noise
 *                        (nullptr for none).
 * @param  bbox           Bounding box of the generated region (xmin, xmax,
 *                        ymin, ymax).
 *
 * @return                A 2D Array containing the generated mountain stump
 *                        heightmap.
 *
 * **Example**
 * @include ex_mountain_stump.cpp
 *
 * **Result**
 * @image html ex_mountain_stump.png
 * */
Array mountain_stump(glm::ivec2    shape,
                     std::uint32_t seed,
                     float         scale = 1.f,
                     int           octaves = 8,
                     float         peak_kw = 6.f,
                     float         rugosity = 0.f,
                     float         angle = 45.f,
                     float         k_smoothing = 0.f,
                     float         gamma = 0.25f,
                     bool          add_deposition = true,
                     float         ridge_amp = 0.75f,
                     float         base_noise_amp = 0.1f,
                     glm::vec2     center = {0.5f, 0.5f},
                     const Array  *p_noise_x = nullptr,
                     const Array  *p_noise_y = nullptr,
                     glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a synthetic "Tibesti" mountain heightmap.
 *
 * This function procedurally creates a 2D elevation field inspired by the
 * Tibesti mountains, a desert volcanic range. The terrain is built by combining
 * Gabor wavelets, fractal simplex noise, and a Gaussian envelope, producing
 * sharp volcanic peaks modulated by erosion-like patterns.
 *
 * @param  shape              Output array shape (resolution in x and y).
 * @param  seed               Random seed used for noise and wavelet generation.
 * @param  scale              Global scaling factor controlling the overall
 *                            mountain size.
 * @param  octaves            Number of octaves used in the fractal noise and
 *                            Gabor FBM.
 * @param  peak_kw            Base frequency of the peak features (scaled
 *                            internally by @p scale).
 * @param  rugosity           Controls roughness of the fractal noise (higher =
 *                            more irregular).
 * @param  angle              Orientation angle (degrees) for Gabor waves.
 * @param  angle_spread_ratio Spread ratio of the Gabor wave orientation
 *                            (controls anisotropy).
 * @param  gamma              Gamma correction factor applied to auxiliary noise
 *                            fields.
 * @param  add_deposition     If true, applies a smoothing fill step simulating
 *                            sediment deposition.
 * @param  bulk_amp           Amplitude of bulk uplift applied to the peaks
 *                            (affects normalization with noise weighting).
 * @param  base_noise_amp     Amplitude of the base displacement noise.
 * @param  center             Center of the mountain in normalized coordinates.
 * @param  p_noise_x          Optional pointer to external displacement noise
 *                            field (X-axis).
 * @param  p_noise_y          Optional pointer to external displacement noise
 *                            field (Y-axis).
 * @param  bbox               Bounding box of the generation domain in
 *                            normalized coordinates.
 *
 * @return                    Array containing the generated Tibesti mountain
 *                            heightmap.
 *
 * **Example**
 * @include ex_mountain_tibesti.cpp
 *
 * **Result**
 * @image html ex_mountain_tibesti.png
 */
Array mountain_tibesti(glm::ivec2    shape,
                       std::uint32_t seed,
                       float         scale = 1.f,
                       int           octaves = 8,
                       float         peak_kw = 20.f,
                       float         rugosity = 0.f,
                       float         angle = 30.f,
                       float         angle_spread_ratio = 0.25f,
                       float         gamma = 1.f,
                       bool          add_deposition = true,
                       float         bulk_amp = 1.f,
                       float         base_noise_amp = 0.1f,
                       glm::vec2     center = {0.5f, 0.5f},
                       const Array  *p_noise_x = nullptr,
                       const Array  *p_noise_y = nullptr,
                       glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generate a tectonic plate–like heightfield using Voronoi FBM and
 * directional talus projection.
 *
 * Combines Voronoi edge-distance FBM with a base noise displacement, then
 * applies a talus projection along a given direction and blends the result with
 * the original field.
 *
 * @param  shape          Grid resolution.
 * @param  kw             Base wavelength/frequency.
 * @param  seed           Random seed.
 * @param  talus          Talus strength for directional projection.
 * @param  direction      Projection direction.
 * @param  mix_ratio      Blend factor between raw Voronoi and projected plates.
 * @param  base_noise_amp Amplitude of the displacement noise.
 * @param  kw_multiplier  Frequency multiplier for base noise.
 * @param  rugosity       Noise roughness.
 * @param  bbox           Bounding box for noise evaluation.
 * @return                Generated plate heightfield.
 *
 * **Example**
 * @include ex_plates.cpp
 *
 * **Result**
 * @image html ex_plates.png
 */
Array plates(glm::ivec2    shape,
             glm::vec2     kw,
             std::uint32_t seed,
             float         talus,
             int           direction = 0,
             float         mix_ratio = 0.9f,
             float         base_noise_amp = 0.05f,
             float         kw_multiplier = 2.f,
             int           octaves = 8,
             float         rugosity = 0.f,
             glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Generates a synthetic "shattered peak" terrain heightmap.
 *
 * This function procedurally creates a 2D elevation field resembling a sharp,
 * fractured mountain peak. It combines a Gaussian pulse envelope, fractal noise
 * displacements, and Voronoi-based primitives with edge-distance shaping. The
 * result is a rugged peak structure with broken ridges and steep slopes.
 *
 * @param  shape          Output array shape (resolution in x and y).
 * @param  seed           Random seed used for noise and Voronoi generation.
 * @param  scale          Global scaling factor controlling the overall peak
 *                        size.
 * @param  peak_kw        Base frequency of the peak features (scaled internally
 *                        by @p scale).
 * @param  rugosity       Controls roughness of the fractal noise (higher = more
 *                        irregular).
 * @param  angle          Orientation angle (degrees) for directional
 *                        displacements.
 * @param  gamma          Gamma correction factor applied to the final
 *                        heightmap.
 * @param  add_deposition If true, applies a smoothing fill step simulating
 *                        sediment deposition.
 * @param  bulk_amp       Amplitude of bulk uplift applied to the peak
 *                        (internally overridden to 0.5f for normalization).
 * @param  base_noise_amp Amplitude of the base displacement noise.
 * @param  k_smoothing    Voronoi smoothing parameter (controls ridge
 *                        sharpness).
 * @param  center         Center of the peak in normalized coordinates.
 * @param  p_noise_x      Optional pointer to external displacement noise field
 *                        (X-axis).
 * @param  p_noise_y      Optional pointer to external displacement noise field
 *                        (Y-axis).
 * @param  bbox           Bounding box of the generation domain in normalized
 *                        coordinates.
 *
 * @return                Array containing the generated shattered peak
 *                        heightmap.
 *
 * **Example**
 * @include ex_shattered_peak.cpp
 *
 * **Result**
 * @image html ex_shattered_peak.png
 */
Array shattered_peak(glm::ivec2    shape,
                     std::uint32_t seed,
                     float         scale = 1.f,
                     int           octaves = 8,
                     float         peak_kw = 4.f,
                     float         rugosity = 0.f,
                     float         angle = 30.f,
                     float         gamma = 1.f,
                     bool          add_deposition = true,
                     float         bulk_amp = 0.3f,
                     float         base_noise_amp = 0.1f,
                     float         k_smoothing = 0.f,
                     glm::vec2     center = {0.5f, 0.5f},
                     const Array  *p_noise_x = nullptr,
                     const Array  *p_noise_y = nullptr,
                     glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

} // namespace hmap::gpu
