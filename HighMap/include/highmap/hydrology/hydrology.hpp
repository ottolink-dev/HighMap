/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file hydrology.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for hydrological modeling functions and utilities.
 *
 * This header file declares functions and utilities for hydrological modeling,
 * including tools for computing flow directions, flow accumulation, and
 * identifying flow sinks within heightmaps. It supports multiple flow direction
 * models and the D8 model for flow direction and accumulation calculations.
 *
 * Key functionalities include:
 * - Computation of flow directions and flow accumulations using various models.
 * - Identification of flow sinks in heightmaps.
 * - Support for multiple flow direction models with customizable flow-partition
 * exponents.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once
#include <cstdint>
#include <functional>
#include <limits>

#include "highmap/array.hpp"
#include "highmap/geometry/path.hpp"
#include "highmap/hydrology/drainage_basin_cell_based.hpp"
#include "highmap/math/profiles.hpp"

namespace hmap
{

/**
 * @brief Label drainage basins using a priority-flood algorithm.
 *
 * Starting from boundary cells, this function floods a 2D scalar field inward
 * using a height-ordered priority queue and assigns each cell to a basin. Each
 * boundary cell initializes a unique basin.
 *
 * Connectivity is 8-neighbor. Spill elevations between basins are detected but
 * not yet recorded.
 *
 * @param  z Input 2D elevation array.
 * @return   Array of basin identifiers with the same shape as @p z.
 *
 * **Example**
 * @include ex_basin_id.cpp
 *
 * **Result**
 * @image html ex_basin_id.png
 */
Array basin_id(const Array        &z,
               FlowDirectionMethod fd_method = FlowDirectionMethod::FDM_D8,
               bool                remove_lakes = true);

/**
 * @brief Blends and carves smooth riverbanks around modified riverbed profiles.
 *
 * Compares an initial elevation terrain @p z with an updated/carved riverbed @p
 * z_river, identifies modified cells, expands riverbanks using talus slopes,
 * smooths the channel bottom, and blends the carved river with the original
 * terrain using distance-transform based exponential falloff and optional
 * domain warping.
 *
 * @param  z                     Original terrain elevation array.
 * @param  z_river               Carved riverbed elevation array (or stream
 * bottom).
 * @param  talus_riverbank       Talus slope angle/ratio for riverbank
 * expansion.
 * @param  smooth_river_bottom   Whether to apply Laplace smoothing along the
 * river bottom.
 * @param  merging_distance      Transition distance (in pixels) for blending
 * riverbanks into terrain.
 * @param  seed                  Random seed for riverbank noise.
 * @param  riverbank_noise_ratio Noise ratio for talus expansion roughness.
 * @param  p_noise_x             Optional noise array for domain warping in X.
 * @param  p_noise_y             Optional noise array for domain warping in Y.
 * @return                       Merged terrain array with smoothly carved river
 * channels.
 */
Array carve_riverbed(const Array  &z,
                     const Array  &z_river,
                     float         talus_riverbank = 0.01f,
                     bool          smooth_river_bottom = true,
                     float         merging_distance = 8.f,
                     std::uint32_t seed = 0,
                     float         riverbank_noise_ratio = 0.f,
                     const Array  *p_noise_x = nullptr,
                     const Array  *p_noise_y = nullptr);

/**
 * @brief Computes the number of drainage paths for each cell based on the D8
 * flow direction model.
 *
 * This function calculates the number of neighboring cells that have flow
 * directions pointing to the current cell according to the D8 model. The result
 * is an array where each cell contains the count of its incoming flow
 * directions.
 *
 * @param  d8 Input array representing the flow directions according to the D8
 *            model. Each cell value indicates the direction of flow according
 *            to the D8 convention.
 * @return    Array An array where each cell contains the number of incoming
 *            flow directions pointing to it.
 */
Array d8_compute_ndip(const Array &d8);

/**
 * @brief Identifies flow apex (source) cells using D8 flow routing.
 *
 * A flow apex is defined as a cell with zero inflow neighbors (no upstream
 * contributors) in the D8 flow graph. These cells correspond to drainage
 * sources, ridge heads, or local flow origins.
 *
 * @param z  Input elevation field.
 * @param is Output x-coordinates of flow apex cells.
 * @param js Output y-coordinates of flow apex cells.
 */
void find_flow_apex(const Array &z, std::vector<int> &is, std::vector<int> &js);

/**
 * @brief Identifies the indices of flow sinks within the heightmap.
 *
 * This function locates the cells in the heightmap that are flow sinks (cells
 * with no outflow) and returns their indices. Flow sinks are cells that do not
 * direct flow to any other cell.
 *
 * @param z  Input array representing the heightmap values.
 * @param is Output vector containing the row indices `i` of the flow sinks.
 * @param js Output vector containing the column indices `j` of the flow sinks.
 *
 * **Example**
 * @include ex_find_flow_sinks.cpp
 */
void find_flow_sinks(const Array      &z,
                     std::vector<int> &is,
                     std::vector<int> &js);

/**
 * @brief Find interior flow sinks in a height field.
 *
 * A sink is a cell whose 8 neighbors all have strictly higher elevation. Border
 * cells are excluded.
 *
 * @param  z Input elevation grid.
 * @return   List of sink cell coordinates (i, j).
 */
std::vector<glm::ivec2> find_flow_sinks(const Array &z);

/**
 * @brief Find flow sinks located on the domain border.
 *
 * A border cell is considered a sink if all its valid neighbors (within bounds)
 * have strictly higher elevation.
 *
 * @param  z Input elevation grid.
 * @return   List of border sink cell coordinates (i, j).
 */
std::vector<glm::ivec2> find_flow_sinks_border(const Array &z);

/**
 * @brief Compute water depth for a uniform flooding level.
 *
 * Subtracts terrain elevation from a reference water level and clamps negative
 * values to zero.
 *
 * @param  z    Input elevation array.
 * @param  zref Reference water level.
 * @return      Water depth array (0 where terrain is above zref).
 *
 * **Example**
 * @include ex_flooding_from_point.cpp
 *
 * **Result**
 * @image html ex_flooding_from_point.png
 */
Array flooding_uniform_level(const Array &z, float zref);

/**
 * @brief Compute flooding starting from the lowest boundary points.
 *
 * Finds the lowest elevation on the selected boundaries and simulates flooding
 * from those points up to a reference water level.
 *
 * @param  z          Input elevation array.
 * @param  zref       Reference water level.
 * @param  from_east  Include east boundary.
 * @param  from_west  Include west boundary.
 * @param  from_north Include north boundary.
 * @param  from_south Include south boundary.
 * @return            Water depth array.
 *
 * **Example**
 * @include ex_flooding_from_point.cpp
 *
 * **Result**
 * @image html ex_flooding_from_point.png
 */
Array flooding_from_boundaries(const Array &z,
                               float        zref,
                               bool         from_east = true,
                               bool         from_west = true,
                               bool         from_north = true,
                               bool         from_south = true);

/**
 * @brief Flood terrain starting from a single seed point.
 *
 * Fills areas below a reference level by propagating from the given cell until
 * higher ground is reached.
 *
 * @param  z         Input elevation array.
 * @param  i         Seed point X index.
 * @param  j         Seed point Y index.
 * @param  depth_min Water depth at the source point (if max, uses 0).
 * @return           Water depth array.
 *
 * **Example**
 * @include ex_flooding_from_point.cpp
 *
 * **Result**
 * @image html ex_flooding_from_point.png
 */
Array flooding_from_point(const Array &z,
                          int          i,
                          int          j,
                          float depth_min = std::numeric_limits<float>::max());

/**
 * @brief Flood terrain starting from multiple seed points.
 *
 * Merges flooding results from each point below the reference level.
 *
 * @param  z         Input elevation array.
 * @param  i         List of X indices for seeds.
 * @param  j         List of Y indices for seeds.
 * @param  depth_min Water depth at the source point (if max, uses 0).
 * @return           Water depth array.
 *
 * **Example**
 * @include ex_flooding_from_point.cpp
 *
 * **Result**
 * @image html ex_flooding_from_point.png
 */
Array flooding_from_point(const Array            &z,
                          const std::vector<int> &i,
                          const std::vector<int> &j,
                          float depth_min = std::numeric_limits<float>::max());

/**
 * @brief Estimate lake water depths on a terrain by filling depressions.
 *
 * This function identifies depressions in a terrain elevation model and
 * simulates the flooding of these areas to produce a lake system. It uses a
 * rough depression-filling algorithm to compute the water surface, then
 * subtracts the original elevations to obtain the water depth at each cell.
 *
 * @param  z                 Input 2D array representing terrain elevations
 *                           (height field).
 * @param  surface_threshold The minimum number of pixels a component must have
 *                           to be retained. Components smaller than this
 *                           threshold will be removed. The default value is 0
 *                           (no filtering).
 *
 * @return                   A 2D array representing the water depth for each
 *                           cell. Values are zero where no lake is present and
 *                           positive where water accumulates in depressions.
 *
 * @see                      depression_filling
 *
 * **Example**
 * @include ex_flooding_lake_system.cpp
 *
 * **Result**
 * @image html ex_flooding_lake_system.png
 */
Array flooding_lake_system(const Array &z, float surface_threshold = 0);

/**
 * @brief Computes the flow accumulation for each cell using the D8 flow
 * direction model.
 *
 * This function calculates the flow accumulation for each cell in the heightmap
 * based on the D8 flow direction model \cite Zhou2019. Flow accumulation
 * represents the total amount of flow that accumulates at each cell from
 * upstream cells.
 *
 * @param  z Input array representing the heightmap values.
 * @return   Array An array where each cell contains the computed flow
 *           accumulation.
 *
 * **Example**
 * @include ex_flow_accumulation_d8.cpp
 *
 * **Result**
 * @image html ex_flow_accumulation_d80.png
 * @image html ex_flow_accumulation_d81.png
 *
 * @see      flow_direction_d8
 */
Array flow_accumulation_d8(const Array &z);

/**
 * @brief Computes flow accumulation using a stochastic (Monte-Carlo) transport
 * estimator.
 *
 * Particles are spawned uniformly, advected along the downhill gradient field
 * with exact voxel traversal, attenuated by an optional decay term, and deposit
 * flux into every cell they exit. A final analytic normalization yields the
 * steady-state flux, including for cells no particle visited. Compared to
 * flow_accumulation_d8, the result is smooth and free of axis-aligned
 * artifacts.
 *
 * Port of the reference implementation of "Stochastic Geomorphological
 * Transport for Terrain Erosion Simulation" (N. McDonald, G. Cordonnier),
 * https://github.com/erosiv/geotransport (MIT).
 *
 * Results are reproducible up to floating-point summation order (deposits
 * accumulate via unordered atomic adds).
 *
 * @param  z         Input array representing the heightmap values.
 * @param  n_samples Number of Monte-Carlo samples (particles).
 * @param  seed      Random seed number (deterministic for fixed inputs).
 * @param  p_source  Optional per-cell source term (default: uniform 1).
 * @param  p_decay   Optional per-cell decay rate (default: 0, pure
 *                   accumulation).
 * @return           Array Per-cell steady-state flux.
 *
 * **Example**
 * @include ex_flow_accumulation_stochastic.cpp
 *
 * **Result**
 * @image html ex_flow_accumulation_stochastic0.png
 * @image html ex_flow_accumulation_stochastic1.png
 * @image html ex_flow_accumulation_stochastic2.png
 *
 * @see              flow_accumulation_d8
 */
Array flow_accumulation_stochastic(const Array  &z,
                                   int           n_samples = 1 << 19,
                                   std::uint32_t seed = 0,
                                   const Array  *p_source = nullptr,
                                   const Array  *p_decay = nullptr);

/**
 * @brief Computes the flow accumulation for each cell using the Multiple Flow
 * Direction (MFD) model.
 *
 * This function calculates the flow accumulation for each cell based on the MFD
 * model \cite Qin2007. Flow accumulation represents the total amount of flow
 * that accumulates at each cell from upstream cells. The flow-partition
 * exponent is locally defined using a reference talus value, where smaller
 * values of `talus_ref` will lead to narrower flow streams. The maximum talus
 * value of the heightmap can be used as a reference.
 *
 * @param  z         Input array representing the heightmap values.
 * @param  talus_ref Reference talus used to locally define the flow-partition
 *                   exponent. Small values will result in thinner flow streams.
 * @return           Array An array where each cell contains the computed flow
 *                   accumulation.
 *
 * **Example**
 * @include ex_flow_accumulation_dinf.cpp
 *
 * **Result**
 * @image html ex_flow_accumulation_dinf0.png
 * @image html ex_flow_accumulation_dinf1.png
 * @image html ex_flow_accumulation_dinf2.png
 *
 * @see              flow_direction_dinf, flow_accumulation_d8
 */
Array flow_accumulation_dinf(const Array &z, float talus_ref);

/**
 * @brief Computes averaged D-infinity flow accumulation with terrain
 * perturbations.
 *
 * This function estimates flow accumulation using the D∞ (D-infinity) method by
 * perturbing the input elevation map multiple times with Perlin noise. For each
 * sample, a noise field is added to the terrain, the flow accumulation is
 * computed, and results are accumulated and averaged.
 *
 * The perturbation introduces stochastic variability, which can help reduce
 * numerical artifacts and improve robustness in flat or ambiguous flow regions.
 *
 * @param  z         Input elevation array.
 * @param  talus_ref Reference talus angle (or slope threshold) used by the D∞
 *                   algorithm.
 * @param  nsamples  Number of noise-perturbed samples to average.
 * @param  kw        Noise frequency scaling factor.
 * @param  seed      Base random seed (incremented for each sample).
 * @param  amp       Amplitude of the noise perturbation.
 * @param  bbox      Bounding box used for noise generation.
 *
 * @return           Averaged flow accumulation array.
 *
 * **Example**
 * @include ex_flow_accumulation_dinf.cpp
 *
 * **Result**
 * @image html ex_flow_accumulation_dinf0.png
 * @image html ex_flow_accumulation_dinf1.png
 * @image html ex_flow_accumulation_dinf2.png
 */

Array flow_accumulation_dinf_perturbed(const Array  &z,
                                       float         talus_ref,
                                       int           nsamples,
                                       glm::vec2     kw,
                                       std::uint32_t seed,
                                       float         amp,
                                       const Array *p_perturb_scaling = nullptr,
                                       glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Computes the flow direction from each cell to its downslope neighbor
 * using the D8 model.
 *
 * This function calculates the direction of flow for each cell in the heightmap
 * using the D8 flow direction model @cite Greenlee1987. The D8 model defines
 * eight possible flow directions as follows:
 * @verbatim 5 6 7 4 . 0 3 2 1
 * @endverbatim
 *
 * @param  z Input array representing the heightmap values.
 * @return   Array An array where each cell contains the flow direction
 *           according to the D8 nomenclature.
 *
 * @see      flow_accumulation_d8
 */
Array flow_direction_d8(const Array &z);

/**
 * @brief Computes the flow direction and weights for each direction using the
 * Multiple Flow Direction (MFD) model.
 *
 * This function calculates the flow direction for each cell and provides the
 * weight for each possible flow direction using the MFD model \cite Qin2007.
 * The flow-partition exponent is determined using a reference talus value.
 * Smaller values of `talus_ref` will lead to thinner flow streams. The maximum
 * talus value of the heightmap can be used as a reference.
 *
 * @param  z         Input array representing the heightmap values.
 * @param  talus_ref Reference talus used to locally define the flow-partition
 *                   exponent. Smaller values will result in thinner flow
 *                   streams.
 * @return           std::vector<Array> A vector of arrays, each containing the
 *                   weights for flow directions at every cell.
 */
std::vector<Array> flow_direction_dinf(const Array &z, float talus_ref);

std::vector<float> flow_direction_dinf_flat(const Array &z, float talus_ref);

/**
 * @brief Computes the flow direction using the Multiple Flow Direction (MFD)
 * model.
 */
Array flow_direction_dinf_angle(const Array &z, float talus_ref);

Array flow_fixing(const Array &z,
                  float        riverbed_talus = 0.f,
                  int          iterations = 5,
                  int          prefilter_ir = 8,
                  bool         carve_riverbed = true,
                  float        merging_distance = 8.f, // pixels
                  const Array *p_noise_r = nullptr);

/**
 * @brief Fixes flow paths and unwanted upslopes using a cell-based drainage
 * basin tree network and elevation updates.
 *
 * This function constructs a stream tree network on the heightmap, resolves
 * internal depressions/lakes towards boundary outlets, and iteratively adjusts
 * terrain elevations along upstream-downstream paths to guarantee unbroken flow
 * while preserving realistic terrain slope gradients.
 *
 * @param  z                Input elevation array.
 * @param  iterations       Number of iterative tree-building and elevation
 *                          update passes.
 * @param  min_slope        Minimum allowed downslope gradient along flow paths.
 * @param  max_slope        Maximum allowed slope gradient.
 * @param  uplift_rate      Rate of elevation adjustment / uplift per unit
 *                          response time.
 * @param  m_exp            Drainage area exponent for response time
 *                          calculations.
 * @param  seed             Random seed for stochastic receiver variations.
 * @param  noise_strength   Strength of random perturbation during receiver
 *                          computation.
 * @param  carve_riverbed   Whether to apply riverbank carving and smoothing
 *                          along altered paths.
 * @param  talus_riverbank  Talus value used when expanding and carving the
 *                          riverbed.
 * @param  merging_distance Distance (in pixels) for blending modified flow
 *                          paths.
 * @param  p_noise_x        Optional noise array for riverbank domain warping in
 *                          X.
 * @param  p_noise_y        Optional noise array for riverbank domain warping in
 *                          Y.
 * @return                  Array with unbroken flow paths.
 */
Array flow_fixing_drainage_basin(
    const Array        &z,
    FlowDirectionMethod fd_method = FlowDirectionMethod::FDM_D8,
    float               riverbed_talus = 1e-4f,
    int                 iterations = 3,
    bool                carve_riverbed = false,
    float               talus_riverbank = 0.01f,
    float               merging_distance = 8.f,
    std::uint32_t       seed = 0,
    float               noise_strength = 0.f,
    const Array        *p_noise_x = nullptr,
    const Array        *p_noise_y = nullptr);

/**
 * @brief Resolves flow sinks and unwanted upslopes using Dijkstra pathfinding
 * and a Minimum Spanning Tree (MST) over sinks and boundary outlets.
 *
 * This function identifies all flow sinks and boundary outlets, computes
 * candidate optimal path connections via multi-source Dijkstra search,
 * constructs a Minimum Spanning Forest (Kruskal's algorithm) guaranteeing that
 * each sink connects to an outlet through the lowest-effort saddle pass, and
 * carves monotonic riverbeds along the MST paths using continuous trench
 * profiles.
 *
 * @param  z                     Input elevation array.
 * @param  riverbed_talus        Minimum talus (slope) along carved riverbeds.
 * @param  elevation_ratio       Balance factor between elevation level and
 *                               slope in Dijkstra cost.
 * @param  distance_exponent     Exponent applied to elevation differences in
 *                               Dijkstra cost.
 * @param  upward_penalization   Penalty factor for uphill moves in Dijkstra
 *                               search.
 * @param  valley_affinity       Weight for valley/concavity affinity in cost.
 * @param  prefilter_ir          Radius of Gaussian/cpulse prefilter applied
 *                               before sink detection.
 * @param  minimum_depth         Minimum incision depth below initial terrain
 *                               elevation.
 * @param  carve_riverbed        Whether to apply riverbank carving and
 *                               smoothing along altered paths.
 * @param  merging_distance      Distance (in pixels) for blending modified flow
 *                               paths.
 * @param  radial_profile        Radial profile cross-section for trench
 *                               carving.
 * @param  radial_profile_parameter Shape parameter for the radial profile.
 * @param  p_noise_r             Optional radial noise array for trench width
 *                               perturbation.
 * @param  fractalize_iterations Number of recursive midpoint displacement
 *                               iterations (0 to disable, > 0 to enable).
 * @param  fractalize_sigma      Displacement amplitude standard deviation
 *                               relative to edge length.
 * @param  decimate_target_points Number of target points to decimate paths down
 *                               to before fractalization (0 or <= 0 to skip
 *                               decimation).
 * @param  fractalize_seed       Random seed for path fractalization.
 * @return                       Array with unbroken flow paths.
 */
Array flow_fixing_mst(
    const Array  &z,
    float         riverbed_talus = 0.f,
    float         elevation_ratio = 0.95f,
    float         distance_exponent = 2.f,
    float         upward_penalization = 50.f,
    float         valley_affinity = 0.5f,
    int           prefilter_ir = 8,
    float         minimum_depth = 1e-4f,
    bool          carve_riverbed = true,
    float         merging_distance = 8.f,
    RadialProfile radial_profile = RadialProfile::RP_SMOOTHSTEP_UPPER,
    float         radial_profile_parameter = 2.f,
    const Array  *p_noise_r = nullptr,
    int           fractalize_iterations = 2,
    float         fractalize_sigma = 0.2f,
    int           decimate_target_points = 0,
    std::uint32_t fractalize_seed = 0);

/**
 * @brief Computes the optimal flow path from a starting point to the boundary
 * of a given elevation array.
 *
 * This function finds the flow path on a grid represented by the input array
 * `z`, starting from the given point `ij_start`. It identifies the best path to
 * the boundary, minimizing upward elevation penalties while accounting for
 * distance and elevation factors.
 *
 * @param  z                   The input 2D array representing elevation values.
 * @param  ij_start            The starting point as a 2D vector of indices (i,
 *                             j) within the array.
 * @param  elevation_ratio     Weight for elevation difference in the cost
 *                             function (default: 0.5).
 * @param  distance_exponent   Exponent for the distance term in the cost
 *                             function (default: 2.0).
 * @param  upward_penalization Penalty factor for upward elevation changes
 *                             (default: 100.0).
 * @return                     A Path object representing the optimal flow path
 *                             with normalized x and y coordinates and
 *                             corresponding elevations.
 *
 * The output path consists of:
 * - Normalized x-coordinates along the path.
 * - Normalized y-coordinates along the path.
 * - Elevation values corresponding to each point on the path.
 *
 * **Example**
 * @include ex_flow_stream.cpp
 *
 * **Result**
 * @image html ex_flow_stream.png
 */
Path flow_stream(const Array     &z,
                 const glm::ivec2 ij_start,
                 const float      elevation_ratio = 0.5f,
                 const float      distance_exponent = 2.f,
                 const float      upward_penalization = 100.f);

/**
 * @brief Generates a 2D array representing a riverbed based on a specified
 * path.
 *
 * This function calculates a scalar depth field (`dz`) for a riverbed shape
 * using a path, which can optionally be smoothed with Bezier curves. It
 * supports noise perturbation and post-filtering to adjust the riverbed's
 * features.
 *
 * @param  path                 The input path defining the riverbed's
 *                              trajectory.
 * @param  shape                The dimensions of the output array (width,
 *                              height).
 * @param  bbox                 The bounding box for the output grid in world
 *                              coordinates.
 * @param  bezier_smoothing     Flag to enable or disable Bezier smoothing of
 *                              the path.
 * @param  depth_start          The depth at the start of the riverbed.
 * @param  depth_end            The depth at the end of the riverbed.
 * @param  slope_start          The slope multiplier at the start of the
 *                              riverbed.
 * @param  slope_end            The slope multiplier at the end of the riverbed.
 * @param  shape_exponent_start The shape exponent at the start of the riverbed.
 * @param  shape_exponent_end   The shape exponent at the end of the riverbed.
 * @param  k_smoothing          The smoothing factor for the riverbed shape
 *                              adjustments.
 * @param  post_filter_ir       The radius of the post-filtering operation for
 *                              smoothing the output.
 * @param  p_noise_x            Optional pointer to a noise array for perturbing
 *                              the x-coordinates.
 * @param  p_noise_y            Optional pointer to a noise array for perturbing
 *                              the y-coordinates.
 * @param  p_noise_r            Optional pointer to a noise array for perturbing
 *                              the radial function.
 * @return                      A 2D array representing the calculated riverbed
 *                              depth field.
 *
 * @note The function requires the path to have at least two points. If the path
 * has fewer points, an empty array is returned with the given shape.
 *
 * **Example**
 * @include ex_generate_riverbed.cpp
 *
 * **Result**
 * @image html ex_generate_riverbed.png
 */
Array generate_riverbed(const Path &path,
                        glm::ivec2  shape,
                        glm::vec4   bbox = {0.f, 1.f, 0.f, 1.f},
                        bool        bezier_smoothing = false,
                        float       depth_start = 0.01f,
                        float       depth_end = 1.f,
                        float       slope_start = 64.f,
                        float       slope_end = 32.f,
                        float       shape_exponent_start = 1.f,
                        float       shape_exponent_end = 10.f,
                        float       k_smoothing = 0.5f,
                        int         post_filter_ir = 0,
                        Array      *p_noise_x = nullptr,
                        Array      *p_noise_y = nullptr,
                        Array      *p_noise_r = nullptr);

/**
 * @brief Merge two water depth fields.
 *
 * Computes the maximum (or smoothed maximum) of two depth arrays.
 *
 * @param  depth1   First water depth array.
 * @param  depth2   Second water depth array.
 * @param  k_smooth Smoothing parameter (0 = sharp max).
 * @return          Combined water depth array.
 */
Array merge_water_depths(const Array &depth1,
                         const Array &depth2,
                         float        k_smooth = 0.f);

/**
 * @brief Generate a spatial snow melting map from terrain features.
 *
 * Computes a melting factor in [0,1] from elevation, slope, and sun exposure.
 * Low elevations, steep slopes, and sun-facing orientations increase melting.
 *
 * @param  z                    Input terrain elevation.
 * @param  melt_start_elevation Elevation where melting starts (no melt below).
 * @param  melt_end_elevation   Elevation where melting reaches full effect.
 * @param  elevation_exp        Exponent controlling elevation melt falloff.
 * @param  sun_azimuth          Sun azimuth angle in degrees.
 * @param  sun_zenith           Sun zenith angle in degrees.
 * @param  aspect_strength      Weight of sun exposure contribution.
 * @param  slope_exp            Exponent controlling slope sensitivity.
 * @param  slope_strength       Weight of slope contribution.
 *
 * @return                      Melting map with values in [0,1], where higher
 *                              values indicate stronger melting.
 *
 * **Example**
 * @include ex_snow_simulation.cpp
 *
 * **Result**
 * @image html ex_snow_simulation.png
 */
Array snow_melting_map(const Array &z,
                       float        melt_start_elevation = 0.f,
                       float        melt_end_elevation = 0.5f,
                       float        elevation_strength = 1.f,
                       float        elevation_exp = 1.f,
                       float        sun_azimuth = 0.f,
                       float        sun_zenith = 60.f,
                       float        aspect_strength = 0.f,
                       float        slope_exp = 1.f,
                       float        slope_strength = 0.f);

/**
 * @brief Compute water depth over a masked terrain using harmonic
 * interpolation.
 *
 * This function estimates the water depth above a terrain surface by solving a
 * Laplace equation on the domain defined by @p mask. The solution is obtained
 * using the harmonic interpolation method with Successive Over-Relaxation
 * (SOR). The resulting water depth is given by the difference between the
 * interpolated surface and the original terrain elevation.
 *
 * @param  z              Input 2D array representing the terrain elevations
 *                        (height field).
 * @param  mask           Mask array of the same shape as @p z. Values greater
 *                        than @p mask_threshold define regions where water can
 *                        accumulate, while lower values represent boundaries or
 *                        fixed terrain.
 * @param  mask_threshold Threshold used to convert @p mask into a binary field
 *                        (0 or 1) for identifying water/terrain boundaries.
 * @param  iterations_max Maximum number of SOR iterations used in the harmonic
 *                        interpolation.
 * @param  tolerance      Convergence criterion: the algorithm stops if the
 *                        maximum absolute update between iterations is less
 *                        than this value.
 *
 * @return                A 2D array containing the computed water depth at each
 *                        grid cell. Depth values are non-negative where water
 *                        is present and zero where the mask indicates no water.
 *
 * @see                   harmonic_interpolation
 *
 * **Example**
 * @include ex_water_depth_from_mask.cpp
 *
 * **Result**
 * @image html ex_water_depth_from_mask.png
 */
Array water_depth_from_mask(const Array &z,
                            const Array &mask,
                            float        mask_threshold = 0.f,
                            int          iterations_max = 10000,
                            float        tolerance = 1e-2f);

/**
 * @brief Apply a drying factor to a water depth field.
 *
 * This function reduces the water depth values in-place by multiplying them by
 * a given @p dry_out_ratio. If a mask is provided, the drying is applied only
 * where the mask is non-zero, allowing selective drying of specific regions.
 *
 * @param water_depth   Reference to the 2D array containing water depth values.
 *                      The array is modified in-place.
 * @param dry_out_ratio Multiplicative factor (typically between 0 and 1) used
 *                      to scale down the water depth. A value of 1 leaves the
 *                      depth unchanged, while 0 removes all water.
 * @param p_mask        Optional pointer to a mask array of the same shape as @p
 *                      water_depth. Drying is applied only where the mask has
 *                      non-zero values. If `nullptr`, the factor is applied
 *                      uniformly to all cells.
 * @param depth_max     Maximum water depth, computed automatically by default.
 *
 * **Example**
 * @include ex_water_depth_dry_out.cpp
 *
 * **Result**
 * @image html ex_water_depth_dry_out.png
 */
void water_depth_dry_out(Array       &water_depth,
                         float        dry_out_ratio = 0.5f,
                         const Array *p_mask = nullptr,
                         float depth_max = std::numeric_limits<float>::max());

/**
 * @brief Simulates a local increase in water depth without global reflooding.
 *
 * This function propagates an additional water depth over a terrain elevation
 * map (`z`), starting only from cells that already contain water in the input
 * `water_depth` array. The propagation is restricted to neighboring cells of
 * higher elevation (8-connectivity), preventing the formation of entirely new
 * flooded basins disconnected from the original water region.
 *
 * Unlike water_depth_increase_with_flooding(), this function does not perform a
 * full flooding simulation of the domain. It is primarily intended for
 * expanding or softening an existing water mask slightly beyond the current
 * shoreline.
 *
 * @param  water_depth      Input array representing the initial water depth.
 * @param  z                Elevation array corresponding to the same grid as
 *                         `water_depth`.
 * @param  additional_depth Additional water depth to propagate from the
 *                          existing water region.
 * @return                  An Array containing the updated water depth after
 *                          local propagation.
 *
 * @note Water spreads only upward toward neighboring cells with higher
 * elevation. Existing disconnected dry basins are not flooded.
 *
 * @note This function is useful for generating extended shoreline masks,
 * erosion influence regions, or visually expanding water boundaries without
 * performing a physically complete flooding simulation.
 *
 * **Example**
 * @include ex_water_depth_increase.cpp
 *
 * **Result**
 * @image html ex_water_depth_increase.png
 */
Array water_depth_increase(const Array &water_depth,
                           const Array &z,
                           float        additional_depth);

/**
 * @brief Simulates a physical rise in water level with full domain flooding.
 *
 * This function increases the water level by `additional_depth` and recomputes
 * flooding over the terrain elevation map (`z`). Unlike water_depth_increase(),
 * the propagation is not restricted to the vicinity of the original water
 * region: water may spread through connected terrain and flood previously dry
 * areas if they become reachable under the new water level.
 *
 * The simulation behaves more like a physical flooding process where rising
 * water can overflow barriers, connect basins, and eventually inundate large
 * portions of the domain depending on terrain topology.
 *
 * @param  water_depth      Input array representing the initial water depth.
 * @param  z                Elevation array corresponding to the same grid as
 *                         `water_depth`.
 * @param  additional_depth Additional water depth used to raise the global
 *                          water level.
 * @return                  An Array containing the updated flooded state after
 *                          recomputing water propagation.
 *
 * @note Flood propagation may extend far beyond the original water mask and can
 * potentially inundate the entire connected terrain domain.
 *
 * @note This function is more physically plausible than water_depth_increase(),
 * and is suitable for terrain flooding simulations, overflow analysis, or
 * hydrological studies.
 *
 * **Example**
 * @include ex_water_depth_increase.cpp
 *
 * **Result**
 * @image html ex_water_depth_increase.png
 */
Array water_depth_increase_with_flooding(const Array &water_depth,
                                         const Array &z,
                                         float        additional_depth);
/**
 * @brief Compute the curvature of the water interface from a signed distance
 * field.
 *
 * The interface curvature is evaluated from the level set built from the water
 * depth mask. Optionally, curvature values can be extended away from the
 * interface using closest-point propagation.
 *
 * @param  water_depth                  Water depth field.
 * @param  prefilter_ir                 Prefilter radius used for curvature
 *                                      evaluation.
 * @param  extend_values_from_interface Extend interface curvature values to the
 *                                      full domain using closest boundary
 *                                      points.
 * @return                              Curvature field.
 *
 * **Example**
 * @include ex_water_frontier_curvature.cpp
 *
 * **Result**
 * @image html ex_water_frontier_curvature.png
 */
Array water_frontier_curvature(const Array &water_depth,
                               int          prefilter_ir,
                               bool extend_values_from_interface = false);

/**
 * @brief Generates a binary mask representing water presence.
 *
 * This version of the function converts the given water depth array into a
 * binary mask where non-zero values indicate the presence of water.
 *
 * @param  water_depth Input array representing water depth values.
 * @return             A binary Array where each cell is 1 if water is present,
 *                     0 otherwise.
 *
 * **Example**
 * @include ex_water_mask.cpp
 *
 * **Result**
 * @image html ex_water_mask.png
 */
Array water_mask(const Array &water_depth);

/**
 * @brief Computes a gradient-based water mask using an extended water depth
 * model.
 *
 * This function computes a smooth mask indicating regions that would be newly
 * flooded when the water depth is artificially increased by a specified
 * additional amount. It uses `water_depth_increase()` to simulate the spread of
 * water over the terrain.
 *
 * @param  water_depth      Input array representing the current water depth.
 * @param  z                Elevation array corresponding to the same grid as
 *                          `water_depth`.
 * @param  additional_depth The amount of additional water depth to simulate.
 * @return                  An Array representing the normalized water extension
 *                          mask, where values range from 0 to 1.
 *
 * **Example**
 * @include ex_water_mask.cpp
 *
 * **Result**
 * @image html ex_water_mask.png
 */
Array water_mask(const Array &water_depth,
                 const Array &z,
                 float        additional_depth);

} // namespace hmap

namespace hmap::gpu
{

/**
 * @brief Computes the omnidirectional coastal fetch for each cell.
 *
 * Casts @p ndirections evenly-spaced rays from each cell and averages the
 * distance to the first land obstacle (or domain boundary) across all
 * directions. High values indicate open, wave-exposed coastline;
 * low values indicate sheltered or enclosed areas such as harbors.
 *
 * @param  z              Input elevation array.
 * @param  ndirections    Number of ray directions (uniformly distributed over
 *                        2π). Higher values improve accuracy at the cost of
 *                        performance.
 * @param  p_compute_mask Optional boolean mask; fetch is only evaluated where
 *                        the mask is true. Pass nullptr to process the entire
 *                        domain.
 * @return                Array of mean fetch values (in grid cells).
 *
 * **Example**
 * @include ex_coastal_fetch.cpp
 *
 * **Result**
 * @image html ex_coastal_fetch.png
 */
Array coastal_fetch(const Array &z,
                    int          ndirections,
                    const Array *p_compute_mask = nullptr);

/**
 * @brief Computes directional coastal fetch weighted by alignment with a
 * prevailing wind or wave direction.
 *
 * Each ray's fetch contribution is weighted by
 * @f$ \max(0,\, \cos(\theta - \alpha))^e @f$, where @f$\theta@f$ is the ray
 * azimuth, @f$\alpha@f$ is @p angle, and
 * @f$e@f$ is @p directional_exp. Rays in the back hemisphere are discarded.
 * Increasing @p directional_exp narrows the effective lobe around the
 * prevailing direction.
 *
 * @param  z               Input elevation array.
 * @param  angle           Prevailing direction (radians, trigonometric
 *                         convention: 0 = East, π/2 = North).
 * @param  directional_exp Exponent applied to the cosine weight. Use 1 for a
 *                         full cosine lobe, 2 to approximate wind-energy
 *                         weighting, or higher values for a narrower
 *                         directional response.
 * @param  ndirections     Number of ray directions sampled over 2π.
 * @param  p_compute_mask  Optional boolean mask. Pass nullptr to process the
 *                         entire domain.
 * @return                 Array of directionally-weighted fetch values (in grid
 *                         cells).
 *
 * **Example**
 * @include ex_coastal_fetch.cpp
 *
 * **Result**
 * @image html ex_coastal_fetch.png
 */
Array coastal_fetch_directional(const Array &z,
                                float        angle,
                                float        directional_exp,
                                int          ndirections,
                                const Array *p_compute_mask = nullptr);

/**
 * @brief Approximates flow accumulation from a 2D velocity field.
 *
 * Advects a scalar accumulation field using an Eulerian transport scheme driven
 * by the velocity components (u, v). Repeated iterations approximate upstream
 * contributing area and highlight flow convergence (thalwegs).
 *
 * @param  u          X-component of the flow velocity field.
 * @param  v          Y-component of the flow velocity field.
 * @param  iterations Number of transport iterations to perform.
 *
 * @return            Flow accumulation field.
 */
Array flow_accumulation_from_velocity_field(const Array &u,
                                            const Array &v,
                                            int          iterations);

/**
 * @brief GPU (OpenCL) variant of hmap::flow_accumulation_stochastic — see the
 * CPU declaration for the algorithm description and parameters.
 */
Array flow_accumulation_stochastic(const Array  &z,
                                   int           n_samples = 1 << 19,
                                   std::uint32_t seed = 0,
                                   const Array  *p_source = nullptr,
                                   const Array  *p_decay = nullptr);

/*! @brief See hmap::flow_direction_d8 */
Array flow_direction_d8(const Array &z);

/**
 * @brief GPU hydraulic flow simulation using a virtual-pipes model; simulates
 * shallow-water transport over a height field using iterative flux computation
 * and water transport passes, with optional flux diffusion and post-simulation
 * depth clamping to remove thin residual water layers.
 *
 * @param  z                       Terrain height field.
 * @param  water_height            Global water scaling factor.
 * @param  depth_map               Initial relative water distribution.
 * @param  iterations              Number of simulation steps.
 * @param  dt                      Time step size.
 * @param  flux_diffusion          Enables flux smoothing.
 * @param  flux_diffusion_strength Strength of flux diffusion.
 * @param  dry_out_ratio           Ratio for removing thin remaining water
 *                                 layers.
 *
 * @return                         Array Final simulated water depth map.
 */
Array flow_simulation(const Array &z,
                      float        water_depth,
                      const Array &depth_map,
                      int          iterations,
                      float        dt = 0.5f,
                      bool         flux_diffusion = true,
                      float        flux_diffusion_strength = 0.01f,
                      float        dry_out_ratio = 0.f);

Array flow_simulation_viscous(const Array &z,
                              float        water_depth,
                              const Array &depth_map,
                              int          iterations,
                              float        dt = 0.5f,
                              float        dry_out_ratio = 0.f,
                              float        viscosity = 1.f,
                              float        power = 2.5f);

/*! @brief See hmap::generate_riverbed */
Array generate_riverbed(const Path &path,
                        glm::ivec2  shape,
                        glm::vec4   bbox = {0.f, 1.f, 0.f, 1.f},
                        bool        bezier_smoothing = false,
                        float       depth_start = 0.01f,
                        float       depth_end = 1.f,
                        float       slope_start = 64.f,
                        float       slope_end = 32.f,
                        float       shape_exponent_start = 1.f,
                        float       shape_exponent_end = 10.f,
                        float       k_smoothing = 0.5f,
                        int         post_filter_ir = 0,
                        Array      *p_noise_x = nullptr,
                        Array      *p_noise_y = nullptr,
                        Array      *p_noise_r = nullptr);

/**
 * @brief Simulate snow accumulation and redistribution over a terrain.
 *
 * Models snow fall, depth-dependent avalanching (above repose angle) and
 * thermal creep smoothing (below repose angle), with optional melting and
 * post-filtering. Snow progressively stabilizes into smooth, rounded forms.
 *
 * @param  z                     Heightfield of the underlying terrain.
 * @param  snow_depth            Total amount of snow to deposit.
 * @param  fall_map              Spatial modulation of snow deposition.
 * @param  melting_map           Spatial melting factor in [0,1].
 * @param  talus                 Local repose angle (talus slope).
 * @param  iterations            Number of simulation steps.
 * @param  dt                    Time step.
 * @param  fall_iterations_ratio Fraction of iterations used for snowfall.
 * @param  k_snow                Avalanche transport coefficient.
 * @param  k_visc                Thermal creep (sub-repose smoothing)
 *                               coefficient.
 * @param  k_melt_factor         Global melting strength.
 * @param  k_depth_ratio         Controls snow stiffening with depth.
 * @param  k_depth_slope_ratio   Depth influence on repose angle.
 * @param  post_filter           Apply final smoothing pass if true.
 * @param  thermal_talus_ratio   Controls creep strength relative to talus.
 *
 * @return                       Final snow depth field.
 *
 * **Example**
 * @include ex_snow_simulation.cpp
 *
 * **Result**
 * @image html ex_snow_simulation.png
 */
Array snow_simulation(const Array &z,
                      float        snow_depth,
                      const Array &fall_map,
                      const Array &melting_map,
                      const Array &talus,
                      int          iterations,
                      float        dt = 0.5f,
                      float        fall_iterations_ratio = 1.f,
                      float        k_snow = 0.5f,
                      float        k_visc = 0.1f,
                      float        k_melt_factor = 0.8f,
                      float        k_depth_ratio = 1.f,
                      float        k_depth_slope_ratio = 1.f,
                      bool         post_filter = true,
                      float        thermal_talus_ratio = 0.2f);

/**
 * @brief Filters water depth values using elevation data and a given radius.
 *
 * @param depth Water depth array to filter (in/out).
 * @param z     Elevation array used for filtering.
 * @param ir    Filter radius.
 *
 * **Example**
 * @include ex_water_depth_filter.cpp
 *
 * **Result**
 * @image html ex_water_depth_filter.png
 */
void water_depth_filter(Array       &depth,
                        const Array &z,
                        int          ir,
                        const Array *p_water_mask = nullptr,
                        bool         smooth_contour = false,
                        float        transition_ratio = 0.1f);

/*! @brief See hmap::water_frontier_curvature */
Array water_frontier_curvature(const Array &water_depth,
                               int          prefilter_ir,
                               bool extend_values_from_interface = false);

/*! @brief See hmap::water_depth_from_mask */
Array water_depth_from_mask(const Array &z,
                            const Array &mask,
                            float        mask_threshold = 0.f,
                            int          iterations_max = 10000,
                            float        tolerance = 1e-2f);

} // namespace hmap::gpu
