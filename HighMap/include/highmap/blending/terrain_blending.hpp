/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

#pragma once
#include "highmap/array.hpp"
#include "highmap/coord_frame.hpp"
#include "highmap/virtual_array.hpp"

namespace hmap
{

/**
 * @brief Blending modes for combining terrains with CoordFrame transforms.
 */
enum class TerrainBlendMode : int
{
  TBM_ADD,
  TBM_AVERAGE,
  TBM_DISTANCE,
  TBM_HEIGHT,
  TBM_MAX,
  TBM_MAX_SMOOTH,
  TBM_MIN,
  TBM_MIN_SMOOTH,
  TBM_MULTISCALE,
  TBM_POISSON,
  TBM_SEAMLESS,
  TBM_SLOPE,
  TBM_SUBTRACT,
};

/**
 * @brief Blend an upper terrain onto a lower terrain in-place, taking their
 * CoordFrame transforms into account.
 *
 * @param lower       Lower (base) terrain array, modified in-place.
 * @param frame_lower Coordinate frame of the lower terrain.
 * @param upper       Upper (overlay) terrain array.
 * @param frame_upper Coordinate frame of the upper terrain.
 * @param mode        Blending mode.
 * @param param       Mode-specific parameter (smoothness k for MIN_SMOOTH /
 *                    MAX_SMOOTH, normalized feathering width in (0, 0.5] for
 *                    DISTANCE).
 */
void blend_terrains(Array            &lower,
                    const CoordFrame &frame_lower,
                    const Array      &upper,
                    const CoordFrame &frame_upper,
                    TerrainBlendMode  mode,
                    float             param = 0.2f);

/**
 * @brief Blend an upper VirtualArray terrain onto a lower VirtualArray terrain
 * in-place, taking their CoordFrame transforms into account.
 *
 * @param lower       Lower (base) VirtualArray, modified in-place.
 * @param frame_lower Coordinate frame of the lower terrain.
 * @param upper       Upper (overlay) VirtualArray.
 * @param frame_upper Coordinate frame of the upper terrain.
 * @param mode        Blending mode.
 * @param param       Mode-specific parameter (smoothness k for MIN_SMOOTH /
 *                    MAX_SMOOTH, normalized feathering width in (0, 0.5] for
 *                    DISTANCE).
 * @param cm          Compute mode for tile execution.
 */
void blend_terrains(VirtualArray       &lower,
                    const CoordFrame   &frame_lower,
                    const VirtualArray &upper,
                    const CoordFrame   &frame_upper,
                    TerrainBlendMode    mode,
                    float               param = 0.2f,
                    const ComputeMode  &cm = ComputeMode{});

/**
 * @brief Height-aware terrain blending for Array inputs.
 */
void blend_terrains_height(Array            &lower,
                           const CoordFrame &frame_lower,
                           const Array      &upper,
                           const CoordFrame &frame_upper,
                           float             feather_width = 0.2f,
                           float             height_offset = 0.f,
                           float             height_range = 0.2f);

/**
 * @brief Height-aware terrain blending for VirtualArray inputs.
 */
void blend_terrains_height(VirtualArray       &lower,
                           const CoordFrame   &frame_lower,
                           const VirtualArray &upper,
                           const CoordFrame   &frame_upper,
                           float               feather_width = 0.2f,
                           float               height_offset = 0.f,
                           float               height_range = 0.2f,
                           const ComputeMode  &cm = ComputeMode{});

/**
 * @brief Multiscale (frequency-band) terrain blending for Array inputs.
 */
void blend_terrains_multiscale(Array            &lower,
                               const CoordFrame &frame_lower,
                               const Array      &upper,
                               const CoordFrame &frame_upper,
                               float             feather_width = 0.25f,
                               int               n_levels = 3,
                               int               base_radius = 16,
                               float             detail_gain = 1.f);

/**
 * @brief Multiscale (frequency-band) terrain blending for VirtualArray inputs
 * using a coarse global proxy for low-frequency bands and tile-local detail
 * blending with halo synchronization.
 */
void blend_terrains_multiscale(VirtualArray       &lower,
                               const CoordFrame   &frame_lower,
                               const VirtualArray &upper,
                               const CoordFrame   &frame_upper,
                               float               feather_width = 0.25f,
                               int                 n_levels = 3,
                               int                 base_radius = 16,
                               float               detail_gain = 1.f,
                               int                 coarse_size = 128,
                               const ComputeMode  &cm = ComputeMode{});

/**
 * @brief Poisson / gradient-domain terrain blending for Array inputs.
 */
void blend_terrains_poisson(Array            &lower,
                            const CoordFrame &frame_lower,
                            const Array      &upper,
                            const CoordFrame &frame_upper,
                            int               iterations = 400,
                            float             feather_width = 0.15f,
                            bool              mixed_gradients = true);

/**
 * @brief Poisson / gradient-domain terrain blending for VirtualArray inputs
 * using a two-level coarse-proxy + tile-local Red-Black SOR solver with halo
 * synchronization.
 */
void blend_terrains_poisson(VirtualArray       &lower,
                            const CoordFrame   &frame_lower,
                            const VirtualArray &upper,
                            const CoordFrame   &frame_upper,
                            int                 iterations = 400,
                            float               feather_width = 0.15f,
                            bool                mixed_gradients = true,
                            int                 coarse_size = 128,
                            int                 sync_interval = 25,
                            const ComputeMode  &cm = ComputeMode{});

/**
 * @brief Seamless boundary blending (harmonic MVC membrane) for Array inputs.
 */
void blend_terrains_seamless(Array            &lower,
                             const CoordFrame &frame_lower,
                             const Array      &upper,
                             const CoordFrame &frame_upper,
                             float             feather_width = 0.2f,
                             int               boundary_samples = 64);

/**
 * @brief Seamless boundary blending (harmonic MVC membrane) for VirtualArray
 * inputs.
 */
void blend_terrains_seamless(VirtualArray       &lower,
                             const CoordFrame   &frame_lower,
                             const VirtualArray &upper,
                             const CoordFrame   &frame_upper,
                             float               feather_width = 0.2f,
                             int                 boundary_samples = 64,
                             const ComputeMode  &cm = ComputeMode{});

/**
 * @brief Slope-aware terrain blending for Array inputs.
 */
void blend_terrains_slope(Array            &lower,
                          const CoordFrame &frame_lower,
                          const Array      &upper,
                          const CoordFrame &frame_upper,
                          float             feather_width = 0.2f,
                          float             slope_scale = 0.1f,
                          int               ir_slope = 2);

/**
 * @brief Slope-aware terrain blending for VirtualArray inputs.
 */
void blend_terrains_slope(VirtualArray       &lower,
                          const CoordFrame   &frame_lower,
                          const VirtualArray &upper,
                          const CoordFrame   &frame_upper,
                          float               feather_width = 0.2f,
                          float               slope_scale = 0.1f,
                          int                 ir_slope = 2,
                          const ComputeMode  &cm = ComputeMode{});

/**
 * @brief Transfer fine-scale terrain features from multiple source exemplar
 * layers onto a target terrain in-place, accounting for each layer's spatial
 * scale via its CoordFrame.
 *
 * Source layers are treated purely as feature exemplars (without spatial
 * alignment). Each source is scale-matched so its pixel resolution corresponds
 * to the target's physical cell size (`CoordFrame::get_size() / shape`), and
 * contributes atoms to a memory-bounded sparse-coding patch dictionary
 * (`max_atoms`). High-resolution features are synthesized across the entire
 * target terrain and integrated via high-pass detail transfer (`transfer`),
 * preserving the target's large-scale shape and elevation structure.
 *
 * @param target              Target terrain array, modified in-place.
 * @param frame_target        Coordinate frame of the target terrain.
 * @param sources             Source terrain exemplar arrays.
 * @param frames_source       Coordinate frames of the source terrains.
 * @param ir                  Smoothing radius (in target pixels) separating
 *                            large-scale structure from transferred details.
 * @param amplitude           Global amplitude scaling for transferred details.
 * @param target_prefiltering If true, smooth the target with radius @p ir
 *                            before adding transferred features.
 * @param source_weights      Optional per-source feature amplitude weights
 *                            (defaults to 1.0 for each source if empty).
 * @param factor              Super-resolution factor (>= 1).
 * @param patch_size          Low-resolution square patch size (>= 2).
 * @param analysis_stride     Patch extraction stride on downsampled sources.
 * @param synthesis_stride    Patch synthesis stride on downsampled target.
 * @param sparsity            Maximum atoms per patch in OMP encoding (>= 1).
 * @param max_atoms           Maximum total dictionary atoms across all sources.
 * @param seed                Random seed for reservoir sampling when candidate
 *                            patches exceed @p max_atoms.
 */
void transfer_terrain_features(
    Array                                 &target,
    const CoordFrame                      &frame_target,
    const std::vector<const Array *>      &sources,
    const std::vector<const CoordFrame *> &frames_source,
    int                                    ir = 16,
    float                                  amplitude = 1.f,
    bool                                   target_prefiltering = false,
    const std::vector<float>              &source_weights = {},
    int                                    factor = 4,
    int                                    patch_size = 8,
    int                                    analysis_stride = 4,
    int                                    synthesis_stride = 4,
    int                                    sparsity = 3,
    int                                    max_atoms = 1024,
    std::uint32_t                          seed = 0);

/**
 * @brief Single-source convenience overload of transfer_terrain_features for
 * Array inputs.
 */
void transfer_terrain_features(Array            &target,
                               const CoordFrame &frame_target,
                               const Array      &source,
                               const CoordFrame &frame_source,
                               int               ir = 16,
                               float             amplitude = 1.f,
                               bool              target_prefiltering = false,
                               int               factor = 4,
                               int               patch_size = 8,
                               int               analysis_stride = 4,
                               int               synthesis_stride = 4,
                               int               sparsity = 3,
                               int               max_atoms = 1024,
                               std::uint32_t     seed = 0);

/**
 * @brief Transfer fine-scale terrain features from multiple VirtualArray
 * source exemplar layers onto a target VirtualArray terrain in-place.
 *
 * Sources are streamed tile-by-tile to build a memory-bounded patch dictionary
 * (`max_atoms`), and target tiles are synthesized and detail-transferred
 * cooperatively with overlap synchronization via `ICommunicator`.
 */
void transfer_terrain_features(
    VirtualArray                            &target,
    const CoordFrame                        &frame_target,
    const std::vector<const VirtualArray *> &sources,
    const std::vector<const CoordFrame *>   &frames_source,
    int                                      ir = 16,
    float                                    amplitude = 1.f,
    bool                                     target_prefiltering = false,
    const std::vector<float>                &source_weights = {},
    int                                      factor = 4,
    int                                      patch_size = 8,
    int                                      analysis_stride = 4,
    int                                      synthesis_stride = 4,
    int                                      sparsity = 3,
    int                                      max_atoms = 1024,
    std::uint32_t                            seed = 0,
    const ComputeMode                       &cm = ComputeMode{});

/**
 * @brief Single-source convenience overload of transfer_terrain_features for
 * VirtualArray inputs.
 */
void transfer_terrain_features(VirtualArray       &target,
                               const CoordFrame   &frame_target,
                               const VirtualArray &source,
                               const CoordFrame   &frame_source,
                               int                 ir = 16,
                               float               amplitude = 1.f,
                               bool                target_prefiltering = false,
                               int                 factor = 4,
                               int                 patch_size = 8,
                               int                 analysis_stride = 4,
                               int                 synthesis_stride = 4,
                               int                 sparsity = 3,
                               int                 max_atoms = 1024,
                               std::uint32_t       seed = 0,
                               const ComputeMode  &cm = ComputeMode{});

} // namespace hmap
