/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file synthesis.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Provides functions for synthesizing new heightmaps from input arrays.
 *
 * This header file defines functions for synthesizing new heightmaps by using
 * various techniques such as non-parametric sampling, patch quilting, and
 * blending. These functions allow for the generation of larger or modified
 * heightmaps based on input data, using methods like patch stitching and random
 * sampling.
 *
 * Functions provided:
 * - `non_parametric_sampling`: Generates a new heightmap using non-parametric
 * sampling.
 * - `quilting`: Synthesizes a new heightmap by stitching patches from a set of
 * input heightmaps.
 * - `quilting_blend`: Similar to `quilting`, but includes blending of patches.
 * - `quilting_expand`: Synthesizes a larger heightmap by expanding and
 * stitching patches from an input heightmap.
 * - `quilting_shuffle`: Similar to `quilting_expand`, but reshuffles patches to
 * generate a new heightmap with the same shape.
 * - `terrain_super_resolution`: Amplifies a low-resolution heightmap using
 * sparse patch coding against a high-resolution exemplar.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */

#pragma once

#include "highmap/array.hpp"

namespace hmap
{

/**
 * @brief Synthesize a new heightmap based on an input array using a
 * non-parametric sampling method.
 *
 * This method generates a new heightmap by sampling patches from the input
 * array non-parametrically. It is a slow process and is based on the technique
 * described in @cite Efros1999.
 *
 * @param  array           Input array from which patches are sampled.
 * @param  patch_shape     Shape of the patches used for sampling.
 * @param  seed            Random seed number for patch selection.
 * @param  error_threshold Threshold for patch selection based on error.
 * @return                 Array Resulting synthesized heightmap.
 *
 * **Example**
 * @include ex_non_parametric_sampling.cpp
 *
 * **Result**
 * @image html ex_non_parametric_sampling.png
 */
Array non_parametric_sampling(const Array  &array,
                              glm::ivec2    patch_shape,
                              std::uint32_t seed,
                              float         error_threshold = 0.1f);

/**
 * @brief Synthesize a new heightmap by stitching together small patches from
 * input heightmaps.
 *
 * This function creates a new heightmap by stitching patches from a set of
 * input heightmaps. The stitching process allows for flipping, rotating, and
 * transposing patches, and includes options for smoothing based on the overlap
 * between patches. This technique is based on @cite Efros2001.
 *
 * @param  p_arrays           Vector of pointers to input heightmaps.
 * @param  patch_base_shape   Shape of the patches to be used.
 * @param  tiling             Number of patches in each direction (x, y).
 * @param  overlap            Overlap between patches, in the range ]0, 1[.
 * @param  seed               Random seed number for patch placement.
 * @param  patch_flip         Allow flipping of patches (up-down, left-right).
 * @param  patch_rotate       Allow 90-degree rotation of patches (for square
 *                            patches).
 * @param  patch_transpose    Allow transposition of patches (for square
 *                            patches).
 * @param  filter_width_ratio Ratio of filter width to overlap length for
 *                            smoothing.
 * @return                    Array Resulting synthesized heightmap.
 *
 * **Example**
 * @include ex_quilting.cpp
 *
 * **Result**
 * @image html ex_quilting0.png
 * @image html ex_quilting1.png
 * @image html ex_quilting2.png
 * @image html ex_quilting3.png
 * @image html ex_quilting4.png
 */
Array quilting(const std::vector<const Array *> &p_arrays,
               glm::ivec2                        patch_base_shape,
               glm::ivec2                        tiling,
               float                             overlap,
               std::uint32_t                     seed,
               std::vector<Array *>              secondary_arrays = {},
               bool                              patch_flip = true,
               bool                              patch_rotate = true,
               bool                              patch_transpose = true,
               float                             filter_width_ratio = 0.25f);

/**
 * @brief Synthesize a new heightmap by stitching together small patches from a
 * list of input heightmaps.
 *
 * This function creates a new heightmap by stitching patches from a set of
 * input heightmaps, similar to `quilting`, but the patches are blended
 * together. The blending options include flipping, rotating, and transposing
 * patches, with smoothing based on the overlap between patches.
 *
 * @param  p_arrays           Vector of pointers to input heightmaps.
 * @param  patch_base_shape   Shape of the patches to be used.
 * @param  overlap            Overlap between patches, in the range ]0, 1[.
 * @param  seed               Random seed number for patch placement.
 * @param  patch_flip         Allow flipping of patches (up-down, left-right).
 * @param  patch_rotate       Allow 90-degree rotation of patches (for square
 *                            patches).
 * @param  patch_transpose    Allow transposition of patches (for square
 *                            patches).
 * @param  filter_width_ratio Ratio of filter width to overlap length for
 *                            smoothing.
 * @return                    Array Resulting synthesized heightmap.
 *
 * **Example**
 * @include ex_quilting.cpp
 *
 * **Result**
 * @image html ex_quilting0.png
 * @image html ex_quilting1.png
 * @image html ex_quilting2.png
 * @image html ex_quilting3.png
 * @image html ex_quilting4.png
 */
Array quilting_blend(const std::vector<const Array *> &p_arrays,
                     glm::ivec2                        patch_base_shape,
                     float                             overlap,
                     std::uint32_t                     seed,
                     bool                              patch_flip = true,
                     bool                              patch_rotate = true,
                     bool                              patch_transpose = true,
                     float filter_width_ratio = 0.25f);

/**
 * @brief Synthesize a new heightmap by expanding the input heightmap and
 * stitching patches.
 *
 * This function generates a larger heightmap by expanding the input heightmap
 * and stitching patches based on the given parameters. The expansion ratio
 * determines the new size of the heightmap, and the function includes options
 * for flipping, rotating, and transposing patches. The output can be either the
 * same shape as the input or expanded based on the given ratio.
 *
 * @param  array              Input heightmap to be expanded.
 * @param  expansion_ratio    Ratio by which to expand the heightmap (e.g., 2
 *                            for doubling the size).
 * @param  patch_base_shape   Shape of the patches to be used.
 * @param  overlap            Overlap between patches, in the range ]0, 1[.
 * @param  seed               Random seed number for patch placement.
 * @param  keep_input_shape   If true, the output heightmap retains the input
 *                            shape.
 * @param  patch_flip         Allow flipping of patches (up-down, left-right).
 * @param  patch_rotate       Allow 90-degree rotation of patches (for square
 *                            patches).
 * @param  patch_transpose    Allow transposition of patches (for square
 *                            patches).
 * @param  filter_width_ratio Ratio of filter width to overlap length for
 *                            smoothing.
 * @return                    Array Resulting synthesized heightmap.
 *
 * **Example**
 * @include ex_quilting.cpp
 *
 * **Result**
 * @image html ex_quilting0.png
 * @image html ex_quilting1.png
 * @image html ex_quilting2.png
 * @image html ex_quilting3.png
 * @image html ex_quilting4.png
 */
Array quilting_expand(const Array         &array,
                      float                expansion_ratio,
                      glm::ivec2           patch_base_shape,
                      float                overlap,
                      std::uint32_t        seed,
                      std::vector<Array *> secondary_arrays = {},
                      bool                 keep_input_shape = false,
                      bool                 patch_flip = true,
                      bool                 patch_rotate = true,
                      bool                 patch_transpose = true,
                      float                filter_width_ratio = 0.25f);

/**
 * @brief Synthesize a new heightmap by reshuffling patches of the input
 * heightmap.
 *
 * This function generates a new heightmap by reshuffling patches of the input
 * heightmap, effectively creating a new pattern while keeping the same shape as
 * the input heightmap. The function includes options for patch flipping,
 * rotating, and transposing.
 *
 * @param  array              Input heightmap to be reshuffled.
 * @param  patch_base_shape   Shape of the patches to be used.
 * @param  overlap            Overlap between patches, in the range ]0, 1[.
 * @param  seed               Random seed number for patch placement.
 * @param  patch_flip         Allow flipping of patches (up-down, left-right).
 * @param  patch_rotate       Allow 90-degree rotation of patches (for square
 *                            patches).
 * @param  patch_transpose    Allow transposition of patches (for square
 *                            patches).
 * @param  filter_width_ratio Ratio of filter width to overlap length for
 *                            smoothing.
 * @return                    Array Resulting synthesized heightmap with
 *                            reshuffled patches.
 *
 * **Example**
 * @include ex_quilting.cpp
 *
 * **Result**
 * @image html ex_quilting0.png
 * @image html ex_quilting1.png
 * @image html ex_quilting2.png
 * @image html ex_quilting3.png
 * @image html ex_quilting4.png
 */
Array quilting_shuffle(const Array         &array,
                       glm::ivec2           patch_base_shape,
                       float                overlap,
                       std::uint32_t        seed,
                       std::vector<Array *> secondary_arrays = {},
                       bool                 patch_flip = true,
                       bool                 patch_rotate = true,
                       bool                 patch_transpose = true,
                       float                filter_width_ratio = 0.25f);

/**
 * @brief Synthesize a high-resolution heightmap from a low-resolution input
 * using an exemplar-based sparse patch representation (terrain
 * amplification).
 *
 * The exemplar is downsampled by `factor` and every `patch_size` x
 * `patch_size` patch of the downsampled exemplar (mean removed, radial mask
 * applied, unit norm) becomes an atom of a dictionary, paired with the
 * corresponding `patch_size * factor` high-resolution patch. Each patch of the
 * input is then sparse-coded against the low-resolution atoms with orthogonal
 * matching pursuit, and the same coefficients are applied to the paired
 * high-resolution atoms to synthesize the output, which is blended from the
 * overlapping patches using the mask as weight. This technique is based on
 * @cite Guerin2016.
 *
 * With `sparsity = 1` (the setting of the reference implementation), each
 * input patch is replaced by the detail of its best-correlated exemplar
 * patch. The strides control the density of the patch grids: smaller strides
 * give more atoms and more overlap at a higher cost. The patch mask vanishes
 * at the patch corners, so `synthesis_stride` should not exceed about
 * `0.7 * patch_size` to keep the output fully covered; uncovered pixels, if
 * any, are filled with a bicubic upsampling of the input.
 *
 * Memory usage scales with the number of atoms times `(patch_size *
 * factor)^2`: large exemplars with small `analysis_stride` can require a lot
 * of memory.
 *
 * @param  array            Low-resolution input heightmap.
 * @param  exemplar         High-resolution exemplar heightmap providing the
 *                          detail.
 * @param  factor           Amplification factor (>= 1); the output shape is
 *                          `array.shape * factor`.
 * @param  patch_size       Side of the square patches on the low-resolution
 *                          grid (>= 2).
 * @param  analysis_stride  Stride of the patch grid over the downsampled
 *                          exemplar, in `[1, patch_size]`.
 * @param  synthesis_stride Stride of the patch grid over the input, in `[1,
 *                          patch_size]`.
 * @param  sparsity         Maximum number of atoms used to code each patch
 *                          (>= 1).
 * @return                  Array Synthesized heightmap of shape `array.shape *
 *                          factor`, or an empty array if the inputs are
 *                          invalid (e.g. the downsampled exemplar is smaller
 *                          than a patch, or the exemplar is flat).
 *
 * **Example**
 * @include ex_terrain_super_resolution.cpp
 *
 * **Result**
 * @image html ex_terrain_super_resolution.png
 */
Array terrain_super_resolution(const Array &array,
                               const Array &exemplar,
                               int          factor,
                               int          patch_size = 16,
                               int          analysis_stride = 8,
                               int          synthesis_stride = 8,
                               int          sparsity = 1);

} // namespace hmap
