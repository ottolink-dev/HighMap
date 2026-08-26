/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file filters.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for filter-related functions and classes.
 *
 * This header file contains declarations for various filter functions and
 * classes used within the project. These filters are typically used for
 * processing data streams or collections according to specific criteria.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once

#include "highmap/array.hpp"
#include "highmap/primitives.hpp"

namespace hmap
{

/**
 * @enum neighborhood
 * @brief Enum representing different types of neighborhood lattices.
 *
 * This enumeration defines the different types of neighborhoods that can be
 * used in a lattice-based system. These neighborhoods determine how cells or
 * nodes are connected to their immediate surroundings.
 */
// clang-format off
enum neighborhood : int
{
	MOORE, ///< Moore neighborhood: includes all eight surrounding cells.
	VON_NEUMANN, ///< Von Neumann neighborhood: includes only the four
	// orthogonal neighbors (N, S, E, W).
	CROSS  ///< Cross-shaped neighborhood: includes only the diagonal
	// neighbors.
};
// clang-format on

/**
 * @brief Enhances local terrain relief adaptively while preventing overshoot.
 *
 * Applies an adaptive relief enhancement filter inspired by edge-adaptive
 * sharpening. Detail is enhanced proportionally to local activity, while steep
 * cliffs/high gradients are protected from over-sharpening artifacts and values
 * are bounded by local neighborhood min/max limits.
 *
 * @param array       Input/output heightmap array to enhance in-place.
 * @param strength    Sharpening/enhancement strength (default: 1.0f).
 * @param clamp_ratio Ratio controlling how strictly overshoot beyond local
 *                    extrema is clamped [0, 1] (1.0 = strict clamping to local
 *                    min/max).
 *
 * **Example**
 * @include ex_adaptive_relief.cpp
 *
 * **Result**
 * @image html ex_adaptive_relief.png
 */
void adaptive_relief(Array &array,
                     float  strength = 1.f,
                     float  clamp_ratio = 1.f);

/**
 * @brief Enhances local terrain relief adaptively with a mask.
 *
 * @param array       Input/output heightmap array to enhance in-place.
 * @param p_mask      Optional filter mask controlling local effect blending.
 * @param strength    Sharpening/enhancement strength (default: 1.0f).
 * @param clamp_ratio Ratio controlling how strictly overshoot beyond local
 *                    extrema is clamped [0, 1].
 *
 * @overload
 */
void adaptive_relief(Array       &array,
                     const Array *p_mask,
                     float        strength = 1.f,
                     float        clamp_ratio = 1.f);

/**
 * @brief Applies a primitive-based displacement to an array.
 *
 * Adds a scaled primitive pattern (optionally perturbed by noise) to the input
 * array.
 *
 * @param  z              Input array.
 * @param  primitive_type Type of primitive to generate.
 * @param  amp            Amplitude scaling factor.
 * @param  p_noise_x      Optional noise array for X perturbation (can be
 *                        nullptr).
 * @param  p_noise_y      Optional noise array for Y perturbation (can be
 *                        nullptr).
 * @param  center         Center of the primitive.
 * @param  bbox           Bounding box of the primitive.
 *
 * @return                Resulting array after applying the displacement.
 *
 * **Example**
 * @include ex_bulkify.cpp
 *
 * **Result**
 * @image html ex_bulkify.png
 */
Array bulkify(const Array         &z,
              const PrimitiveType &primitive_type,
              float                amp = 1.f,
              const Array         *p_noise_x = nullptr,
              const Array         *p_noise_y = nullptr,
              glm::vec2            center = {0.5f, 0.5f},
              glm::vec4            bbox = {0.f, 1.f, 0.f, 1.f});

/**
 * @brief Reduce quantization artifacts using dithering and smoothing.
 *
 * Adds white noise dithering to the input array and applies iterative Laplacian
 * filtering to smooth quantization steps.
 *
 * @param  array                Input quantized array.
 * @param  seed                 Random seed used for dithering noise generation.
 * @param  dither_amplitude     Amplitude of the added white noise.
 * @param  filtering_iterations Number of Laplacian smoothing iterations.
 * @return                      Dequantized and smoothed array.
 *
 * **Example**
 * @include ex_dequantize.cpp
 *
 * **Result**
 * @image html ex_dequantize.png
 */
Array dequantize(const Array  &array,
                 std::uint32_t seed,
                 float         dither_amplitude = 0.01f,
                 int           filtering_iterations = 1);

/**
 * @brief Apply histogram equalization to the array values.
 *
 * This function performs histogram equalization on the input array to enhance
 * the contrast of the image or data by adjusting the intensity distribution.
 *
 * @param array Input array to be equalized.
 */
void equalize(Array &array);

/**
 * @brief Apply histogram equalization to the array values with a mask.
 *
 * This function performs histogram equalization on the input array while
 * considering the specified mask. Only the elements of the array corresponding
 * to non-zero values in the mask are equalized.
 *
 * @param array  Input array to be equalized.
 * @param p_mask Optional mask array. Only the elements where the mask is
 *               non-zero are equalized.
 *
 * @overload
 */
void equalize(Array &array, const Array *p_mask);

/**
 * @brief Apply expansion, or "inflation", to emphasize the bulk of the terrain.
 *
 * This function expands the features of the terrain represented by the input
 * array. The expansion operation, often referred to as "inflation", increases
 * the prominence of elevated regions, effectively making them more pronounced.
 *
 * @param array  Input array representing the terrain or data to be expanded.
 * @param ir     Filter radius, which determines the extent of the expansion.
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               The mask controls which parts of the array are affected by the
 *               expansion.
 *
 * **Example**
 * @include ex_expand.cpp
 *
 * **Result**
 * @image html ex_expand.png
 */
void expand(Array &array, int ir, const Array *p_mask, int iterations = 1);

/**
 * @brief Apply expansion to emphasize the bulk of the terrain using a filter
 * radius.
 *
 * This overload applies expansion without the use of a mask, utilizing only the
 * filter radius.
 *
 * @param array Input array representing the terrain or data to be expanded.
 * @param ir    Filter radius, which determines the extent of the expansion.
 *
 * @overload
 */
void expand(Array &array, int ir, int iterations = 1);

/**
 * @brief Apply expansion using a custom kernel.
 *
 * This overload allows the use of a custom kernel for expansion, offering more
 * control over the operation compared to the filter radius method.
 *
 * @param array  Input array representing the terrain or data to be expanded.
 * @param kernel Custom kernel array used to perform the expansion.
 *
 * @overload
 */
void expand(Array &array, const Array &kernel, int iterations = 1);

/**
 * @brief Apply expansion using a custom kernel with an optional mask.
 *
 * This overload allows the use of a custom kernel for expansion, along with an
 * optional mask to control which parts of the array are affected.
 *
 * @param array  Input array representing the terrain or data to be expanded.
 * @param kernel Custom kernel array used to perform the expansion.
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               The mask controls which parts of the array are affected by the
 *               expansion.
 *
 * @overload
 */
void expand(Array       &array,
            const Array &kernel,
            const Array *p_mask,
            int          iterations = 1);

/**
 * @brief Apply expansion, or "inflation", to emphasize the bulk of the terrain,
 * using a directional kernel.
 *
 * This function expands features of the terrain represented by the input array,
 * with emphasis applied directionally based on the specified angle. The
 * directional kernel allows for anisotropic expansion, making the expansion
 * more pronounced in a specified direction.
 *
 * @param array        Input array representing the terrain or data to be
 *                     expanded.
 * @param ir           Filter radius, determining the extent of the expansion.
 * @param angle        Angle (in degrees) that sets the direction of the
 *                     expansion.
 * @param aspect_ratio Pulse aspect ratio, determining the shape of the
 *                     expansion.
 * @param anisotropy   Pulse width ratio between upstream and downstream sides,
 *                     default is 1.0 (isotropic expansion).
 * @param p_mask       Optional filter mask with values expected in the range
 *                     [0, 1]. The mask controls which parts of the array are
 *                     affected by the expansion.
 *
 * **Example**
 * @include ex_expand_directional.cpp
 *
 * **Result**
 * @image html ex_expand_directional.png
 */
void expand_directional(Array       &array,
                        int          ir,
                        float        angle,
                        float        aspect_ratio,
                        float        anisotropy = 1.f,
                        const Array *p_mask = nullptr);

/**
 * @brief Expand heights outward from a masked region with a talus slope limit.
 *
 * Starting from the non-zero cells of `mask`, heights are propagated outward
 * using a min-heap so lower cells are processed first. Neighbor heights are
 * computed from the source height plus `dist * talus`, optionally perturbed
 * with noise.
 *
 * Intuition: grows terrain outward from seed cells, filling empty space while
 * respecting a maximum slope.
 *
 * @param z           Height array modified in-place.
 * @param mask        Initial region to expand from.
 * @param talus       Maximum slope per unit distance.
 * @param seed        RNG seed for noise.
 * @param ir          Neighborhood radius.
 * @param noise_ratio Amplitude of slope noise.
 */
void expand_talus(Array        &z,
                  const Array  &mask,
                  float         talus,
                  std::uint32_t seed,
                  int           ir = 1,
                  float         noise_ratio = 0.2f);
/**
 * @brief Enforce a talus slope constraint on a height field.
 *
 * Cells are processed from highest to lowest using a max-heap. Each cell may
 * raise neighboring cells so that the slope relative to the current cell does
 * not exceed the local talus value.
 *
 * Intuition: stabilizes an existing terrain by redistributing height so slopes
 * do not exceed the talus limit (similar to talus erosion / relaxation but by
 * adding material).
 *
 * @param z           A reference to the 2D array representing the terrain
 *                    heights. The function modifies this array in place to
 *                    introduce talus slopes.
 * @param talus       The critical slope angle that determines where material
 *                    will move from higher elevations to lower ones. Slopes
 *                    steeper than this value will be flattened by material
 *                    transport.
 * @param seed        The seed for the random number generator, ensuring
 *                    reproducibility of the noise effects in the talus
 *                    formation process. The same seed will produce the same
 *                    terrain modifications.
 * @param noise_ratio A parameter that controls the amount of randomness or
 *                    noise introduced in the talus formation process. The
 *                    default value is 0.2.
 *
 * **Example**
 * @include ex_fill_talus.cpp
 *
 * **Result**
 * @image html ex_fill_talus.png
 *
 * See unit tests: @ref test_fill_talus.cpp
 */
void fill_talus(Array        &z,
                float         talus,
                std::uint32_t seed,
                int           ir = 1,
                float         noise_ratio = 0.2f,
                const Array  *p_seed_mask = nullptr);

void fill_talus(Array        &z,
                const Array  &talus,
                std::uint32_t seed,
                int           ir = 1,
                float         noise_ratio = 0.2f,
                const Array  *p_seed_mask = nullptr);

/**
 * @brief Apply a "folding" filter (successive absolute values) to the array
 * elements.
 *
 * This function applies a folding filter to the input array, where the values
 * are successively folded (using absolute values) within a specified range,
 * defined by `vmin` and `vmax`. The process can be repeated for a number of
 * iterations, with an optional smoothing parameter `k` to control the intensity
 * of the folding effect.
 *
 * @param array      Input array to which the folding filter is applied.
 * @param vmin       Minimum reference value used as a lower bound for the
 *                   folding operation.
 * @param vmax       Maximum reference value used as an upper bound for the
 *                   folding operation.
 * @param iterations Number of iterations for applying the folding filter. The
 *                   default is 3.
 * @param k          Smoothing parameter for the absolute value operation,
 *                   expected to be greater than 0. It controls the degree of
 *                   smoothing applied during folding. The default value is
 *                   0.05.
 *
 * **Example**
 * @include ex_fold.cpp
 *
 * **Result**
 * @image html ex_fold.png
 */
void fold(Array &array,
          float  vmin,
          float  vmax,
          int    iterations = 3,
          float  k = 0.05f);

/**
 * @brief Apply a "folding" filter with default reference values and parameters.
 *
 * This overload applies the folding filter to the input array with the default
 * number of iterations and smoothing parameter `k`. The folding is done without
 * explicit minimum and maximum reference values.
 *
 * @param array      Input array to which the folding filter is applied.
 * @param iterations Number of iterations for applying the folding filter. The
 *                   default is 3.
 * @param k          Smoothing parameter for the absolute value operation,
 *                   expected to be greater than 0. It controls the degree of
 *                   smoothing applied during folding. The default value is
 *                   0.05.
 *
 * @overload
 */
void fold(Array &array, int iterations = 3, float k = 0.05f);

/**
 * @brief Apply a gain correction to the array elements.
 *
 * This function applies a gain correction to the input array using a power law.
 * The correction enhances or reduces the intensity of the values based on the
 * specified gain factor.
 *
 * @param array  Input array to which the gain correction is applied.
 * @param factor Gain factor, expected to be greater than 0. It determines the
 *               strength of the gain correction.
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               The mask controls which parts of the array are affected by the
 *               gain correction.
 *
 * @warning Array values are expected to be in the range [0, 1]. Applying the
 * gain correction outside of this range may lead to unexpected results.
 *
 * **Example**
 * @include ex_gain.cpp
 *
 * **Result**
 * @image html ex_gain.png
 */
void gain(Array &array, float factor, const Array *p_mask);

/**
 * @brief Apply a gain correction to the array elements without a mask.
 *
 * This overload applies the gain correction to the entire array without using a
 * mask. The gain correction is performed using the specified gain factor.
 *
 * @param array  Input array to which the gain correction is applied.
 * @param factor Gain factor, expected to be greater than 0. It determines the
 *               strength of the gain correction.
 *
 * @overload
 */
void gain(Array &array, float factor);

/**
 * @brief Apply gamma correction to the input array.
 *
 * This function applies gamma correction to the input array, which adjusts the
 * brightness and contrast based on a specified gamma factor. Gamma correction
 * is used to correct or adjust the tonal range of an image or data to achieve a
 * desired visual effect.
 *
 * @param array  Input array to which the gamma correction is applied.
 * @param gamma  Gamma factor, expected to be greater than 0. It determines the
 *               level of correction applied. A gamma value less than 1 will
 *               lighten the image, while a value greater than 1 will darken it.
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               The mask controls which parts of the array are affected by the
 *               gamma correction.
 *
 * @warning Array values are expected to be in the range [0, 1]. Applying gamma
 * correction outside of this range may result in incorrect or unintended
 * modifications.
 *
 * **Example**
 * @include ex_gamma_correction.cpp
 *
 * **Result**
 * @image html ex_gamma_correction.png
 */
void gamma_correction(Array &array, float gamma, const Array *p_mask);

/**
 * @brief Apply gamma correction to the input array without a mask.
 *
 * This overload applies gamma correction to the entire array without using a
 * mask. The gamma correction is performed using the specified gamma factor.
 *
 * @param array Input array to which the gamma correction is applied.
 * @param gamma Gamma factor, expected to be greater than 0. It determines the
 *              level of correction applied. A gamma value less than 1 will
 *              lighten the image, while a value greater than 1 will darken it.
 *
 * @overload
 */
void gamma_correction(Array &array, float gamma);

/**
 * @brief Apply a "local" gamma correction to the input array.
 *
 * This function performs gamma correction by first normalizing the values
 * locally within a specified neighborhood around each pixel. The normalization
 * is achieved using minimum and maximum filters before applying the gamma
 * correction based on the specified gamma factor. This method allows for
 * localized adjustments of contrast and brightness in the input array.
 *
 * @param array  Input array to which the local gamma correction is applied.
 * @param gamma  Gamma factor, expected to be greater than 0. It determines the
 *               level of correction applied after local normalization. A gamma
 *               value less than 1 will lighten the local regions, while a value
 *               greater than 1 will darken them.
 * @param ir     Filter radius that defines the size of the local neighborhood
 *               used for normalization.
 * @param k      Smoothing factor used to control the degree of smoothing
 *               applied during normalization. The default value is 0.1.
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               The mask controls which parts of the array are affected by the
 *               local gamma correction.
 *
 * **Example**
 * @include ex_gamma_correction_local.cpp
 *
 * **Result**
 * @image html ex_gamma_correction_local.png
 */
void gamma_correction_local(Array &array, float gamma, int ir, float k = 0.1f);

/**
 * @brief Apply a "local" gamma correction with a mask to the input array.
 *
 * This overload applies local gamma correction to the input array using a mask.
 * The gamma correction is performed after normalizing values locally within the
 * specified neighborhood defined by the filter radius. The optional mask allows
 * for selective application of the correction.
 *
 * @param array  Input array to which the local gamma correction is applied.
 * @param gamma  Gamma factor, expected to be greater than 0. It determines the
 *               level of correction applied after local normalization. A gamma
 *               value less than 1 will lighten the local regions, while a value
 *               greater than 1 will darken them.
 * @param ir     Filter radius that defines the size of the local neighborhood
 *               used for normalization.
 * @param p_mask Filter mask, with values expected in the range [0, 1]. The mask
 *               controls which parts of the array are affected by the local
 *               gamma correction.
 * @param k      Smoothing factor used to control the degree of smoothing
 *               applied during normalization. The default value is 0.1.
 *
 * @overload
 */
void gamma_correction_local(Array       &array,
                            float        gamma,
                            int          ir,
                            const Array *p_mask,
                            float        k = 0.1f);

/**
 * @brief Applies the Kuwahara filter to an array with optional per-pixel
 * masking.
 *
 * This overloaded version of the Kuwahara filter allows applying the filter
 * with a per-pixel mask. The mask determines the blending ratio for each pixel,
 * allowing selective smoothing.
 *
 * @param array     The input/output array to apply the filter on. The array is
 *                  modified in place.
 * @param ir        The radius of the region to analyze for each pixel. Larger
 *                  values result in stronger smoothing.
 * @param mix_ratio A blending factor between the original and filtered arrays
 *                  when `p_mask` is `nullptr`. Ignored if `p_mask` is provided.
 *
 * **Example**
 * @include ex_kuwahara.cpp
 *
 * **Result**
 * @image html ex_kuwahara.png
 */
void kuwahara(Array &array, int ir, float mix_ratio = 1.f);

/**
 * @brief Applies the Kuwahara filter to an array with optional per-pixel
 * masking.
 *
 * This overloaded version of the Kuwahara filter allows applying the filter
 * with a per-pixel mask. The mask determines the blending ratio for each pixel,
 * allowing selective smoothing.
 *
 * @param array     The input/output array to apply the filter on. The array is
 *                  modified in place.
 * @param ir        The radius of the region to analyze for each pixel. Larger
 *                  values result in stronger smoothing.
 * @param p_mask    A pointer to an array representing the mask. The values in
 *                  the mask range from 0.0 to 1.0, specifying the blending
 *                  ratio for each pixel. If
 * `nullptr`, the filter is applied without masking.
 * @param mix_ratio A blending factor between the original and filtered arrays
 *                  when `p_mask` is `nullptr`. Ignored if `p_mask` is provided.
 *
 * If no mask is provided:
 * - The function behaves like the single-argument version of `kuwahara`.
 *
 * **Example**
 * @include ex_kuwahara.cpp
 *
 * **Result**
 * @image html ex_kuwahara.png
 */
void kuwahara(Array &array, int ir, const Array *p_mask, float mix_ratio = 1.f);

/**
 * @brief Applies an iterative 8-neighbor Laplacian smoothing filter on an
 * array.
 *
 * This function performs isotropic Laplacian diffusion on a 2D array using an
 * explicit finite-difference scheme. Each cell is updated according to the
 * local Laplacian computed from its 8 surrounding neighbors:
 *
 * \f[
 * A_{i,j}^{t+1} = A_{i,j}^{t} + \sigma \cdot
 * \left( \sum_{(p,q)\in N_{8}(i,j)} A_{p,q}^{t} - 8A_{i,j}^{t} \right)
 * \f]
 *
 * where \( N_{8}(i,j) \) denotes the 8 neighboring cells.
 *
 * Boundary cells are handled by clamping neighbor indices to the valid range,
 * effectively replicating edge values at the borders.
 *
 * @param array      The 2D array to be smoothed (modified in place).
 * @param sigma      Diffusion coefficient (small positive value, e.g.
 * 0.05–0.25).
 * @param iterations Number of diffusion iterations to apply.
 *
 * **Example**
 * @include ex_laplace.cpp
 *
 * **Result**
 * @image html ex_laplace.png
 */
void laplace(Array &array, float sigma = 0.125f, int iterations = 3);

/**
 * @brief Apply a low-pass Laplace filter with a mask.
 *
 * This overload applies the Laplace filter to the input array using an optional
 * mask. The mask determines which parts of the array are affected by the
 * filter. The filtering intensity and the number of iterations control the
 * extent of the smoothing effect.
 *
 * @param array      Input array to which the Laplace filter is applied.
 * @param p_mask     Optional filter mask, with values expected in the range [0,
 *                   1]. The mask specifies which parts of the array are
 *                   affected by the filter.
 * @param sigma      Filtering intensity, expected to be in the range [0, 1]. It
 *                   controls the strength of the filtering effect. A value
 *                   closer to 1 results in more smoothing. The default value is
 *                   0.2.
 * @param iterations Number of iterations to apply the filter. More iterations
 *                   will increase the smoothing effect. The default value is 3.
 *
 * @overload
 */
void laplace(Array       &array,
             const Array *p_mask,
             float        sigma = 0.125f,
             int          iterations = 3);

/**
 * @brief Apply a low-pass Laplace filter to a vector.
 *
 * This function applies a low-pass Laplace filter to a 1D vector to smooth out
 * high-frequency noise and detail. The filtering intensity and the number of
 * iterations control the extent of the smoothing effect.
 *
 * @param v          Input vector to which the Laplace filter is applied.
 * @param sigma      Filtering intensity, expected to be in the range [0, 1]. It
 *                   determines the strength of the filtering effect. A value
 *                   closer to 1 results in more smoothing.
 * @param iterations Number of iterations to apply the filter. More iterations
 *                   will increase the smoothing effect. The default value is 1.
 */
void laplace1d(std::vector<float> &v, float sigma = 0.5f, int iterations = 1);

/**
 * @brief Apply a low-pass anisotropic Laplace filter to the input array.
 *
 * This function applies an anisotropic Laplace filter to the input array.
 * Anisotropic diffusion, also known as edge-preserving smoothing, helps to
 * reduce noise while preserving edges. The effect of the filter is controlled
 * by the talus limit, filtering intensity, and the number of iterations.
 *
 * @param array      Input array to which the anisotropic Laplace filter is
 *                   applied.
 * @param talus      Talus limit for edge sensitivity. Gradients above this
 *                   value are less affected by the filtering, preserving
 *                   important edges in the data.
 * @param sigma      Filtering intensity, expected to be in the range [0, 1]. It
 *                   controls the strength of the filtering effect. A value
 *                   closer to 1 results in more smoothing. The default value is
 *                   0.2.
 * @param iterations Number of iterations to apply the filter. More iterations
 *                   will increase the smoothing effect. The default value is 3.
 *
 * **Example**
 * @include ex_laplace.cpp
 *
 * **Result**
 * @image html ex_laplace.png
 */
void laplace_edge_preserving(Array &array,
                             float  talus,
                             float  sigma = 0.2f,
                             int    iterations = 3);

/**
 * @brief Apply a low-pass anisotropic Laplace filter with a mask.
 *
 * This overload applies the anisotropic Laplace filter to the input array using
 * an optional mask. The mask determines which parts of the array are affected
 * by the filter. The effect of the filtering is controlled by the talus limit,
 * filtering intensity, and the number of iterations.
 *
 * @param array      Input array to which the anisotropic Laplace filter is
 *                   applied.
 * @param talus      Talus limit for edge sensitivity. Gradients above this
 *                   value are less affected by the filtering, preserving
 *                   important edges in the data.
 * @param p_mask     Optional filter mask, with values expected in the range [0,
 *                   1]. The mask specifies which parts of the array are
 *                   affected by the filter.
 * @param sigma      Filtering intensity, expected to be in the range [0, 1]. It
 *                   controls the strength of the filtering effect. A value
 *                   closer to 1 results in more smoothing. The default value is
 *                   0.2.
 * @param iterations Number of iterations to apply the filter. More iterations
 *                   will increase the smoothing effect. The default value is 3.
 *
 * @overload
 */
void laplace_edge_preserving(Array       &array,
                             float        talus,
                             const Array *p_mask,
                             float        sigma = 0.2f,
                             int          iterations = 3);

/**
 * @brief Apply a low-pass high-order filter to the input array.
 *
 * This function applies a high-order low-pass filter to the input array. The
 * filter can be of 5th, 7th, or 9th order, providing increasingly sharp
 * filtering effects. The filtering intensity is controlled by the `sigma`
 * parameter, which determines the strength of the filtering.
 *
 * @param array Input array to which the high-order low-pass filter is applied.
 * @param order Filter order, which can be 5, 7, or 9. Higher orders provide
 *              more refined filtering.
 * @param sigma Filtering intensity, expected to be in the range [0, 1]. It
 *              controls the strength of the filtering effect. A value closer to
 *              1 results in more smoothing. The default value is 1.0.
 *
 * **Example**
 * @include ex_low_pass_high_order.cpp
 *
 * **Result**
 * @image html ex_low_pass_high_order.png
 *
 * @see         {@link laplace}
 */
void low_pass_high_order(Array &array, int order = 9, float sigma = 1.f);

/**
 * @brief Convert array values to binary using a threshold.
 *
 * This function applies a binary filter to the input array, converting all
 * values to either 0 or 1 based on a specified threshold. Values above the
 * threshold are set to 1, while values equal to or below the threshold are set
 * to 0.
 *
 * @param array     Input array to be converted to binary.
 * @param threshold Threshold value. Any array value greater than this threshold
 *                  is set to 1, while values less than or equal to the
 *                  threshold are set to 0. The default threshold is 0.0.
 */
void make_binary(Array &array, float threshold = 0.f);

/**
 * @brief Transform the input array elevation to match the histogram of a
 * reference array.
 *
 * This function performs histogram matching on the input array to adjust its
 * elevation values so that the resulting histogram closely matches the
 * histogram of a reference array. This process is useful for aligning the
 * statistical properties of the input data with those of the reference data.
 *
 * @param array           Input array whose elevation values are to be
 *                        transformed.
 * @param array_reference Reference array whose histogram is used as the target
 *                        for matching.
 *
 * **Example**
 * @include ex_match_histogram.cpp
 *
 * **Result**
 * @image html ex_match_histogram.png
 */
void match_histogram(Array &array, const Array &array_reference);

/**
 * @brief Applies the mean shift algorithm to the input array.
 *
 * The mean shift algorithm iteratively adjusts each value in the input array by
 * averaging nearby values within a specified radius (`ir`), using either simple
 * or weighted mean computation. The process stops after the specified number of
 * iterations or if convergence criteria are met.
 *
 * @param  array          The input array to process.
 * @param  ir             The radius of the neighborhood to consider for mean
 *                        computation.
 * @param  threshold      The value threshold for considering neighboring
 *                        elements. Only elements with differences below this
 *                        threshold are included.
 * @param  iterations     The number of iterations to perform the mean shift
 *                        process.
 * @param  talus_weighted If true, uses weighted mean based on differences
 *                        between values and their neighbors. If false, uses
 *                        unweighted mean.
 * @return                A new array containing the result of the mean shift
 *                        process.
 *
 * **Example**
 * @include ex_mean_shift.cpp
 *
 * **Result**
 * @image html ex_mean_shift.png
 */
Array mean_shift(const Array &array,
                 int          ir,
                 float        talus,
                 int          iterations = 1,
                 bool         talus_weighted = true);

Array mean_shift(const Array &array,
                 int          ir,
                 float        talus,
                 const Array *p_mask,
                 int          iterations = 1,
                 bool         talus_weighted = true);

/**
 * @brief Apply a 3x3 median filter to the input array.
 *
 * This function applies a 3x3 median filter to the input array to reduce noise
 * while preserving edges. The median filter replaces each pixel with the median
 * value of the 3x3 neighborhood surrounding it. An optional mask can be used to
 * control which parts of the array are filtered.
 *
 * @param array  Input array to which the median filter is applied.
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               The mask determines which parts of the array are affected by
 *               the median filter.
 *
 * **Example**
 * @include ex_median_3x3.cpp
 *
 * **Result**
 * @image html ex_median_3x3.png
 */
void median_3x3(Array &array, const Array *p_mask);

/**
 * @brief Apply a 3x3 median filter to the input array without a mask.
 *
 * This overload applies a 3x3 median filter to the input array, affecting all
 * pixels. The median filter replaces each pixel with the median value of the
 * 3x3 neighborhood surrounding it.
 *
 * @param array Input array to which the median filter is applied.
 *
 * @overload
 */
void median_3x3(Array &array);

/**
 * @brief Computes a fast pseudo-median approximation of a local neighborhood in
 * an array.
 *
 * This function approximates the effect of a median filter by computing the
 * local minimum, maximum, and mean values within a square neighborhood of
 * radius `ir`, and averaging the results:
 *
 * \f[
 * \text{pseudo\_median}(x, y) = \frac{\text{min} + \text{max} + \text{mean}}{3}
 * \f]
 *
 * This method is computationally cheaper than a true median filter and can
 * approximate its noise-reduction and edge-preserving properties to some
 * extent.
 *
 * @param  array Input array (e.g., image or 2D signal).
 * @param  ir    Radius of the square neighborhood (kernel size will be \f$2
 *               \cdot ir + 1\f$).
 * @return       An Array containing the pseudo-median filtered output.
 *
 * @note This method works best on images with impulsive noise (e.g.,
 * salt-and-pepper), but is only an approximation and may behave differently
 * from a true median filter on complex structures or edges.
 *
 * **Example**
 * @include ex_median_pseudo.cpp
 *
 * **Result**
 * @image html ex_median_pseudo.png
 *
 * @see          local_min, local_max, local_mean
 */
Array median_pseudo(const Array &array, int ir);

/**
 * @brief Apply a displacement to the terrain along the normal direction.
 *
 * This function displaces the terrain along its normal direction, simulating
 * effects like erosion or terrain modification. The amount of displacement and
 * direction (normal or reversed) can be specified. Optionally, a mask can be
 * used to control which parts of the terrain are affected.
 *
 * @param array   Input array representing the terrain.
 * @param p_mask  Optional filter mask, with values expected in the range [0,
 *                1]. The mask specifies which parts of the array are affected
 *                by the displacement.
 * @param amount  Amount of displacement to apply. The default value is 0.1.
 * @param ir      Pre-filtering radius. The radius used to smooth the terrain
 *                before applying displacement. The default value is 0.
 * @param reverse If true, the displacement direction is reversed. The default
 *                value is false.
 *
 * **Example**
 * @include ex_normal_displacement.cpp
 *
 * **Result**
 * @image html ex_normal_displacement.png
 */
void normal_displacement(Array &array,
                         float  amount = 0.1f,
                         int    ir = 0,
                         bool   reverse = false);

void normal_displacement(Array       &array,
                         const Array *p_mask,
                         float        amount = 0.1f,
                         int          ir = 0,
                         bool         reverse = false); ///< @overload

/**
 * @brief Apply a plateau-shape filter to the input array.
 *
 * This function applies a plateau-shape filter to the input array, modifying
 * the terrain to create flat, plateau-like regions. The filter's radius and
 * flatness can be adjusted to control the extent and shape of the plateaus. An
 * optional mask can be used to control which parts of the array are affected.
 *
 * @param array  Input array to which the plateau-shape filter is applied.
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               The mask specifies which parts of the array are affected by the
 *               filter.
 * @param ir     Plateau radius. The size of the area over which the plateau
 *               effect is applied.
 * @param factor Gain factor that controls the flatness of the plateau. The
 *               higher the factor, the flatter the resulting plateau. The
 *               default value is 1.0.
 *
 * **Example**
 * @include ex_plateau.cpp
 *
 * **Result**
 * @image html ex_plateau.png
 */
void plateau(Array &array, const Array *p_mask, int ir, float factor);

void plateau(Array &array, int ir, float factor); ///< @overload

/**
 * @brief Quantize array values into discrete levels.
 *
 * Values are clamped to the range [vmin, vmax] and remapped to `nlevels`
 * uniformly distributed discrete values.
 *
 * @param  array   Input array.
 * @param  nlevels Number of quantization levels.
 * @param  vmin    Minimum clamp value.
 * @param  vmax    Maximum clamp value.
 * @return         Quantized array.
 *
 * **Example**
 * @include ex_dequantize.cpp
 *
 * **Result**
 * @image html ex_dequantize.png
 */
Array quantize(const Array &array, int nlevels, float vmin, float vmax);

/**
 * @brief Quantize array values using the array min/max range.
 *
 * @param  array   Input array.
 * @param  nlevels Number of quantization levels.
 * @return         Quantized array.
 *
 * **Example**
 * @include ex_dequantize.cpp
 *
 * **Result**
 * @image html ex_dequantize.png
 */
Array quantize(const Array &array, int nlevels);

/**
 * @brief Transform heightmap to give a "billow" like appearance.
 *
 * This function modifies the heightmap to produce a "billow" effect. The
 * transformation is applied based on a reference elevation (`vref`) where the
 * elevation is folded. The smoothing coefficient (`k`) controls the degree of
 * smoothing applied during the transformation.
 *
 * @param array Input array representing the heightmap to be transformed.
 * @param vref  Reference elevation where the heightmap values are folded. This
 *              value determines the baseline around which the billow effect is
 *              centered.
 * @param k     Smoothing coefficient that influences the degree of smoothness
 *              in the billow transformation. Higher values result in more
 *              smoothing.
 *
 * **Example**
 * @include ex_recast.cpp
 *
 * **Result**
 * @image html ex_recast.png
 */
void recast_billow(Array &array, float vref, float k);

/**
 * @brief Transform heightmap to give a "canyon" like appearance.
 *
 * This function modifies the heightmap to create a "canyon" effect. The
 * transformation is based on a canyon top elevation (`vcut`). The effect can be
 * adjusted with a gamma factor to control the intensity of the canyon effect.
 * Optional parameters include a filter mask for selective application and a
 * noise array for additional variation.
 *
 * @param array   Input array representing the heightmap to be transformed.
 * @param vcut    Canyon top elevation. This value defines the threshold for
 *                canyon formation.
 * @param gamma   Gamma factor (> 0) that adjusts the intensity of the canyon
 *                effect. A higher value increases the effect. The default value
 *                is 4.0.
 * @param p_mask  Optional filter mask, with values expected in the range [0,
 *                1]. The mask specifies which parts of the heightmap are
 *                affected by the transformation.
 * @param p_noise Optional noise array used to introduce variation in the canyon
 *                effect. If not provided, no additional noise is applied.
 *
 * **Example**
 * @include ex_recast.cpp
 *
 * **Result**
 * @image html ex_recast.png
 */
void recast_canyon(Array &array, const Array &vcut, float gamma = 4.f);

void recast_canyon(Array       &array,
                   const Array &vcut,
                   const Array *p_mask,
                   float        gamma = 4.f); ///< @overload

void recast_canyon(Array       &array,
                   float        vcut,
                   const Array *p_mask,
                   float        gamma = 4.f,
                   const Array *p_noise = nullptr); ///< @overload

void recast_canyon(Array       &array,
                   float        vcut,
                   float        gamma = 4.f,
                   const Array *p_noise = nullptr); ///< @overload

/**
 * @brief Transform heightmap to add cliffs where gradients are steep enough.
 *
 * This function modifies the heightmap to introduce cliffs based on a reference
 * talus angle. Cliffs are added where the gradient exceeds the specified talus
 * angle. The amplitude of the cliffs and the gain factor for steepness are
 * adjustable. Additionally, a filter mask can be applied to control which parts
 * of the heightmap are affected.
 *
 * @param array     Input array representing the heightmap to be modified.
 * @param talus     Reference talus angle. This angle determines the threshold
 *                  above which cliffs are formed.
 * @param ir        Filter radius used to smooth the heightmap before applying
 *                  the cliff effect.
 * @param amplitude Amplitude of the cliffs. This value controls the height of
 *                  the cliffs.
 * @param gain      Gain factor for the gain filter, influencing the steepness
 *                  of the cliffs. Higher values result in steeper cliffs. The
 *                  default value is 2.0.
 * @param p_mask    Optional filter mask, with values expected in the range [0,
 *                  1]. The mask specifies which parts of the heightmap are
 *                  affected by the cliff transformation.
 *
 * **Example**
 * @include ex_recast.cpp
 *
 * **Result**
 * @image html ex_recast.png
 */
void recast_cliff(Array &array,
                  float  talus,
                  int    ir,
                  float  amplitude,
                  float  gain = 2.f);

void recast_cliff(Array       &array,
                  float        talus,
                  int          ir,
                  float        amplitude,
                  const Array *p_mask,
                  float        gain = 2.f); ///< @overload

/**
 * @brief Transform heightmap to add directional cliffs where gradients are
 * steep enough.
 *
 * This function modifies the heightmap to introduce cliffs in a specific
 * direction, based on a reference talus angle. The cliffs are added where the
 * gradient exceeds the specified talus angle, with the direction of the cliffs
 * controlled by the specified angle. The amplitude of the cliffs and the gain
 * factor for steepness are also adjustable. A filter mask can be used to
 * specify which parts of the heightmap are affected.
 *
 * @param array     Input array representing the heightmap to be modified.
 * @param talus     Reference talus angle. This angle determines the threshold
 *                  above which cliffs are formed.
 * @param ir        Filter radius used to smooth the heightmap before applying
 *                  the cliff effect.
 * @param amplitude Amplitude of the cliffs. This value controls the height of
 *                  the cliffs.
 * @param angle     Angle (in degrees) determining the direction of the cliffs.
 * @param gain      Gain factor for the gain filter, influencing the steepness
 *                  of the cliffs. Higher values result in steeper cliffs. The
 *                  default value is 2.0.
 * @param p_mask    Optional filter mask, with values expected in the range [0,
 *                  1]. The mask specifies which parts of the heightmap are
 *                  affected by the cliff transformation.
 *
 * **Example**
 * @include ex_recast.cpp
 *
 * **Result**
 * @image html ex_recast.png
 */
void recast_cliff_directional(Array &array,
                              float  talus,
                              int    ir,
                              float  amplitude,
                              float  angle,
                              float  gain = 2.f); ///< @overload

void recast_cliff_directional(Array       &array,
                              float        talus,
                              int          ir,
                              float        amplitude,
                              float        angle,
                              const Array *p_mask,
                              float        gain = 2.f); ///< @overload

void recast_cracks(Array &array,
                   float  cut_min = 0.05f,
                   float  cut_max = 0.5f,
                   float  k_smoothing = 0.01f,
                   float  vmin = 0.f,
                   float  vmax = -1.f);

/**
 * @brief Applies an escarpment effect to the given 2D array, modifying its
 * values based on cumulative displacement with optional directional and
 * transpositional transformations.
 *
 * This function calculates cumulative displacement along the x-axis based on
 * relative elevation differences between adjacent cells. The displacement is
 * scaled and optionally reversed, then used to apply a warping effect to the
 * array, simulating an escarpment feature. An optional global scaling factor
 * can further adjust the displacement effect intensity.
 *
 * @param array            Reference to the 2D array where the escarpment effect
 *                         will be applied.
 * @param ir               Radius for the smoothing kernel used on the
 *                         displacement, controlling the smoothness of the
 *                         effect. Larger values result in smoother transitions.
 * @param ratio            The ratio influencing displacement; values > 1.0
 *                         increase displacement sensitivity to height
 *                         differences.
 * @param scale            The scaling factor for the cumulative displacement,
 *                         affecting the intensity of the effect.
 * @param reverse          If true, reverses the direction of the displacement
 *                         effect.
 * @param transpose_effect If true, transposes the array before and after
 *                         applying the effect.
 * @param global_scaling   An additional scaling factor for the displacement; if
 *                         set to 0, a default value is computed based on the
 *                         array's range and size. Higher values increase the
 *                         overall effect.
 *
 * **Example**
 * @include ex_recast.cpp
 *
 * **Result**
 * @image html ex_recast.png
 */
void recast_escarpment(Array &array,
                       int    ir = 16,
                       float  ratio = 0.1f,
                       float  scale = 1.f,
                       bool   reverse = false,
                       bool   transpose_effect = false,
                       float  global_scaling = 0.f);

/**
 * @brief Applies an escarpment effect to the given 2D array, with an optional
 * mask to blend the effect.
 *
 * This overload allows for a blending mask that controls the intensity of the
 * escarpment effect at each point in the array. The mask values range between 0
 * and 1, where values closer to 1 fully apply the effect, and values closer to
 * 0 reduce it. An optional global scaling factor can further adjust the
 * displacement effect intensity.
 *
 * @param array            Reference to the 2D array where the escarpment effect
 *                         will be applied.
 * @param p_mask           Pointer to an optional mask array for blending the
 *                         effect, where values range from 0 to 1. A nullptr
 *                         applies the effect fully without blending.
 * @param ir               Radius for the smoothing kernel used on the
 *                         displacement, controlling the smoothness of the
 *                         effect. Larger values result in smoother transitions.
 * @param ratio            The ratio influencing displacement; values > 1.0
 *                         increase displacement sensitivity to height
 *                         differences.
 * @param scale            The scaling factor for the cumulative displacement,
 *                         affecting the intensity of the effect.
 * @param reverse          If true, reverses the direction of the displacement
 *                         effect.
 * @param transpose_effect If true, transposes the array before and after
 *                         applying the effect.
 * @param global_scaling   An additional scaling factor for the displacement; if
 *                         set to 0, a default value is computed based on the
 *                         array's range and size. Higher values increase the
 *                         overall effect.
 *
 * **Example**
 * @include ex_recast.cpp
 *
 * **Result**
 * @image html ex_recast.png
 */
void recast_escarpment(Array       &array,
                       const Array *p_mask,
                       int          ir = 16,
                       float        ratio = 0.1f,
                       float        scale = 1.f,
                       bool         reverse = false,
                       bool         transpose_effect = false,
                       float        global_scaling = 0.f);

/**
 * @brief Transform heightmap to give a "peak" like appearance.
 *
 * This function modifies the heightmap to create a "peak" effect, where the
 * elevations are adjusted to emphasize peak-like features. The transformation
 * is applied based on a filter radius (`ir`) and involves gamma correction and
 * smoothing parameters. The effect can be selectively applied using a filter
 * mask.
 *
 * @warning Array values are expected to be in the range [0, 1].
 *
 * @param array  Input array representing the heightmap to be transformed.
 * @param ir     Filter radius used to smooth the heightmap before applying the
 *               peak effect.
 * @param gamma  Gamma factor (> 0) that adjusts the intensity of the peak
 *               effect. A higher value increases the effect. The default value
 *               is 2.0. For details on gamma correction, see {@link
 *               gamma_correction}.
 * @param k      Smoothing parameter (> 0) that controls the degree of smoothing
 *               applied during the transformation. Higher values result in more
 *               smoothing. The default value is 0.1.
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               The mask specifies which parts of the heightmap are affected by
 *               the peak transformation.
 *
 * **Example**
 * @include ex_recast.cpp
 *
 * **Result**
 * @image html ex_recast.png
 */
void recast_peak(Array &array, int ir, float gamma = 2.f, float k = 0.1f);

void recast_peak(Array       &array,
                 int          ir,
                 const Array *p_mask,
                 float        gamma = 2.f,
                 float        k = 0.1f); ///< @overload

/**
 * @brief Transform heightmap by adding "rock-like" features at higher slopes.
 *
 * This function modifies the heightmap to introduce rock-like features,
 * particularly at locations with higher slopes. The transformation is
 * controlled by parameters defining talus limit, filter radius, amplitude, and
 * noise characteristics. The rock features are generated using a specified
 * random seed and can be customized with optional noise and gamma correction. A
 * filter mask can be applied to control which parts of the heightmap are
 * affected.
 *
 * @param array     Input array representing the heightmap to be transformed.
 * @param talus     Talus limit that determines the threshold above which
 *                  rock-like features are introduced based on slope steepness.
 * @param ir        Filter radius used to smooth the heightmap and define the
 *                  area of influence for adding rock features.
 * @param amplitude Amplitude of the rock features. This value controls the
 *                  intensity of the rock-like appearance.
 * @param seed      Random number seed used to generate noise patterns for the
 *                  rock features. Different seeds produce different rock
 *                  patterns.
 * @param kw        Noise wavenumber with respect to a unit domain, influencing
 *                  the frequency of rock features. Higher values result in more
 *                  detailed rock textures.
 * @param gamma     Gamma correction coefficient applied to the rock features.
 *                  This parameter adjusts the contrast of the rock features.
 *                  The default value is 0.5.
 * @param p_noise   Optional reference to an input noise array used for
 *                  generating rock features. If provided, it overrides the
 *                  default noise generator.
 * @param p_mask    Optional filter mask, with values expected in the range [0,
 *                  1]. The mask specifies which parts of the heightmap are
 *                  affected by the rock feature transformation.
 * @param bbox      Bounding box defining the region of interest within the
 *                  heightmap. This parameter specifies the coordinates of the
 *                  bounding box in the format {xmin, xmax, ymin, ymax}. The
 *                  default bounding box covers the entire domain.
 *
 * **Example**
 * @include ex_recast.cpp
 *
 * **Result**
 * @image html ex_recast.png
 */
void recast_rocky_slopes(Array        &array,
                         float         talus,
                         int           ir,
                         float         amplitude,
                         std::uint32_t seed,
                         float         kw,
                         float         gamma = 0.5f,
                         const Array  *p_noise = nullptr,
                         glm::vec4     bbox = {0.f, 1.f, 0.f, 1.f});

void recast_rocky_slopes(Array        &array,
                         float         talus,
                         int           ir,
                         float         amplitude,
                         std::uint32_t seed,
                         float         kw,
                         const Array  *p_mask,
                         float         gamma = 0.5f,
                         const Array  *p_noise = nullptr,
                         glm::vec4 bbox = {0.f, 1.f, 0.f, 1.f}); ///< @overload

/**
 * @brief Transform heightmap to give a "cliff" like appearance.
 *
 * This function modifies the heightmap to create a "cliff" effect by folding
 * the elevations around a reference height (`vref`). The transformation is
 * influenced by a smoothing coefficient (`k`) and can be selectively applied
 * using a filter mask.
 *
 * @param array  Input array representing the heightmap to be transformed.
 * @param vref   Reference elevation where the heightmap is folded to create the
 *               cliff-like appearance. Elevations near this value will be
 *               prominently altered.
 * @param k      Smoothing coefficient that controls the degree of smoothing
 *               applied during the transformation. A higher value results in
 *               smoother cliffs. The default value is not specified in this
 *               function; it must be provided by the user.
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               The mask specifies which parts of the heightmap are affected by
 *               the cliff transformation. If not provided, the transformation
 *               is applied to the entire heightmap.
 *
 * **Example**
 * @include ex_recast.cpp
 *
 * **Result**
 * @image html ex_recast.png
 */
void recast_sag(Array &array, float vref, float k);

void recast_sag(Array       &array,
                float        vref,
                float        k,
                const Array *p_mask); ///< @overload

/**
 * @brief Apply a curve adjustment filter to the array.
 *
 * This function applies a curve adjustment to the heightmap by mapping input
 * values to output values according to a given correction curve. The curve is
 * monotonically interpolated, and any values outside the defined input range
 * are clipped. The filter can be selectively applied using a mask.
 *
 * @param array  Input array representing the heightmap to be adjusted.
 * @param t      Vector of input values for the correction curve. These values
 *               define the x-coordinates (input values) of the curve.
 * @param v      Vector of output values for the correction curve. These values
 *               define the y-coordinates (output values) of the curve
 *               corresponding to the input values in `t`.
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               The mask specifies which parts of the heightmap are affected by
 *               the curve adjustment. If not provided, the curve adjustment is
 *               applied to the entire heightmap.
 *
 * **Example**
 * @include ex_recurve.cpp
 *
 * **Result**
 * @image html ex_recurve.png
 */
void recurve(Array                    &array,
             const std::vector<float> &t,
             const std::vector<float> &v);

void recurve(Array                    &array,
             const std::vector<float> &t,
             const std::vector<float> &v,
             const Array              *p_mask); ///< @overload

/**
 * @brief Apply a curve adjustment filter using a "bumpy exponential-shape"
 * curve.
 *
 * This function adjusts the heightmap using a "bumpy exponential-shape" curve,
 * where the curve has an exponential decay defined by the parameter `tau`. The
 * curve transformation is applied to the entire heightmap or selectively using
 * a filter mask.
 *
 * @warning Array values are expected to be in the range [0, 1].
 *
 * @param array  Input array representing the heightmap to be adjusted.
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               The mask specifies which parts of the heightmap are affected by
 *               the curve adjustment. If not provided, the adjustment is
 *               applied to the entire heightmap.
 * @param tau    Exponential decay parameter that defines the shape of the
 *               "bumpy"
 * curve. Higher values of `tau` result in a sharper curve. The default value is
 * 0.5.
 *
 * **Example**
 * @include ex_recurve_xxx.cpp
 *
 * **Result**
 * @image html ex_recurve_xxx.png
 */
void recurve_bexp(Array &array, float tau = 0.5f);

void recurve_bexp(Array       &array,
                  const Array *p_mask,
                  float        tau = 0.5f); ///< @overload

/**
 * @brief Apply a curve adjustment filter using a "sharp exponential-shape"
 * curve.
 *
 * This function adjusts the heightmap using a "sharp exponential-shape" curve,
 * where the curve has an exponential decay defined by the parameter `tau`. The
 * curve transformation can be applied to the entire heightmap or selectively
 * using a filter mask.
 *
 * @warning Array values are expected to be in the range [0, 1].
 *
 * @param array  Input array representing the heightmap to be adjusted.
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               The mask specifies which parts of the heightmap are affected by
 *               the curve adjustment. If not provided, the adjustment is
 *               applied to the entire heightmap.
 * @param tau    Exponential decay parameter that defines the shape of the
 *               "sharp"
 * curve. Higher values of `tau` result in a steeper curve. The default value is
 * 0.5.
 *
 * **Example**
 * @include ex_recurve_xxx.cpp
 *
 * **Result**
 * @image html ex_recurve_xxx.png
 */
void recurve_exp(Array &array, float tau = 0.5f);

void recurve_exp(Array       &array,
                 const Array *p_mask,
                 float        tau = 0.5f); ///< @overload

/**
 * @brief Apply a curve adjustment filter using Kumaraswamy's cumulative
 * distribution function (CDF).
 *
 * This function transforms the heightmap using Kumaraswamy's CDF, which is
 * defined by two shape parameters, `a` and `b`. The transformation applies a
 * curve adjustment that can shape the data closer to 0 or 1 depending on the
 * values of `a` and `b`. The adjustment can be applied to the entire heightmap
 * or selectively using a filter mask.
 *
 * @param array  Input array representing the heightmap to be adjusted.
 * @param a      Parameter 'a' of Kumaraswamy's CDF, which drives the curve
 *               shape towards 0.
 * @param b      Parameter 'b' of Kumaraswamy's CDF, which drives the curve
 *               shape towards 1.
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               If provided, the adjustment is applied according to this mask.
 *               If not provided, the entire heightmap is adjusted.
 *
 * **Example**
 * @include ex_recurve_xxx.cpp
 *
 * **Result**
 * @image html ex_recurve_xxx.png
 */
void recurve_kura(Array &array, float a, float b);

void recurve_kura(Array       &array,
                  float        a,
                  float        b,
                  const Array *p_mask); ///< @overload

/**
 * @brief Apply a curve adjustment filter using a smooth "S-shape" curve.
 *
 * This function adjusts the heightmap using a smooth "S-shape" curve, which
 * transforms the input values in a manner that creates an S-shaped
 * distribution. The transformation can be applied to the entire heightmap or
 * selectively using a filter mask.
 *
 * @warning Array values are expected to be in the range [0, 1].
 *
 * @param array  Input array representing the heightmap to be adjusted.
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               If provided, the adjustment is applied according to this mask.
 *               If not provided, the entire heightmap is adjusted.
 *
 * **Example**
 * @include ex_recurve_xxx.cpp
 *
 * **Result**
 * @image html ex_recurve_xxx.png
 */
void recurve_s(Array &array);

void recurve_s(Array &array, const Array *p_mask); ///< @overload

/**
 * @brief Apply a curve adjustment filter using an nth-order smoothstep curve.
 *
 * This function applies a curve adjustment to the heightmap using a smoothstep
 * function of order `n`. The smoothstep function provides a smooth transition
 * between 0 and 1, where `n` determines the smoothness of the transition. The
 * adjustment can be applied to the entire heightmap or selectively using a
 * filter mask.
 *
 * @warning Array values are expected to be in the range [0, 1].
 *
 * @param array  Input array representing the heightmap to be adjusted.
 * @param n      Smoothstep order, which determines the degree of smoothness for
 *               the transition. `n` should be a non-negative value (in [0, ∞[).
 * @param p_mask Optional filter mask, with values expected in the range [0, 1].
 *               If provided, the adjustment is applied according to this mask.
 *               If not provided, the entire heightmap is adjusted.
 *
 * **Example**
 * @include ex_recurve_xxx.cpp
 *
 * **Result**
 * @image html ex_recurve_xxx.png
 */
void recurve_smoothstep_rational(Array &array, float n);

void recurve_smoothstep_rational(Array       &array,
                                 float        n,
                                 const Array *p_mask); ///< @overload

/**
 * @brief Applies a smooth reversal of values above a given threshold.
 *
 * This function modifies the elements of an array by "reversing" values that
 * exceed a given threshold value, optionally blending the effect over a
 * transition range.
 *
 * For values above the threshold:
 *  - The difference between the value and the threshold is scaled down
 * proportionally to `scaling` and modulated by a smooth transition (using
 * smoothstep3) based on `transition_extent`.
 *
 * For values below the threshold:
 *  - The values are linearly blended (`lerp`) with the threshold value using a
 * smooth transition factor, gradually reducing the reversal effect.
 *
 * @param array             Array to modify in place.
 * @param threshold         Array containing the per-element threshold values.
 * @param scaling           Multiplier applied to the difference above the
 *                          threshold.
 * @param transition_extent Width of the smooth transition zone for blending.
 *
 * **Example**
 * @include ex_reverse_above_theshold.cpp
 *
 * **Result**
 * @image html ex_reverse_above_theshold.png
 */
void reverse_above_theshold(Array       &array,
                            const Array &threshold,
                            float        scaling = 1.f,
                            float        transition_extent = 0.f);

void reverse_above_theshold(Array &array,
                            float  threshold,
                            float  scaling = 1.f,
                            float  transition_extent = 0.f); ///< @overload

void reverse_above_theshold(Array       &array,
                            const Array &threshold,
                            const Array *p_mask,
                            float        scaling = 1.f,
                            float transition_extent = 0.f); ///< @overload

void reverse_above_theshold(Array       &array,
                            float        threshold,
                            const Array *p_mask,
                            float        scaling = 1.f,
                            float transition_extent = 0.f); ///< @overload

/**
 * @brief Saturate the array values based on the input interval [vmin, vmax]
 * (the output amplitude is not modified).
 *
 * This function adjusts the heightmap values to fit within the specified range
 * [vmin, vmax]. Values outside the range [from_min, from_max] are clamped to
 * [vmin, vmax]. The output amplitude of the heightmap remains unchanged.
 *
 * @param array    Input array representing the heightmap to be saturated.
 * @param vmin     The lower bound of the range to remap to.
 * @param vmax     The upper bound of the range to remap to.
 * @param from_min The lower bound of the range to remap from.
 * @param from_max The upper bound of the range to remap from.
 *
 * **Example**
 * @include ex_saturate.cpp
 *
 * **Result**
 * @image html ex_saturate.png
 *
 * See unit tests: @ref test_saturate.cpp
 */
void saturate(Array &array,
              float  vmin,
              float  vmax,
              float  from_min,
              float  from_max,
              float  k = 0.f);

void saturate(Array &array, float vmin, float vmax, float k = 0.f); ///<
// @overload

/**
 * @brief Saturate values based on percentile range.
 *
 * Computes low/high bounds from percentiles and applies saturation scaling.
 *
 * @param array           Input array (modified in place).
 * @param percentile_low  Lower percentile [0–1].
 * @param percentile_high Upper percentile [0–1].
 * @param k               Saturation strength.
 *
 * **Example**
 * @include ex_saturate.cpp
 *
 * **Result**
 * @image html ex_saturate.png
 *
 * See unit tests: @ref test_saturate.cpp
 */
void saturate_percentile(Array &array,
                         float  percentile_low,
                         float  percentile_high,
                         float  k = 0.f);

/**
 * @brief Apply a sharpening filter based on the Laplace operator.
 *
 * This function applies a sharpening filter to the input array using the
 * Laplace operator. The sharpening effect can be controlled by adjusting the
 * `ratio` parameter, which determines the balance between the sharpened and
 * non-sharpened output. A `ratio` of 1 applies full sharpening, while a `ratio`
 * of 0 applies no sharpening.
 *
 * @param array  Input array representing the heightmap or image to be
 *               sharpened.
 * @param p_mask Optional filter mask, expected in the range [0, 1]. If
 *               provided, sharpening is applied according to this mask. If not
 *               provided, the entire array is sharpened.
 * @param ratio  Ratio of the sharpening effect. A value of 1 applies full
 *               sharpening, while a value of 0 applies no sharpening. Default
 *               is 1.
 *
 * **Example**
 * @include ex_sharpen.cpp
 *
 * **Result**
 * @image html ex_sharpen.png
 */
void sharpen(Array &array, float ratio = 1.f);

void sharpen(Array       &array,
             const Array *p_mask,
             float        ratio = 1.f); ///< @overload

/**
 * @brief Apply a sharpening filter based on a smooth cone filter.
 *
 * This function applies a sharpening filter to the input array using a smooth
 * cone filter. The sharpening effect can be controlled by the `scale`
 * parameter, and the size of the smoothing operation is determined by the `ir`
 * parameter. The optional filter mask allows selective application of the
 * sharpening effect.
 *
 * @param array  Input array representing the heightmap or image to be
 *               sharpened.
 * @param p_mask Optional filter mask, expected in the range [0, 1]. If
 *               provided, sharpening is applied according to this mask. If not
 *               provided, the entire array is sharpened.
 * @param ir     Filter radius, which determines the size of the smoothing
 *               operation.
 * @param scale  Sharpening scale. Adjusts the intensity of the sharpening
 *               effect. Default is 0.5.
 *
 * **Example**
 * @include ex_sharpen_cone.cpp
 *
 * **Result**
 * @image html ex_sharpen_cone.png
 */
void sharpen_cone(Array &array, int ir, float intensity = 0.5f);

void sharpen_cone(Array       &array,
                  const Array *p_mask,
                  int          ir,
                  float        scale = 0.5f); ///< @overload

/**
 * @brief Apply shrinking, or "deflating", to emphasize the ridges in the
 * heightmap.
 *
 * This function applies a shrinking filter to the input array, which
 * accentuates the ridges by shrinking or deflating the terrain features. The
 * filter radius `ir` determines the extent of the deflation. You can also apply
 * the shrinking effect selectively using an optional filter mask. Additionally,
 * an alternative method using a custom kernel is available.
 *
 * @param array  Input array representing the heightmap to be processed.
 * @param ir     Filter radius, which controls the extent of the shrinking
 *               effect.
 * @param p_mask Optional filter mask, expected in the range [0, 1]. If
 *               provided, shrinking is applied according to this mask. If not
 *               provided, the entire array is processed.
 * @param kernel Optional custom kernel to be used for the shrinking operation.
 *               If provided, this kernel will override the default filter
 *               radius.
 *
 * **Example**
 * @include ex_expand.cpp
 *
 * **Result**
 * @image html ex_expand.png
 *
 * @see          {@link ex_expand}
 */
void shrink(Array &array, int ir, int iterations = 1);
void shrink(Array       &array,
            int          ir,
            const Array *p_mask,
            int          iterations = 1); ///< @overload
void shrink(Array       &array,
            const Array &kernel,
            int          iterations = 1); ///< @overload
void shrink(Array       &array,
            const Array &kernel,
            const Array *p_mask,
            int          iterations = 1); ///< @overload

/**
 * @brief Apply directional shrinking, or "deflating", to emphasize the ridges
 * in the terrain.
 *
 * This function applies a directional shrinking filter to the input array,
 * enhancing the ridges by shrinking or deflating the terrain features with
 * respect to a specified direction. Parameters include filter radius `ir`,
 * direction angle, aspect ratio, and anisotropy to control the shape and
 * orientation of the filter. An optional filter mask can be used for selective
 * application.
 *
 * @param array        Input array representing the heightmap to be processed.
 * @param ir           Filter radius, which controls the extent of the shrinking
 *                     effect.
 * @param angle        Angle (in degrees) of the directional filter.
 * @param aspect_ratio Pulse aspect ratio, which influences the shape of the
 *                     filter.
 * @param anisotropy   Pulse width ratio between upstream and downstream sides.
 *                     Default is 1.
 * @param p_mask       Optional filter mask, expected in the range [0, 1]. If
 *                     provided, shrinking is applied according to this mask. If
 *                     not provided, the entire array is processed.
 *
 * **Example**
 * @include ex_expand_directional.cpp
 *
 * **Result**
 * @image html ex_expand_directional.png
 */
void shrink_directional(Array       &array,
                        int          ir,
                        float        angle,
                        float        aspect_ratio,
                        float        anisotropy = 1.f,
                        const Array *p_mask = nullptr);

/**
 * @brief Apply a convolution filter with a cone kernel to smooth the array.
 *
 * This function performs convolution on the input array using a cone-shaped
 * kernel. The cone radius `ir` determines the extent of the smoothing effect.
 * An optional filter mask can be used for selective smoothing.
 *
 * @param array  Input array representing the heightmap or image to be smoothed.
 * @param ir     Cone radius, which controls the extent of the smoothing effect.
 *               The cone's half-width is half this radius.
 * @param p_mask Optional filter mask, expected in the range [0, 1]. If
 *               provided, smoothing is applied according to this mask. If not
 *               provided, the entire array is processed.
 *
 * **Example**
 * @include ex_smooth_cone.cpp
 *
 * **Result**
 * @image html ex_smooth_cone.png
 */
void smooth_cone(Array &array, int ir);
void smooth_cone(Array &array, int ir, const Array *p_mask); ///< @overload

/**
 * @brief Apply filtering to the array using convolution with a cubic pulse
 * kernel.
 *
 * This function performs convolution on the input array using a cubic pulse
 * kernel. It serves as an alternative to Gaussian smoothing, offering a more
 * compact support with potentially faster computations. For direct comparison
 * with Gaussian smoothing, the pulse radius `ir` should be set to twice the
 * desired Gaussian half-width.
 *
 * @param array  Input array representing the heightmap or image to be smoothed.
 * @param ir     Pulse radius, where the half-width of the cubic pulse kernel is
 *               half of this value.
 * @param p_mask Optional filter mask, expected in the range [0, 1]. If
 *               provided, filtering is applied according to this mask. If not
 *               provided, the entire array is processed.
 *
 * **Example**
 * @include ex_smooth_cpulse.cpp
 *
 * **Result**
 * @image html ex_smooth_cpulse.png
 *
 * See unit tests: @ref test_smooth_cpulse.cpp
 */
void smooth_cpulse(Array &array, int ir);
void smooth_cpulse(Array &array, int ir, const Array *p_mask); ///< @overload

/**
 * @brief Smooths an array while attenuating edges using a gradient-based pulse.
 *
 * This function first computes the gradient norm of the input array, applies a
 * sigmoid filter to highlight or attenuate regions based on the local gradient
 * magnitude, and then performs a constrained pulse smoothing pass that respects
 * those gradient-based weights. This is useful for smoothing heightmaps or
 * scalar fields while preserving significant edges such as ridges or cliffs.
 *
 * @param array       Reference to the input array to be smoothed. The array is
 *                    modified in-place.
 * @param talus       Threshold value used in the sigmoid function to control
 *                    where edge attenuation starts. Lower values affect more
 *                    regions.
 * @param talus_width Width (steepness) of the sigmoid transition zone around
 *                    the talus threshold.
 * @param ir          Radius of the smoothing kernel used in the constrained
 *                    pulse smoothing step.
 */
void smooth_cpulse_edge_removing(Array &array,
                                 float  talus,
                                 float  talus_width,
                                 int    ir);

/**
 * @brief Applies a smoothing average filter to the given 2D array in both
 * dimensions.
 *
 * This function creates a smoothing kernel of size \(2 \times \text{ir} + 1\)
 * with uniform weights, then applies a 1D convolution along both the i (rows)
 * and j (columns) dimensions of the array to achieve a 2D smoothing effect.
 *
 * @param array Reference to the 2D array to be smoothed.
 * @param ir    Radius of the smoothing kernel, determining its size as \(2
 *              \times
 * \text{ir} + 1\). Larger values produce more smoothing.
 */
void smooth_flat(Array &array, int ir);

/**
 * @brief Apply Gaussian filtering to the array.
 *
 * This function performs convolution on the input array using a Gaussian
 * kernel, which smooths the array by averaging values within a
 * Gaussian-weighted neighborhood. The `ir` parameter specifies the Gaussian
 * half-width, which influences the extent of the smoothing effect. An optional
 * filter mask can be provided for selective filtering.
 *
 * @param array  Input array representing the heightmap or image to be smoothed.
 * @param ir     Gaussian half-width, which determines the extent of the
 *               smoothing effect.
 * @param p_mask Optional filter mask, expected in the range [0, 1]. If
 *               provided, filtering is applied according to this mask. If not
 *               provided, the entire array is processed.
 *
 * **Example**
 * @include ex_smooth_gaussian.cpp
 *
 * **Result**
 * @image html ex_smooth_gaussian.png
 */
void smooth_gaussian(Array &array, int ir);
void smooth_gaussian(Array &array, int ir, const Array *p_mask); ///< @overload
/**
 * @brief Apply cubic pulse smoothing to fill lower flat regions while
 * preserving some sharpness.
 *
 * This function applies a cubic pulse smoothing technique to the input array,
 * specifically designed to fill in lower flat regions without significantly
 * impacting sharp features. It provides a more computationally efficient
 * alternative to the `thermal_auto_bedrock` method. The `p_deposition_map`
 * output field captures the amount of smoothing applied.
 *
 * @param array            Input array representing the heightmap or image to be
 *                         processed.
 * @param ir               Pulse radius, which determines the extent of the
 *                         smoothing effect.
 * @param k                Transition smoothing parameter in the range [0, 1].
 *                         It controls the balance between smoothing and
 *                         preserving sharpness.
 * @param p_mask           Optional filter mask, expected in the range [0, 1].
 *                         If provided, smoothing is applied according to this
 *                         mask. If not provided, the entire array is processed.
 * @param p_deposition_map [out] Optional reference to the deposition map. This
 *                         output field records the amount of deposition or
 *                         smoothing applied.
 *
 * **Example**
 * @include ex_smooth_fill.cpp
 *
 * **Result**
 * @image html ex_smooth_fill.png
 *
 * @see                    {@link smooth_cpulse}, {@link thermal_auto_bedrock}
 */
void smooth_fill(Array &array,
                 int    ir,
                 float  k = 0.1f,
                 Array *p_deposition_map = nullptr);

void smooth_fill(Array       &array,
                 int          ir,
                 const Array *p_mask,
                 float        k = 0.1f,
                 Array       *p_deposition_map = nullptr); ///< @overload

/**
 * @brief Apply smoothing to fill holes (elliptic concave surfaces).
 *
 * This function applies a smoothing filter to the input array to fill in holes
 * or concave surfaces. The goal is to smooth out regions that are concave,
 * providing a more uniform appearance.
 *
 * @param array  Input array representing the heightmap or image to be
 *               processed.
 * @param ir     Filter radius, which controls the extent of the smoothing
 *               effect.
 * @param p_mask Optional filter mask, expected in the range [0, 1]. If
 *               provided, smoothing is applied according to this mask. If not
 *               provided, the entire array is processed.
 *
 * **Example**
 * @include ex_smooth_fill_holes.cpp
 *
 * **Result**
 * @image html ex_smooth_fill_holes.png
 *
 * @see          {@link smooth_fill_smear_peaks}
 */
void smooth_fill_holes(Array &array, int ir);
void smooth_fill_holes(Array       &array,
                       int          ir,
                       const Array *p_mask); ///< @overload

/**
 * @brief Apply smoothing to smear peaks (elliptic convex surfaces).
 *
 * This function applies a smoothing filter designed to reduce the prominence of
 * peaks or convex surfaces in the input array. It smooths out convex regions,
 * creating a more level surface.
 *
 * @param array  Input array representing the heightmap or image to be
 *               processed.
 * @param ir     Filter radius, which controls the extent of the smoothing
 *               effect.
 * @param p_mask Optional filter mask, expected in the range [0, 1]. If
 *               provided, smoothing is applied according to this mask. If not
 *               provided, the entire array is processed.
 *
 * **Example**
 * @include ex_smooth_fill_holes.cpp
 *
 * **Result**
 * @image html ex_smooth_fill_holes.png
 *
 * @see          {@link smooth_fill_holes}
 */
void smooth_fill_smear_peaks(Array &array, int ir);
void smooth_fill_smear_peaks(Array       &array,
                             int          ir,
                             const Array *p_mask); ///< @overload

/**
 * @brief Applies a localized smoothstep operation to the provided array.
 *
 * This function modifies the input array using a localized smoothstep
 * operation. It calculates the local minimum and maximum values within a radius
 * (`ir`) and smooths the values in the array based on these bounds.
 *
 * @param array The 2D array to be smoothed.
 * @param ir    The radius used to compute the local minimum and maximum values.
 *
 * **Example**
 * @include ex_smoothstep_local.cpp
 *
 * **Result**
 * @image html ex_smoothstep_local.png
 */
void smoothstep_local(Array &array, int ir);

/**
 * @brief Applies a localized smoothstep operation to the provided array with an
 * optional mask.
 *
 * If a mask is provided, the function blends the smoothed values with the
 * original array using the mask. Otherwise, it directly applies the localized
 * smoothstep operation.
 *
 * @param array  The 2D array to be smoothed.
 * @param ir     The radius used to compute the local minimum and maximum
 *               values.
 * @param p_mask A pointer to an optional mask array. If provided, the smoothed
 *               array is blended with the original using this mask.
 *
 * **Example**
 * @include ex_smoothstep_local.cpp
 *
 * **Result**
 * @image html ex_smoothstep_local.png
 */
void smoothstep_local(Array &array, int ir, const Array *p_mask);

/**
 * @brief Steepen (or flatten) the array map.
 *
 * This function applies a steepening effect to the input array, enhancing the
 * contrast of the gradients or flattening features based on the given filter
 * amplitude. The effect is controlled by the `scale` parameter, and the
 * `ir` parameter defines the filtering radius for computing the array
 * gradients.
 *
 * @param array  Input array representing the heightmap or image to be
 *               processed.
 * @param scale  Filter amplitude that determines the extent of steepening or
 *               flattening.
 * @param p_mask Optional filter mask, expected in the range [0, 1]. If
 *               provided, steepening is applied according to this mask. If not
 *               provided, the entire array is processed.
 * @param ir     Filtering radius of the array gradients, which influences the
 *               extent of the effect.
 *
 * **Example**
 * @include ex_steepen.cpp
 *
 * **Result**
 * @image html ex_steepen.png
 */
void steepen(Array &array, float scale, int ir = 8);
void steepen(Array       &array,
             float        scale,
             const Array *p_mask,
             int          ir = 8); ///< @overload

/**
 * @brief Steepen array values by applying a nonlinear convection operator in a
 * given direction.
 *
 * This function applies a nonlinear convection operator to the input array,
 * simulating the effect described by Burger's equation. The steepening effect
 * is applied in the specified direction (angle) and is controlled by several
 * parameters including the number of iterations, smoothing radius, and time
 * step. The input array values are expected to be in the range [-1, 1].
 *
 * @param array      Input array with elements expected to be in the range [-1,
 *                   1].
 * @param angle      Steepening direction in degrees.
 * @param iterations Number of iterations to perform the convection process.
 * @param ir         Smoothing radius applied to the array values before
 *                   differentiation.
 * @param dt         Time step for the convection process, can be chosen smaller
 *                   than 1 for finer tuning of the steepening effect.
 * @param p_mask     Optional filter mask, expected in the range [0, 1]. If
 *                   provided, steepening is applied according to this mask. If
 *                   not provided, the entire array is processed.
 *
 * **Example**
 * @include ex_steepen_convective.cpp
 *
 * **Result**
 * @image html ex_steepen_convective.png
 */
void steepen_convective(Array &array,
                        float  angle,
                        int    iterations = 1,
                        int    ir = 0,
                        float  dt = 0.1f);

void steepen_convective(Array       &array,
                        float        angle,
                        const Array *p_mask,
                        int          iterations = 1,
                        int          ir = 0,
                        float        dt = 0.1f); ///< @overload

/**
 * @brief Applies a terrace effect to the values in an array.
 *
 * This function adjusts the values in the `array` by applying a terrace or
 * stepped effect, often used for terrain generation or other natural-looking
 * height variations. The terrace effect is controlled by several parameters
 * such as the number of levels, gain, noise ratio, and optional min/max range.
 * The noise is applied to levels for added randomness, and a gain function is
 * applied to smooth the transitions.
 *
 * @param array       The array of values to modify with the terrace effect.
 * @param seed        Seed value for random number generation.
 * @param nlevels     Number of terrace levels to apply.
 * @param gain        Gain factor for controlling the sharpness of the terrace
 *                    levels, in [0, 1].
 * @param noise_ratio Ratio of noise applied to each terrace level, except the
 *                    first and last.
 * @param p_noise     Optional noise array to introduce additional variation per
 *                    element.
 * @param vmin        Minimum value for terracing; if less than `vmax`, will be
 *                    auto-determined.
 * @param vmax        Maximum value for terracing; if less than `vmin`, will be
 *                    auto-determined.
 *
 * @note If `p_noise` is provided, each value in `array` is transformed using
 * both the original value and the corresponding noise value from `p_noise`.
 *
 * **Example**
 * @include ex_terrace.cpp
 *
 * **Result**
 * @image html ex_terrace.png
 */
void terrace(Array        &array,
             std::uint32_t seed,
             int           nlevels,
             float         gain = 0.9f,
             float         noise_ratio = 0.f,
             const Array  *p_noise = nullptr,
             float         vmin = 0.f,
             float         vmax = -1.f);

/**
 * @brief Applies a terrace effect to an array with optional masking.
 *
 * This overloaded version of the terrace function modifies `array` based on the
 * terrace levels, gain, noise ratio, and an optional mask array. If a mask is
 * provided, the terrace effect is applied conditionally based on the mask
 * values.
 *
 * @param array       The array of values to modify with the terrace effect.
 * @param seed        Seed value for random number generation.
 * @param nlevels     Number of terrace levels to apply.
 * @param p_mask      Optional mask array. If provided, blends the terrace
 *                    effect with original values based on the mask.
 * @param gain        Gain factor for controlling the sharpness of the terrace
 *                    levels, in [0, 1].
 * @param noise_ratio Ratio of noise applied to each terrace level, except the
 *                    first and last.
 * @param p_noise     Optional noise array to introduce additional variation per
 *                    element.
 * @param vmin        Minimum value for terracing; if less than `vmax`, will be
 *                    auto-determined.
 * @param vmax        Maximum value for terracing; if less than `vmin`, will be
 *                    auto-determined.
 *
 * This function:
 * - If no mask is provided, directly applies the terrace effect using the first
 * terrace overload.
 * - If a mask is provided, creates a temporary copy of `array`, applies the
 * terrace effect to it, and then interpolates between `array` and the modified
 * copy based on mask values.
 *
 * @note The mask array allows for blending the terrace effect with the original
 * array for more localized effects.
 *
 * **Example**
 * @include ex_terrace.cpp
 *
 * **Result**
 * @image html ex_terrace.png
 */
void terrace(Array        &array,
             std::uint32_t seed,
             int           nlevels,
             const Array  *p_mask,
             float         gain = 0.9f,
             float         noise_ratio = 0.f,
             const Array  *p_noise = nullptr,
             float         vmin = 0.f,
             float         vmax = -1.f);

} // namespace hmap

namespace hmap::gpu
{

/*! @brief See hmap::adaptive_relief */
void adaptive_relief(Array &array,
                     float  strength = 1.f,
                     float  clamp_ratio = 1.f);
void adaptive_relief(Array       &array,
                     const Array *p_mask,
                     float        strength = 1.f,
                     float        clamp_ratio = 1.f); ///< @overload

/**
 * @brief Apply a bilateral filter to an array.
 *
 * This function performs edge-preserving smoothing by combining a spatial
 * kernel (2D) and a value-domain kernel (1D). The spatial kernel controls the
 * influence of neighboring samples based on their position, while the value
 * kernel weights samples according to their value difference.
 *
 * @param  array                  Input array to be filtered.
 * @param  kernel2d               2D spatial kernel defining neighborhood
 *                                weights.
 * @param  kernel1d               1D kernel applied in the value domain.
 * @param  kernel1d_value_scaling Scaling factor applied to value differences
 *                                before sampling the 1D kernel.
 *
 * @return                        A new array containing the filtered result.
 *
 * **Example**
 * @include ex_bilateral_filter.png
 *
 * **Result**
 * @image html ex_bilateral_filter.png
 */
Array bilateral_filter(const Array &array,
                       const Array &kernel2d,
                       const Array &kernel1d,
                       float        kernel1d_value_scaling);

/**
 * @brief Apply a bilateral filter to an array using a radius-based Gaussian
 * kernel.
 *
 * This overload builds the spatial kernel implicitly from the given radius. The
 * value-domain kernel behavior is controlled by the same scaling factor as the
 * main overload.
 *
 * @param  array                  Input array to be filtered.
 * @param  ir                     Radius of the spatial neighborhood.
 * @param  kernel1d_value_scaling Scaling factor applied to value differences
 *                                before computing value-domain weights.
 *
 * @return                        A new array containing the filtered result.
 *
 * @overload
 *
 * **Example**
 * @include ex_bilateral_filter.png
 *
 * **Result**
 * @image html ex_bilateral_filter.png
 */
Array bilateral_filter(const Array &array,
                       int          ir,
                       float        kernel1d_value_scaling); ///< @overload

/**
 * @brief Applies a directional blur to the input array.
 *
 * Blurs values along a direction specified per element by the angle array,
 * using a given radius and number of sampling steps.
 *
 * @param array  Input array to be blurred (modified in place).
 * @param radius Blur radius.
 * @param angle  Array of angles (in radians) defining blur direction per
 *               element.
 * @param steps  Number of samples along the blur direction (default: 32).
 *
 * **Example**
 * @include ex_directional_blur.png
 *
 * **Result**
 * @image html ex_directional_blur.png
 */
void directional_blur(Array       &array,
                      float        radius,
                      const Array &angle,
                      int          steps = 32);

void directional_blur(Array       &array,
                      float        radius,
                      const Array &angle,
                      const Array *p_mask,
                      int          steps = 32); ///< @overload

/*! @brief See hmap::expand */
void expand(Array &array, int ir, int iterations = 1);
void expand(Array       &array,
            int          ir,
            const Array *p_mask,
            int          iterations = 1); ///< @overload
void expand(Array       &array,
            const Array &kernel,
            int          iterations = 1); ///< @overload
void expand(Array       &array,
            const Array &kernel,
            const Array *p_mask,
            int          iterations = 1); ///< @overload

/*! @brief See hmap::gamma_correction_local */
void gamma_correction_local(Array &array, float gamma, int ir, float k = 0.1f);
void gamma_correction_local(Array       &array,
                            float        gamma,
                            int          ir,
                            const Array *p_mask,
                            float        k = 0.1f); ///< @overload

/*! @brief See hmap::laplace */
void laplace(Array &array, float sigma = 0.2f, int iterations = 3);
void laplace(Array       &array,
             const Array *p_mask,
             float        sigma = 0.2f,
             int          iterations = 3); ///< @overload

/*! @brief See hmap::mean_shift */
Array mean_shift(const Array &array,
                 int          ir,
                 float        talus,
                 int          iterations = 1,
                 bool         talus_weighted = true);

Array mean_shift(const Array &array,
                 int          ir,
                 float        talus,
                 const Array *p_mask,
                 int          iterations = 1,
                 bool         talus_weighted = true);

/*! @brief See hmap::median_3x3 */
void median_3x3(Array &array);
void median_3x3(Array &array, const Array *p_mask); ///< @overload

/*! @brief See hmap::median_pseudo */
Array median_pseudo(const Array &array, int ir);

/*! @brief See hmap::normal_displacement */
void normal_displacement(Array &array,
                         float  amount = 0.1f,
                         int    ir = 0,
                         bool   reverse = false);
void normal_displacement(Array       &array,
                         const Array *p_mask,
                         float        amount = 0.1f,
                         int          ir = 0,
                         bool         reverse = false); ///< @overload

/*! @brief See hmap::plateau */
void plateau(Array &array, const Array *p_mask, int ir, float factor);
void plateau(Array &array, int ir, float factor); ///< @overload

/**
 * @brief Projects array values along a given direction using talus attenuation.
 *
 * Applies a directional propagation (D8 convention) where values decay with
 * distance according to a talus coefficient. The computation is performed on
 * the GPU using an OpenCL kernel.
 *
 * @param  array     Input 2D array.
 * @param  talus     Talus coefficient controlling distance-based attenuation.
 * @param  direction Propagation direction (D8 convention, range [0–7]).
 *
 * @return           A new Array containing the projected values.
 *
 * @note The array is temporarily shifted to non-negative values to ensure
 * correctness of OpenCL atomic operations on floats.
 *
 * **Example**
 * @include ex_project_talus_along_direction.png
 *
 * **Result**
 * @image html ex_project_talus_along_direction.png
 */
Array project_talus_along_direction(const Array &array,
                                    const Array &talus,
                                    int          direction = 0,
                                    float        vmin = -FLT_MAX);

Array project_talus_along_direction(const Array &array,
                                    const Array &talus,
                                    const Array *p_mask,
                                    int          direction = 0,
                                    float        vmin = -FLT_MAX);

/*! @brief See hmap::shrink */
void shrink(Array &array, int ir, int iterations = 1);
void shrink(Array       &array,
            int          ir,
            const Array *p_mask,
            int          iterations = 1); ///< @overload
void shrink(Array       &array,
            const Array &kernel,
            int          iterations = 1); ///< @overload
void shrink(Array       &array,
            const Array &kernel,
            const Array *p_mask,
            int          iterations = 1); ///< @overload

/*! @brief See hmap::smooth_cpulse */
void smooth_cpulse(Array &array, int ir);

/*! @brief See hmap::smooth_cpulse */
void smooth_cpulse(Array &array, int ir, const Array *p_mask);

/*! @brief See hmap::smooth_cpulse_edge_removing */
void smooth_cpulse_edge_removing(Array &array,
                                 float  talus,
                                 float  talus_width,
                                 int    ir);

/*! @brief See hmap::smooth_fill */
void smooth_fill(Array &array,
                 int    ir,
                 float  k = 0.1f,
                 Array *p_deposition_map = nullptr);
void smooth_fill(Array       &array,
                 int          ir,
                 const Array *p_mask,
                 float        k = 0.1f,
                 Array       *p_deposition_map = nullptr); ///< @overload

/*! @brief See hmap::smooth_fill_holes */
void smooth_fill_holes(Array &array, int ir);
void smooth_fill_holes(Array       &array,
                       int          ir,
                       const Array *p_mask); ///< @overload

/*! @brief See hmap::smooth_fill_smear_peaks */
void smooth_fill_smear_peaks(Array &array, int ir);
void smooth_fill_smear_peaks(Array       &array,
                             int          ir,
                             const Array *p_mask); ///< @overload

/**
 * @brief Spectral equalizer for 2D fields using a multiscale Gaussian pyramid.
 *
 * Decomposes the input array into log-spaced frequency bands using Gaussian
 * smoothing, then recombines them using normalized user-defined weights.
 *
 * @param  array   Input 2D field.
 * @param  weights Spectral gains (mapped from low to high frequencies).
 * @param  ir_min  Minimum Gaussian radius.
 * @param  ir_max  Maximum Gaussian radius.
 * @return         Filtered array with adjusted spectral content.
 *
 * **Example**
 * @include ex_spectral_equalizer.png
 *
 * **Result**
 * @image html ex_spectral_equalizer.png
 */
Array spectral_equalizer(const Array              &array,
                         const std::vector<float> &weights,
                         int                       ir_min,
                         int                       ir_max);

Array spectral_equalizer(const Array              &array,
                         const std::vector<float> &weights,
                         int                       ir_min,
                         int                       ir_max,
                         const Array              *p_mask); ///< @overload

} // namespace hmap::gpu
