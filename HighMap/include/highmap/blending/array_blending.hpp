/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file blending.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Provides functions for various array blending modes.
 *
 * This header file contains functions for blending two arrays using different
 * techniques and modes. The available blending modes include exclusion, negate,
 * overlay, soft light, and gradient-based blending. These functions are useful
 * for image processing and heightmap manipulations where blending effects are
 * required.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */

#pragma once
#include "highmap/array.hpp"
#include "highmap/texture.hpp"

namespace hmap
{

/**
 * @brief Return the 'exclusion' blending of two arrays.
 *
 * The exclusion blend mode creates a blend that is similar to the difference
 * blend mode but with lower contrast. It produces an effect that is often used
 * in image editing for special effects. For more details, see
 * https://en.wikipedia.org/wiki/Blend_modes.
 *
 * @warning Array values are expected to be in the range [0, 1].
 *
 * @param  array1 First input array.
 * @param  array2 Second input array.
 * @return        Array Resulting array after applying the exclusion blend mode.
 *
 * **Example**
 * @include ex_blend.cpp
 *
 * **Result**
 * @image html ex_blend0.png
 * @image html ex_blend1.png
 * @image html ex_blend2.png
 * @image html ex_blend3.png
 * @image html ex_blend4.png
 *
 * @see           {@link blend_exclusion}, {@link blend_negate}, {@link
 *                blend_overlay},
 *      {@link blend_soft}
 */
Array blend_exclusion(const Array &array1, const Array &array2);

/**
 * @brief Return the blending of two arrays based on their gradients.
 *
 * This function blends two arrays by considering their gradients, which can be
 * useful for producing smooth transitions between the arrays based on their
 * directional changes.
 *
 * @param  array1 First input array.
 * @param  array2 Second input array.
 * @param  ir     Filtering radius in pixels.
 * @return        Array Resulting array after blending based on gradients.
 */
Array blend_gradients(const Array &array1, const Array &array2, int ir = 4);

/**
 * @brief Return the 'negate' blending of two arrays.
 *
 * The negate blend mode inverts the colors of one of the arrays relative to the
 * other. It creates a high-contrast effect that can be useful for emphasizing
 * differences between two arrays.
 *
 * @param  array1 First input array.
 * @param  array2 Second input array.
 * @return        Array Resulting array after applying the negate blend mode.
 *
 * @see           {@link blend_exclusion}, {@link blend_negate}, {@link
 *                blend_overlay},
 *      {@link blend_soft}
 */
Array blend_negate(const Array &array1, const Array &array2);

/**
 * @brief Return the 'overlay' blending of two arrays.
 *
 * The overlay blend mode combines the colors of two arrays in a way that
 * enhances contrast and emphasizes details.
 *
 * @warning Array values are expected to be in the range [0, 1].
 *
 * @param  array1 First input array.
 * @param  array2 Second input array.
 * @return        Array Resulting array after applying the overlay blend mode.
 *
 * @see           {@link blend_exclusion}, {@link blend_negate}, {@link
 *                blend_overlay},
 *      {@link blend_soft}
 */
Array blend_overlay(const Array &array1, const Array &array2);

/**
 * @brief Return the 'soft' blending of two arrays.
 *
 * The soft blend mode applies a soft light effect to the blending of two
 * arrays, based on the soft light mode described in
 * <a href=http://www.pegtop.net/delphi/articles/blendmodes/softlight.htm>Pegtop
 * soft light mode</a>. It creates a subtle blend with reduced contrast.
 *
 * @warning Array values are expected to be in the range [0, 1].
 *
 * @param  array1 First input array.
 * @param  array2 Second input array.
 * @return        Array Resulting array after applying the soft light blend
 *                mode.
 *
 * @see           {@link blend_exclusion}, {@link blend_negate}, {@link
 *                blend_overlay},
 *      {@link blend_soft}
 */
Array blend_soft(const Array &array1, const Array &array2);

/**
 * @brief Return the power law blending of a set of arrays.
 *
 * Blends two (or more) arrays using power law averaging:
 * \f[
 * f_\text{blend} = \frac{\sum_i f_i^{\alpha+1}}{\sum_i f_i^\alpha}
 * \f] The behavior varies with the bias \f$\alpha\f$ from a simple average of
 * the heights at each node when \f$\alpha=0\f$, approaching \f$\max(f_i)\f$ as
 * \f$\alpha \to \infty\f$.
 *
 * @param  arrays Array references to blend.
 * @param  alpha  Power law exponent (default: 0.0).
 * @return        Array Resulting array after power law blending.
 *
 * **Example**
 * @include ex_blend_power_law.cpp
 *
 * **Result**
 * @image html ex_blend_power_law0.png
 * @image html ex_blend_power_law1.png
 * @image html ex_blend_power_law2.png
 * @image html ex_blend_power_law3.png
 */
Array blend_power_law(const std::vector<const Array *> &arrays,
                      float                             alpha = 0.f);

/**
 * @brief Return the power law blending of two arrays.
 *
 * @param  array1 First input array.
 * @param  array2 Second input array.
 * @param  alpha  Power law exponent (default: 0.0).
 * @return        Array Resulting array after power law blending.
 */
Array blend_power_law(const Array &array1,
                      const Array &array2,
                      float        alpha = 0.f);

/**
 * @brief Return the mixing of a set of arrays based on a parameter `t`.
 *
 * This function mixes a set of arrays based on a mixing coefficient `t`
 * that determines the contribution of each array. The parameter `t` should be
 * in the range [0, 1], where `t` represents the blending ratio.
 *
 * @warning Values of array `t` are expected to be in the range [0, 1].
 *
 * @param  t      Mixing coefficient, in the range [0, 1].
 * @param  arrays References to the input arrays to be mixed.
 * @return        Array Resulting array after mixing.
 *
 * **Example**
 * @include ex_mixer.cpp
 *
 * **Result**
 * @image html ex_mixer.png
 */
Array mixer(const Array                      &t,
            const std::vector<const Array *> &arrays,
            float                             gain_factor = 1.f);

/**
 * @brief Transfers spatial details from a source array onto a target array.
 *
 * This function applies a high-pass spatial filter to the @p source array
 * (using a smoothing pulse of radius @p ir) to extract its high-frequency
 * content. That content is then scaled by @p amplitude and blended with the
 * @p target array. Optionally, the target can be prefiltered before blending.
 *
 * The result is a new array that combines the base structure of @p target with
 * the high-frequency details from @p source.
 *
 * @param  source              The input array providing high-frequency details.
 * @param  target              The input array providing the base structure.
 * @param  ir                  Radius (in pixels or samples) for the smoothing
 *                             filter used to separate low/high frequencies.
 * @param  amplitude           Scaling factor applied to the high-frequency
 *                             content extracted from @p source.
 * @param  target_prefiltering If true, @p target is smoothed before blending to
 *                             better integrate with the transferred details.
 *
 * @return                     A new Array containing the target data enhanced
 *                             by the source's high-frequency information.
 *
 * **Example**
 * @include ex_transfer.cpp
 *
 * **Result**
 * @image html ex_transfer.png
 */
Array transfer(const Array &source,
               const Array &target,
               int          ir,
               float        amplitude,
               bool         target_prefiltering = false);

} // namespace hmap

namespace hmap::gpu
{

/*! @brief See hmap::blend_gradients */
Array blend_gradients(const Array &array1, const Array &array2, int ir = 4);

/**
 * @brief Blends two arrays using Poisson blending with a brute-force solver.
 *
 * This function performs Poisson blending between @p array1 and @p array2 over
 * a specified number of iterations. Optionally, a mask can be provided to
 * control the blending regions.
 *
 * @param  array1     The first input array.
 * @param  array2     The second input array.
 * @param  iterations The number of iterations for the blending process
 *                    (default:
 *                    500).
 * @param  p_mask     Optional pointer to an array defining the blending mask.
 *                    If null, blending is applied globally.
 * @return            The blended array resulting from the Poisson blending
 *                    operation.
 *
 * **Example**
 * @include ex_blend_poisson_bf.cpp
 *
 * **Result**
 * @image html ex_blend_poisson_bf.png
 */
Array blend_poisson_bf(const Array &array1,
                       const Array &array2,
                       const int    iterations = 500,
                       const Array *p_mask = nullptr);

/**
 * @brief Blends two textures using Poisson blending with a brute-force solver.
 *
 * This function performs Poisson blending between @p texture1 and @p texture2
 * channel by channel over a specified number of iterations. Optionally, a mask
 * can be provided to control the blending regions.
 *
 * @param  texture1   The first input texture.
 * @param  texture2   The second input texture.
 * @param  iterations The number of iterations for the blending process
 *                    (default: 500).
 * @param  p_mask     Optional pointer to an array defining the blending mask.
 *                    If null, blending is applied globally.
 * @return            The blended texture resulting from the Poisson blending
 *                    operation.
 */
Texture blend_poisson_bf(const Texture &texture1,
                         const Texture &texture2,
                         const int      iterations = 500,
                         const Array   *p_mask = nullptr);

/*! @brief See hmap::transfer */
Array transfer(const Array &source,
               const Array &target,
               int          ir,
               float        amplitude,
               bool         target_prefiltering = false);

} // namespace hmap::gpu
