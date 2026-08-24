/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file colorize.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @brief Header file for functions related to applying colorization and
 * hillshading.
 *
 * This file contains the declarations of functions that are used to apply
 * colorization and hillshading effects to images and arrays. The functions
 * allow for creating grayscale, histogram-based, and colored images from
 * elevation data and other input arrays.
 *
 * @copyright Copyright (c) 2023 Otto Link
 */
#pragma once
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/texture.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

/**
 * @enum MixMethod
 * @brief Methods for mixing/blending texture colors.
 *
 * @example ex_mixbox.cpp This example shows the visual differences between
 * linear, square-root, and physical pigment-based (Mixbox) color blending.
 *
 * Output results:
 * @image html ex_mixbox_linear.png "Linear Blending"
 * @image html ex_mixbox_sqrt.png "Square-root Average Blending"
 * @image html ex_mixbox_mixbox.png "Mixbox Blending"
 */
enum MixMethod : int
{
  MM_LINEAR,   ///< Linear RGB mixing.
  MM_SQRT_AVG, ///< Square-root average mixing.
  MM_MIXBOX    ///< Physical pigment-based mixing (Mixbox).
};

enum NormalMapBlendingMethod : int
{
  NMAP_LINEAR,
  NMAP_DERIVATIVE,
  NMAP_UDN,
  NMAP_UNITY,
  NMAP_WHITEOUT
};

static std::map<std::string, int> normal_map_blending_method_as_string = {
    {"Linear", NMAP_LINEAR},
    {"Partial derivative", NMAP_DERIVATIVE},
    {"Unreal Developer Network", NMAP_UDN},
    {"Unity", NMAP_UNITY},
    {"Whiteout", NMAP_WHITEOUT},
};

struct ColorAdjust
{
  float in_min = 0.0f;
  float in_max = 1.0f;
  float exposure = 0.0f;
  float contrast = 1.0f;
  float saturation = 1.0f;
  float temperature = 0.0f;
  float gamma = 1.f;
  float dither_amp = 0.f;
  bool  filmic_tonemap = false;
  bool  aces_tonemap = false;
  bool  agx_tonemap = false;
};

/**
 * @enum Cmap
 * @brief Enumeration for different colormap types.
 *
 * This enumeration is defined in the highmap/colormap.hpp file and is used for
 * selecting the colormap to be applied during colorization.
 */
enum Cmap : int; // highmap/colormap.hpp

/**
 * @brief Apply hillshading to a Tensor image.
 *
 * This function applies a hillshading effect to the provided tensor image using
 * the elevation data in the input array. The effect is controlled by the power
 * exponent and can be scaled by specifying the minimum and maximum values.
 *
 * @param img      Input image represented as a Tensor.
 * @param array    Elevation data used to generate the hillshading effect.
 * @param vmin     Minimum value for scaling the hillshading effect (default is
 *                 0.f).
 * @param vmax     Maximum value for scaling the hillshading effect (default is
 *                 1.f).
 * @param exponent Power exponent applied to the hillshade values (default is
 *                 1.f).
 */
void apply_hillshade(Texture     &img,
                     const Array &array,
                     float        vmin = 0.f,
                     float        vmax = 1.f,
                     float        exponent = 1.f);

/**
 * @brief Apply hillshading to an 8-bit image.
 *
 * This function applies a hillshading effect to a vector of 8-bit image data.
 * The effect is computed based on the elevation data in the input array. The
 * image can be either RGB or RGBA format depending on the is_img_rgba flag.
 *
 * @param img         Input image represented as a vector of 8-bit data.
 * @param array       Elevation data used to generate the hillshading effect.
 * @param vmin        Minimum value for scaling the hillshading effect (default
 *                    is 0.f).
 * @param vmax        Maximum value for scaling the hillshading effect (default
 *                    is 1.f).
 * @param exponent    Power exponent applied to the hillshade values (default is
 *                    1.f).
 * @param is_img_rgba Boolean flag indicating if the input image has an alpha
 *                    channel (default is false).
 */
void apply_hillshade(std::vector<uint8_t> &img,
                     const Array          &array,
                     float                 vmin = 0.f,
                     float                 vmax = 1.f,
                     float                 exponent = 1.f,
                     bool                  is_img_rgba = false);

/**
 * @brief Apply color correction and tonemapping to RGB arrays.
 *
 * Performs per-pixel adjustments including levels, exposure, tonemapping
 * (Reinhard, ACES, or AGX), contrast, saturation, temperature, gamma, and
 * optional dithering.
 *
 * @param r           Red channel array (modified in place).
 * @param g           Green channel array (modified in place).
 * @param b           Blue channel array (modified in place).
 * @param param       Color adjustment parameters.
 * @param dither_seed Seed used for dithering noise.
 *
 * **Example**
 * @include ex_color_adjust.cpp
 *
 * **Result**
 * @image html ex_color_adjust.png
 */
void color_adjust(Array      &r,
                  Array      &g,
                  Array      &b,
                  ColorAdjust param,
                  glm::ivec2  dither_seed);

/**
 * @brief Create a mask of pixels matching a target color within a tolerance.
 *
 * @param  r         Red channel array.
 * @param  g         Green channel array.
 * @param  b         Blue channel array.
 * @param  color     Target RGB color.
 * @param  tolerance Allowed color difference (0 = exact match).
 *
 * @return           Array Mask array (1 for match, 0 otherwise), same size as
 *                   inputs.
 *
 * **Example**
 * @include ex_color_match_mask.cpp
 *
 * **Result**
 * @image html ex_color_match_mask.png
 *
 * See unit tests: @ref test_color_match_mask.cpp
 */
Array color_match_mask(const Array     &r,
                       const Array     &g,
                       const Array     &b,
                       const glm::vec3 &color,
                       float            tolerance = 0.f);

/**
 * @brief Apply colorization to an array.
 *
 * This function applies a colormap to the input array and optionally applies
 * hillshading and noise. The colorization can be reversed, and the function
 * returns a Tensor representing the colorized image.
 *
 * @param  array       Input array to be colorized.
 * @param  vmin        Minimum value for scaling the colormap.
 * @param  vmax        Maximum value for scaling the colormap.
 * @param  cmap        Colormap to be applied (using the Cmap enum).
 * @param  hillshading Boolean flag to apply hillshading.
 * @param  reverse     Boolean flag to reverse the colormap (default is false).
 * @param  p_noise     Optional pointer to a noise array (default is nullptr).
 * @return             Tensor Colorized Tensor image.
 */
Texture colorize(const Array &array,
                 float        vmin,
                 float        vmax,
                 int          cmap,
                 bool         hillshading,
                 bool         reverse = false,
                 const Array *p_noise = nullptr);

/**
 * @brief Colorize a scalar field into a Texture using a custom colormap.
 *
 * @param  array           Input scalar values.
 * @param  vmin            Lower bound for normalization.
 * @param  vmax            Upper bound for normalization.
 * @param  positions       Normalized color positions.
 * @param  colormap_colors RGB(A) colors per position.
 * @param  reverse         Reverse colormap mapping.
 * @param  p_noise         Optional noise for dithering.
 * @return                 Colorized Texture.
 *
 * **Example**
 * @include ex_texture_mix_luminance.cpp
 *
 * @see                    tests/src/test_texture.cpp
 *                         (TextureColorize.CustomColormapAndOperations)
 */
Texture colorize(const Array                  &array,
                 float                         vmin,
                 float                         vmax,
                 const std::vector<float>     &positions,
                 const std::vector<glm::vec3> &colormap_colors,
                 bool                          reverse = false,
                 const Array                  *p_noise = nullptr);

/**
 * @brief Colorize two scalar fields into a Texture using two custom colormaps
 * and a color mixing method.
 *
 * Each scalar field is independently normalized using its corresponding range
 * and mapped to a color using its custom colormap. The two resulting colors are
 * then combined according to the specified mixing method.
 *
 * @param  a1               First input scalar field.
 * @param  a2               Second input scalar field.
 * @param  range1           Lower and upper bounds for normalization of the
 *                          first field.
 * @param  range2           Lower and upper bounds for normalization of the
 *                          second field.
 * @param  positions1       Normalized color positions for the first colormap.
 * @param  positions2       Normalized color positions for the second colormap.
 * @param  colormap_colors1 RGB colors for the first colormap.
 * @param  colormap_colors2 RGB colors for the second colormap.
 * @param  method           Method used to mix the two resulting colors.
 * @param  reverse1         Reverse the first colormap mapping.
 * @param  reverse2         Reverse the second colormap mapping.
 * @param  p_noise1         Optional noise added to the first scalar field
 *                          before normalization.
 * @param  p_noise2         Optional noise added to the second scalar field
 *                          before normalization.
 * @return                  Colorized Texture containing the mixed colors.
 *
 * **Example**
 * @include ex_texture_colorize_bivariate.cpp
 *
 * **Result**
 * @image html ex_colorize_bivariate0.png
 * @image html ex_colorize_bivariate1.png
 *
 * @see                     tests/src/test_texture.cpp
 */
Texture colorize_bivariate(const Array                  &a1,
                           const Array                  &a2,
                           const glm::vec2               range1,
                           const glm::vec2               range2,
                           const std::vector<float>     &positions1,
                           const std::vector<float>     &positions2,
                           const std::vector<glm::vec3> &colormap_colors1,
                           const std::vector<glm::vec3> &colormap_colors2,
                           MixMethod                     method = MM_SQRT_AVG,
                           bool                          reverse1 = false,
                           bool                          reverse2 = false,
                           const Array                  *p_noise1 = nullptr,
                           const Array                  *p_noise2 = nullptr);

/**
 * @brief Convert an array to a grayscale image.
 *
 * This function converts the input array to an 8-bit grayscale Texture image.
 *
 * @param  array Input array.
 * @return       Texture Grayscale Texture image.
 */
Texture colorize_grayscale(const Array &array);

/**
 * @brief Convert an array to a histogram-based grayscale image.
 *
 * This function converts the input array to an 8-bit grayscale Texture image
 * using a histogram-based method for enhanced contrast.
 *
 * @param  array Input array.
 * @return       Texture Grayscale Texture image with histogram-based contrast.
 */
Texture colorize_histogram(const Array &array);

/**
 * @brief Colorizes a slope height heatmap based on the gradient norms of a
 * given array.
 *
 * This function computes a colorized heatmap using a two-dimensional histogram
 * that considers the gradient norm of the input array. It normalizes both the
 * input array values and their corresponding gradient norms, calculates a 2D
 * histogram, and then applies a colormap to visualize the heatmap.
 *
 * @param  array The input array for which the slope height heatmap is computed.
 *               This should be a 2D array representing height values.
 * @param  cmap  An integer representing the colormap to be used for
 *               colorization. Colormap options depend on the colorization
 *               function used internally.
 *
 * @return       Texture representing the colorized heatmap, which visualizes
 *               the distribution of height values and their corresponding
 *               gradient norms.
 *
 * @details
 * - The function normalizes the height values and gradient norms independently
 * using the minimum and maximum values.
 * - A 2D histogram is constructed, where each bin corresponds to a pair of
 * normalized height and gradient values.
 * - The colormap is then applied to this histogram to produce a colorized
 * output.
 *
 * @note If the input array has a constant value (i.e., min == max), no
 * normalization is applied, and the function may not produce meaningful
 * results.
 *
 * @warning Ensure that the input array is non-empty and has valid dimensions.
 *
 * **Example**
 * @include ex_colorize_slope_height_heatmap.cpp
 *
 * **Result**
 * @image html ex_colorize_slope_height_heatmap.png
 */
Texture colorize_slope_height_heatmap(const Array &array, int cmap);

/**
 * @brief Combine two arrays into a colored image.
 *
 * This function takes two input arrays and combines them into a single 8-bit
 * colored Texture image. The resulting image uses the data from both arrays to
 * create a composite color representation.
 *
 * @param  array1 First input array.
 * @param  array2 Second input array.
 * @return        Texture Colorized Texture image.
 *
 * **Example**
 * @include ex_colorize_vec2.cpp
 *
 * **Result**
 * @image html ex_colorize_vec2.png
 */
Texture colorize_vec2(const Array &array1, const Array &array2);

/**
 * @brief Compute luminance from a Texture.
 *
 * Computes a grayscale luminance array from the RGB channels of a texture.
 *
 * @param  tex Input texture (expects at least 3 channels).
 * @return     Luminance Array.
 *
 * **Example**
 * @include ex_texture_mix_luminance.cpp
 *
 * @include ex_mixbox.cpp
 *
 * **Result**
 * @image html ex_mixbox_linear.png
 * @image html ex_mixbox_mixbox.png
 * @image html ex_mixbox_sqrt.png
 *
 * @see        tests/src/test_texture.cpp
 *             (TextureColorize.CustomColormapAndOperations)
 */
Array luminance(const Texture &tex);

/**
 * @brief Mix two textures into an output texture.
 * @param  tex1   First input texture.
 * @param  tex2   Second input texture.
 * @param  method Mixing method to use.
 * @return        Mixed Texture.
 *
 * **Example**
 * @include ex_texture_mix_luminance.cpp
 *
 * @include ex_mixbox.cpp
 *
 * **Result**
 * @image html ex_mixbox_linear.png
 * @image html ex_mixbox_mixbox.png
 * @image html ex_mixbox_sqrt.png
 *
 * @see           tests/src/test_texture.cpp
 *                (TextureColorize.CustomColormapAndOperations)
 */
Texture mix(const Texture &tex1,
            const Texture &tex2,
            MixMethod      method = MM_SQRT_AVG);

/**
 * @brief Mix a list of textures sequentially.
 * @param  texs   Input list of textures.
 * @param  method Mixing method to use.
 * @return        Mixed Texture.
 */
Texture mix(const std::vector<const Texture *> &texs,
            MixMethod                           method = MM_SQRT_AVG);

/**
 * @brief Blend two normal maps into a single output normal map.
 *
 * Combines a base normal map with a detail normal map using the specified
 * blending method and scaling factor.
 *
 * @param  nmap_base       Base normal map texture.
 * @param  nmap_detail     Detail normal map texture.
 * @param  detail_scaling  Strength of the detail normal map.
 * @param  blending_method Normal map blending method.
 * @return                 Blended Texture normal map.
 *
 * @see                    tests/src/test_texture.cpp
 *                         (TextureColorize.CustomColormapAndOperations)
 */
Texture mix_normal_map(const Texture          &nmap_base,
                       const Texture          &nmap_detail,
                       float                   detail_scaling,
                       NormalMapBlendingMethod blending_method);

} // namespace hmap

#include "highmap/virtual_array/virtual_texture.hpp"
