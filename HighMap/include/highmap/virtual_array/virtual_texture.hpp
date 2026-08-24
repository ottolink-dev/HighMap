/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */

/**
 * @file virtual_array.hpp
 * @author Otto Link (otto.link.bv@gmail.com)
 * @copyright Copyright (c) 2025
 */
#pragma once
#include <memory>
#include <string>
#include <vector>

#include <glm/vec2.hpp>
#include <glm/vec4.hpp>

#include "macrologger.h"

#include "highmap/colorize.hpp"
#include "highmap/texture.hpp"
#include "highmap/virtual_array/virtual_array.hpp"
#include "highmap/virtual_array/virtual_texture_storage.hpp"

namespace hmap
{

class VirtualTexture
{
public:
  // --- Ctor

  VirtualTexture() = default;

  VirtualTexture(glm::ivec2                   shape,
                 glm::vec4                    bbox,
                 glm::ivec2                   tile_shape,
                 int                          halo,
                 int                          channels,
                 std::unique_ptr<TileStorage> storage_proto);

  VirtualTexture(glm::ivec2  shape,
                 glm::vec4   bbox,
                 glm::ivec2  tile_shape,
                 int         halo,
                 int         channels,
                 StorageMode storage_mode);

  VirtualTexture(glm::ivec2  shape,
                 glm::ivec2  tile_shape,
                 int         halo,
                 int         channels,
                 StorageMode storage_mode);

  // --- Copy

  void copy_from(VirtualTexture &src, const ComputeMode &cm);

  // --- Channels data

  int                               channels() const;
  VirtualArray                     &channel(int c);
  const VirtualArray               &channel(int c) const;
  std::vector<VirtualArray *>       channels_ptr();
  std::vector<const VirtualArray *> channels_ptr() const;
  std::vector<VirtualArray>        &get_arrays();

  void fill(float value, const ComputeMode &cm);
  void fill(int c, float value, const ComputeMode &cm);

  // --- Converters

  void from_arrays(const std::vector<const Array *> &p_arrays,
                   const ComputeMode                &cm);

  std::vector<uint8_t> to_img_8bit(const glm::ivec2  &img_shape,
                                   const ComputeMode &c,
                                   bool               flip_y = false) const;

  void to_png(const glm::ivec2  &array_shape,
              const std::string &fname,
              const ComputeMode &cm,
              int                depth = CV_8U) const;

  void to_png(const std::string &fname,
              const ComputeMode &cm,
              int                depth = CV_8U) const;

  std::vector<float> to_raw(const ComputeMode &cm, bool flip_y = false);

  Texture to_texture(const glm::ivec2 &img_shape, const ComputeMode &cm) const;

  // --- Members

  glm::ivec2 shape;
  glm::vec4  bbox;
  glm::ivec2 tile_shape;
  int        halo;

private:
  std::vector<VirtualArray> arrays;
};

// functions

VirtualTexture convert_texture_channels(const VirtualTexture &src,
                                        int                   dst_channels,
                                        float                 fill_value,
                                        const ComputeMode    &cm);

// Forward declarations and imports for colorize / image operations on virtual
// textures
enum NormalMapBlendingMethod : int;
enum Cmap : int;

/**
 * @brief Colorize a scalar field into a VirtualTexture using a predefined
 * colormap.
 * @param out     Output virtual texture.
 * @param level   Input scalar values.
 * @param cm      Compute mode (CPU/GPU).
 * @param vmin    Lower bound for normalization.
 * @param vmax    Upper bound for normalization.
 * @param cmap    Colormap identifier.
 * @param p_alpha Optional alpha channel.
 * @param reverse Reverse colormap mapping.
 * @param p_noise Optional noise for dithering.
 *
 * **Example**
 * @include ex_virtual_texture.cpp
 */
void colorize(VirtualTexture    &out,
              VirtualArray      &level,
              const ComputeMode &cm,
              float              vmin,
              float              vmax,
              int                cmap,
              VirtualArray      *p_alpha = nullptr,
              bool               reverse = false,
              VirtualArray      *p_noise = nullptr);

/**
 * @brief Colorize a scalar field into a VirtualTexture using a custom colormap.
 * @param out             Output virtual texture.
 * @param level           Input scalar values.
 * @param cm              Compute mode (CPU/GPU).
 * @param vmin            Lower bound for normalization.
 * @param vmax            Upper bound for normalization.
 * @param positions       Normalized color positions.
 * @param colormap_colors RGB(A) colors per position.
 * @param p_alpha         Optional alpha channel.
 * @param reverse         Reverse colormap mapping.
 * @param p_noise         Optional noise for dithering.
 *
 * **Example**
 * @include ex_virtual_texture.cpp
 */
void colorize(VirtualTexture               &out,
              VirtualArray                 &level,
              const ComputeMode            &cm,
              float                         vmin,
              float                         vmax,
              const std::vector<float>     &positions,
              const std::vector<glm::vec3> &colormap_colors,
              VirtualArray                 *p_alpha = nullptr,
              bool                          reverse = false,
              VirtualArray                 *p_noise = nullptr);

/**
 * @brief Colorize two scalar fields into a VirtualTexture using two custom
 * colormaps and a color mixing method.
 * @param out              Output virtual texture.
 * @param a1               First input scalar field.
 * @param a2               Second input scalar field.
 * @param cm               Compute mode (CPU/GPU).
 * @param range1           Lower and upper bounds for normalization of the first
 *                         field.
 * @param range2           Lower and upper bounds for normalization of the
 *                         second field.
 * @param positions1       Normalized color positions for the first colormap.
 * @param positions2       Normalized color positions for the second colormap.
 * @param colormap_colors1 RGB colors for the first colormap.
 * @param colormap_colors2 RGB colors for the second colormap.
 * @param method           Method used to mix the two resulting colors.
 * @param reverse1         Reverse the first colormap mapping.
 * @param reverse2         Reverse the second colormap mapping.
 * @param p_noise1         Optional noise added to the first scalar field before
 *                         normalization.
 * @param p_noise2         Optional noise added to the second scalar field
 *                         before normalization.
 */
void colorize_bivariate(VirtualTexture               &out,
                        VirtualArray                 &a1,
                        VirtualArray                 &a2,
                        const ComputeMode            &cm,
                        glm::vec2                     range1,
                        glm::vec2                     range2,
                        const std::vector<float>     &positions1,
                        const std::vector<float>     &positions2,
                        const std::vector<glm::vec3> &colormap_colors1,
                        const std::vector<glm::vec3> &colormap_colors2,
                        MixMethod                     method = MM_SQRT_AVG,
                        bool                          reverse1 = false,
                        bool                          reverse2 = false,
                        VirtualArray                 *p_noise1 = nullptr,
                        VirtualArray                 *p_noise2 = nullptr);

/**
 * @brief Compute luminance from a VirtualTexture.
 *
 * Computes a grayscale luminance array from the RGB channels of a virtual
 * texture and stores the result in @p out.
 *
 * @param out Output luminance array.
 * @param tex Input virtual texture (expects at least 3 channels).
 * @param cm  Compute mode (execution and storage behavior).
 */
void luminance(VirtualArray &out, VirtualTexture &tex, const ComputeMode &cm);

/**
 * @brief Mix two virtual textures into an output virtual texture.
 * @param out    Output virtual texture.
 * @param tex1   First input virtual texture.
 * @param tex2   Second input virtual texture.
 * @param cm     Compute mode (CPU/GPU).
 * @param method Mixing method to use.
 *
 * **Example**
 * @include ex_virtual_texture.cpp
 */
void mix(VirtualTexture    &out,
         VirtualTexture    &tex1,
         VirtualTexture    &tex2,
         const ComputeMode &cm,
         MixMethod          method = MM_SQRT_AVG);

/**
 * @brief Mix multiple virtual textures sequentially.
 * @param out    Output virtual texture.
 * @param texs   Vector of input virtual textures.
 * @param cm     Compute mode (CPU/GPU).
 * @param method Mixing method to use.
 */
void mix(VirtualTexture                &out,
         std::vector<VirtualTexture *> &texs,
         const ComputeMode             &cm,
         MixMethod                      method = MM_SQRT_AVG);

/**
 * @brief Blend two normal maps into a single output virtual normal map.
 *
 * @param out             Output normal map texture.
 * @param nmap_base       Base normal map texture.
 * @param nmap_detail     Detail normal map texture.
 * @param cm              Compute mode (execution and storage behavior).
 * @param detail_scaling  Strength of the detail normal map.
 * @param blending_method Normal map blending method.
 */
void mix_normal_map(VirtualTexture         &out,
                    VirtualTexture         &nmap_base,
                    VirtualTexture         &nmap_detail,
                    const ComputeMode      &cm,
                    float                   detail_scaling,
                    NormalMapBlendingMethod blending_method);

template <typename Func>
void for_each_tile(VirtualTexture &tex, Func &&func, const ComputeMode &cm);

template <typename Func>
void for_each_pixel(VirtualTexture &tex, Func &&func, const ComputeMode &cm);

#include "highmap/virtual_array/virtual_texture.inl"

} // namespace hmap
