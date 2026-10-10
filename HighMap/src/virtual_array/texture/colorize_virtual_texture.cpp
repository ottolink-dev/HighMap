/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cstddef>
#include <functional>
#include <future>
#include <type_traits>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/colorize.hpp"
#include "highmap/colormaps.hpp"
#include "highmap/interpolate/interpolate1d.hpp"
#include "highmap/logger.hpp"
#include "highmap/operator.hpp"
#include "highmap/texture.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/virtual_array.hpp"
#include "highmap/virtual_array/virtual_texture.hpp"

#include "mixbox.h"

#include <format>

namespace hmap
{

void colorize(VirtualTexture    &out,
              VirtualArray      &level,
              const ComputeMode &cm,
              float              vmin,
              float              vmax,
              int                cmap,
              VirtualArray      *p_alpha,
              bool               reverse,
              VirtualArray      *p_noise)
{
  std::vector<glm::vec3> colors = get_colormap_data(cmap);

  std::vector<float> positions = linspace(0.f, 1.f, colors.size(), true);

  colorize(out,
           level,
           cm,
           vmin,
           vmax,
           positions,
           colors,
           p_alpha,
           reverse,
           p_noise);
}

void colorize(VirtualTexture               &out,
              VirtualArray                 &level,
              const ComputeMode            &cm,
              float                         vmin,
              float                         vmax,
              const std::vector<float>     &positions,
              const std::vector<glm::vec3> &colormap_colors,
              VirtualArray                 *p_alpha,
              bool                          reverse,
              VirtualArray                 *p_noise)
{
  if (out.channels() < 3)
  {
    hmap::log::error(
        "VirtualTexture must have at least 3 channels to be colorized");
    return;
  }

  // copy to reverse them if requested
  auto cpos = positions;
  auto colors = colormap_colors;

  if (reverse)
  {
    std::reverse(colors.begin(), colors.end());
    std::reverse(cpos.begin(), cpos.end());

    for (auto &p : cpos)
      p = 1.f - p;
  }

  // colorize fct
  auto lambda = [&](std::vector<const Array *> p_arrays_in,
                    std::vector<Array *>       p_arrays_out,
                    const TileRegion          &region)
  {
    const Array &z = *p_arrays_in[0];
    const Array *pa_noise = p_arrays_in[1];
    const Array *pa_alpha = p_arrays_in[2];
    Array       &r = *p_arrays_out[0];
    Array       &g = *p_arrays_out[1];
    Array       &b = *p_arrays_out[2];

    // color interpolators
    std::vector<float> cc_r, cc_g, cc_b;
    for (const auto &col : colors)
    {
      cc_r.push_back(col[0]);
      cc_g.push_back(col[1]);
      cc_b.push_back(col[2]);
    }

    Interpolator1D citp_r = hmap::Interpolator1D(
        cpos,
        cc_r,
        hmap::InterpolationMethod1D::LINEAR);

    Interpolator1D citp_g = hmap::Interpolator1D(
        cpos,
        cc_g,
        hmap::InterpolationMethod1D::LINEAR);

    Interpolator1D citp_b = hmap::Interpolator1D(
        cpos,
        cc_b,
        hmap::InterpolationMethod1D::LINEAR);

    // colorize
    for (int j = 0; j < region.shape.y; ++j)
      for (int i = 0; i < region.shape.x; ++i)
      {
        float v = z(i, j) + (pa_noise ? (*pa_noise)(i, j) : 0.f);
        v = (v - vmin) / (vmax - vmin);
        v = std::clamp(v, 0.f, 1.f);

        r(i, j) = citp_r(v);
        g(i, j) = citp_g(v);
        b(i, j) = citp_b(v);
      }

    // alpha channel if out has at least 4 channels
    if (p_arrays_out.size() >= 4 && p_arrays_out[3])
    {
      if (pa_alpha)
        *p_arrays_out[3] = *pa_alpha;
      else
        *p_arrays_out[3] = 1.f;
    }
  };

  // apply
  std::vector<const VirtualArray *> vas_in = {&level, p_noise, p_alpha};
  std::vector<VirtualArray *>       vas_out = out.channels_ptr();

  for_each_tile(vas_in, vas_out, lambda, cm);
}

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
                        MixMethod                     method,
                        bool                          reverse1,
                        bool                          reverse2,
                        VirtualArray                 *p_noise1,
                        VirtualArray                 *p_noise2)
{
  if (out.channels() < 3)
  {
    hmap::log::error(
        "VirtualTexture must have at least 3 channels to be colorized");
    return;
  }

  // colormap preparation
  auto prepare_colormap = [](const std::vector<float>     &positions,
                             const std::vector<glm::vec3> &colors,
                             bool                          reverse)
  {
    auto cpos = positions;
    auto colors_ = colors;

    if (reverse)
    {
      std::reverse(colors_.begin(), colors_.end());
      std::reverse(cpos.begin(), cpos.end());

      for (auto &p : cpos)
        p = 1.f - p;
    }

    return std::make_pair(std::move(cpos), std::move(colors_));
  };

  auto [cpos1,
        colors1] = prepare_colormap(positions1, colormap_colors1, reverse1);
  auto [cpos2,
        colors2] = prepare_colormap(positions2, colormap_colors2, reverse2);

  // setup color interpolators
  std::vector<float> cc_r1, cc_g1, cc_b1;
  for (const auto &col : colors1)
  {
    cc_r1.push_back(col[0]);
    cc_g1.push_back(col[1]);
    cc_b1.push_back(col[2]);
  }
  Interpolator1D citp_r1(cpos1, cc_r1, InterpolationMethod1D::LINEAR);
  Interpolator1D citp_g1(cpos1, cc_g1, InterpolationMethod1D::LINEAR);
  Interpolator1D citp_b1(cpos1, cc_b1, InterpolationMethod1D::LINEAR);

  std::vector<float> cc_r2, cc_g2, cc_b2;
  for (const auto &col : colors2)
  {
    cc_r2.push_back(col[0]);
    cc_g2.push_back(col[1]);
    cc_b2.push_back(col[2]);
  }
  Interpolator1D citp_r2(cpos2, cc_r2, InterpolationMethod1D::LINEAR);
  Interpolator1D citp_g2(cpos2, cc_g2, InterpolationMethod1D::LINEAR);
  Interpolator1D citp_b2(cpos2, cc_b2, InterpolationMethod1D::LINEAR);

  auto mix_colors = [method](const glm::vec3 &color1, const glm::vec3 &color2)
  {
    switch (method)
    {
    case MixMethod::MM_LINEAR: return 0.5f * (color1 + color2);

    case MixMethod::MM_SQRT_AVG:
      return glm::sqrt(0.5f * (color1 * color1 + color2 * color2));

    case MixMethod::MM_MIXBOX:
    {
      glm::vec3 cmix;
      mixbox_lerp_float(color1.x,
                        color1.y,
                        color1.z,
                        color2.x,
                        color2.y,
                        color2.z,
                        0.5f,
                        &cmix.x,
                        &cmix.y,
                        &cmix.z);
      return cmix;
    }

    default: return color1;
    }
  };

  float denom1 = (range1.y != range1.x) ? (range1.y - range1.x) : 1.f;
  float denom2 = (range2.y != range2.x) ? (range2.y - range2.x) : 1.f;

  auto lambda = [&](std::vector<const Array *> p_arrays_in,
                    std::vector<Array *>       p_arrays_out,
                    const TileRegion          &region)
  {
    const Array &za1 = *p_arrays_in[0];
    const Array &za2 = *p_arrays_in[1];
    const Array *pa_noise1 = p_arrays_in[2];
    const Array *pa_noise2 = p_arrays_in[3];
    Array       &r = *p_arrays_out[0];
    Array       &g = *p_arrays_out[1];
    Array       &b = *p_arrays_out[2];

    for (int j = 0; j < region.shape.y; ++j)
      for (int i = 0; i < region.shape.x; ++i)
      {
        float v1 = za1(i, j) + (pa_noise1 ? (*pa_noise1)(i, j) : 0.f);
        float v2 = za2(i, j) + (pa_noise2 ? (*pa_noise2)(i, j) : 0.f);

        v1 = (v1 - range1.x) / denom1;
        v2 = (v2 - range2.x) / denom2;

        v1 = std::clamp(v1, 0.f, 1.f);
        v2 = std::clamp(v2, 0.f, 1.f);

        glm::vec3 c1{citp_r1(v1), citp_g1(v1), citp_b1(v1)};
        glm::vec3 c2{citp_r2(v2), citp_g2(v2), citp_b2(v2)};

        glm::vec3 color = mix_colors(c1, c2);

        r(i, j) = color[0];
        g(i, j) = color[1];
        b(i, j) = color[2];
      }

    // alpha channel if 4 channels
    if (p_arrays_out.size() >= 4 && p_arrays_out[3])
    {
      *p_arrays_out[3] = 1.f;
    }
  };

  std::vector<const VirtualArray *> vas_in = {&a1, &a2, p_noise1, p_noise2};
  std::vector<VirtualArray *>       vas_out = out.channels_ptr();

  for_each_tile(vas_in, vas_out, lambda, cm);
}

void luminance(VirtualArray &out, VirtualTexture &tex, const ComputeMode &cm)
{
  if (tex.channels() < 3)
  {
    hmap::log::error("inputs mismatch, virtual textures must have 3 channels "
                     "for luminance.");
    return;
  }

  auto lambda = [](std::vector<const Array *> p_arrays_in,
                   std::vector<Array *>       p_arrays_out,
                   const TileRegion &)
  {
    hmap::Array       &lum = *p_arrays_out[0];
    const hmap::Array &r = *p_arrays_in[0];
    const hmap::Array &g = *p_arrays_in[1];
    const hmap::Array &b = *p_arrays_in[2];

    lum = 0.299f * r + 0.587f * g + 0.114f * b;
  };

  std::vector<const VirtualArray *> vas_in;
  for (const auto &ch : tex.get_arrays())
    vas_in.push_back(&ch);
  std::vector<VirtualArray *> vas_out = {&out};

  for_each_tile(vas_in, vas_out, lambda, cm);
}

void mix(VirtualTexture    &out,
         VirtualTexture    &tex1,
         VirtualTexture    &tex2,
         const ComputeMode &cm,
         MixMethod          method)
{
  // --- failsafe

  if (out.channels() != 4 || tex1.channels() != 4 || tex2.channels() != 4 ||
      out.shape != tex1.shape || out.shape != tex2.shape)
  {
    hmap::log::error("inputs mismatch, virtual textures must have 4 channels "
                     "and same shape.");
    return;
  }

  // --- colorize fct

  auto lambda = [method](std::vector<const Array *> p_arrays_in,
                         std::vector<Array *>       p_arrays_out,
                         const TileRegion &)
  {
    // tex1: in 0..3
    // tex2: in 4..7
    // out:  out 0..3

    // Construct temporary textures for the tiles
    Texture t1(*p_arrays_in[0],
               *p_arrays_in[1],
               *p_arrays_in[2],
               *p_arrays_in[3]);
    Texture t2(*p_arrays_in[4],
               *p_arrays_in[5],
               *p_arrays_in[6],
               *p_arrays_in[7]);

    // Mix using the standard mix function
    Texture blended = mix(t1, t2, method);

    // Write back the result to the output tile
    for (int c = 0; c < 4; ++c)
    {
      *p_arrays_out[c] = blended[c];
    }
  };

  // apply
  std::vector<const VirtualArray *> vas_in;
  for (const auto &ch : tex1.get_arrays())
    vas_in.push_back(&ch);
  for (const auto &ch : tex2.get_arrays())
    vas_in.push_back(&ch);

  std::vector<VirtualArray *> vas_out = out.channels_ptr();

  for_each_tile(vas_in, vas_out, lambda, cm);
}

void mix(VirtualTexture                &out,
         std::vector<VirtualTexture *> &texs,
         const ComputeMode             &cm,
         MixMethod                      method)
{
  if (texs.size() == 0) return;

  out.copy_from(*texs.front(), cm);

  for (size_t k = 1; k < texs.size(); k++)
    mix(out, out, *(texs[k]), cm, method);
}

void mix_normal_map(VirtualTexture         &out,
                    VirtualTexture         &nmap_base,
                    VirtualTexture         &nmap_detail,
                    const ComputeMode      &cm,
                    float                   detail_scaling,
                    NormalMapBlendingMethod blending_method)
{
  // output, also used to store first normal map
  out.copy_from(nmap_base, cm);

  // mix and then re-normalize values assuming a RGB channels
  // represent a normal vector
  auto lambda =
      [detail_scaling, blending_method](std::vector<const Array *> p_arrays_in,
                                        std::vector<Array *>       p_arrays_out,
                                        const TileRegion          &region)
  {
    Array *pa_r1 = p_arrays_out[0];
    Array *pa_g1 = p_arrays_out[1];
    Array *pa_b1 = p_arrays_out[2];

    const Array *pa_r2 = p_arrays_in[0];
    const Array *pa_g2 = p_arrays_in[1];
    const Array *pa_b2 = p_arrays_in[2];

    std::function<glm::vec3(glm::vec3 &, glm::vec3 &)> blending_fct;

    switch (blending_method)
    {
    case NormalMapBlendingMethod::NMAP_LINEAR:
    {
      blending_fct = [](glm::vec3 &n1, glm::vec3 &n2) { return n1 + n2; };
    }
    break;
    //
    case NormalMapBlendingMethod::NMAP_DERIVATIVE:
    {
      blending_fct = [](glm::vec3 &n1, glm::vec3 &n2)
      {
        glm::vec3 vn = glm::vec3(n1.x * n2.z + n2.x * n1.z,
                                 n1.y * n2.z + n2.y * n1.z,
                                 n1.z * n2.z);
        return vn;
      };
    }
    break;
    //
    case NormalMapBlendingMethod::NMAP_UDN:
    {
      blending_fct = [](glm::vec3 &n1, glm::vec3 &n2)
      {
        glm::vec3 vn = glm::vec3(n1.x + n2.x, n1.y + n2.y, n1.z);
        return vn;
      };
    }
    break;
      //
    case NormalMapBlendingMethod::NMAP_UNITY:
    {
      blending_fct = [](glm::vec3 &n1, glm::vec3 &n2)
      {
        glm::vec3 m0 = glm::vec3(n1.z, n1.x, -n1.x);
        glm::vec3 m1 = glm::vec3(n1.x, n1.z, -n1.y);
        glm::vec3 m2 = glm::vec3(n1.x, n1.y, n1.z);

        glm::vec3 vn = glm::vec3(n2.x * m0.x + n2.y * m1.x + n2.z * m2.x,
                                 n2.x * m0.y + n2.y * m1.y + n2.z * m2.y,
                                 n2.x * m0.z + n2.y * m1.z + n2.z * m2.z);
        return vn;
      };
    }
    break;
    //
    case NormalMapBlendingMethod::NMAP_WHITEOUT:
    default:
    {
      blending_fct = [](glm::vec3 &n1, glm::vec3 &n2)
      {
        glm::vec3 vn = glm::vec3(n1.x + n2.x, n1.y + n2.y, n1.z * n2.z);
        return vn;
      };
    }
    }

    for (int j = 0; j < region.shape.y; j++)
      for (int i = 0; i < region.shape.x; i++)
      {
        // do some rescaling because RGBA texture expected in [0, 1]
        // but normal vector expected in [-1, 1]

        glm::vec3 v111 = glm::vec3(1.f, 1.f, 1.f);
        glm::vec3 n1 = 2.f * glm::vec3((*pa_r1)(i, j),
                                       (*pa_g1)(i, j),
                                       (*pa_b1)(i, j)) -
                       v111;
        glm::vec3 n2 = 2.f * glm::vec3((*pa_r2)(i, j),
                                       (*pa_g2)(i, j),
                                       (*pa_b2)(i, j)) -
                       v111;

        n2.x *= detail_scaling;
        n2.y *= detail_scaling;
        n2.z *= detail_scaling;

        glm::vec3 vn = blending_fct(n1, n2);
        vn = glm::normalize(vn);

        (*pa_r1)(i, j) = 0.5f * vn.x + 0.5f;
        (*pa_g1)(i, j) = 0.5f * vn.y + 0.5f;
        (*pa_b1)(i, j) = 0.5f * vn.z + 0.5f;
      }
  };

  // apply
  std::vector<const VirtualArray *> vas_in;
  for (const auto &ch : nmap_detail.get_arrays())
    vas_in.push_back(&ch);

  std::vector<VirtualArray *> vas_out = out.channels_ptr();

  for_each_tile(vas_in, vas_out, lambda, cm);
}

} // namespace hmap
