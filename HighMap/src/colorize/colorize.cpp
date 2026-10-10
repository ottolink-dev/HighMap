/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <type_traits>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/colorize.hpp"
#include "highmap/colormaps.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate1d.hpp"
#include "highmap/logger.hpp"
#include "highmap/math/array.hpp"
#include "highmap/range.hpp"
#include "highmap/shadows.hpp"
#include "highmap/texture.hpp"

#include "mixbox.h"

#include <format>

namespace hmap
{

void apply_hillshade(Texture     &color3,
                     const Array &array,
                     float        vmin,
                     float        vmax,
                     float        exponent)
{
  if (!validate_non_empty(array)) return;
  if (!validate_non_empty(color3)) return;
  if (!validate_same_shape(color3, array)) return;

  // compute and scale hillshading
  Array hs = Array(array.shape, 1.f);
  hs = hillshade(array, 180.f, 45.f, 10.f * array.ptp() / (float)array.shape.y);
  remap(hs, vmin, vmax);

  if (exponent != 1.f) hs = pow(hs, exponent);

  clamp(hs);

  // apply to image
  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
      for (int ch = 0; ch < 3 && ch < color3.num_channels(); ch++)
        color3(i, j, ch) *= hs(i, j);
}

void apply_hillshade(std::vector<uint8_t> &img,
                     const Array          &array,
                     float                 vmin,
                     float                 vmax,
                     float                 exponent,
                     bool                  is_img_rgba)
{
  if (!validate_non_empty(array)) return;
  size_t expected_size = (size_t)array.shape.x * (size_t)array.shape.y *
                         (is_img_rgba ? 4 : 3);
  if (img.size() < expected_size)
  {
    hmap::log::warn("Image buffer size ({}) is smaller than required size ({})",
                    img.size(),
                    expected_size);
    return;
  }

  // compute and scale hillshading
  Array hs = Array(array.shape, 1.f);
  hs = hillshade(array, 180.f, 45.f, 10.f * array.ptp() / (float)array.shape.y);
  remap(hs, vmin, vmax);

  if (exponent != 1.f) hs = pow(hs, exponent);

  clamp(hs);

  // apply to image
  int k = 0;

  if (is_img_rgba)
  {
    for (int j = array.shape.y - 1; j > -1; j--)
      for (int i = 0; i < array.shape.x; i++)
      {
        img[k] = (uint8_t)((float)img[k] * hs(i, j));
        img[k + 1] = (uint8_t)((float)img[k + 1] * hs(i, j));
        img[k + 2] = (uint8_t)((float)img[k + 2] * hs(i, j));
        // skip alpha channel
        k += 4;
      }
  }
  else
  {
    for (int j = array.shape.y - 1; j > -1; j--)
      for (int i = 0; i < array.shape.x; i++)
      {
        img[k] = (uint8_t)((float)img[k] * hs(i, j));
        img[k + 1] = (uint8_t)((float)img[k + 1] * hs(i, j));
        img[k + 2] = (uint8_t)((float)img[k + 2] * hs(i, j));
        k += 3;
      }
  }
}

Texture colorize(const Array &array,
                 float        vmin,
                 float        vmax,
                 int          cmap,
                 bool         hillshading,
                 bool         reverse,
                 const Array *p_noise)
{
  if (!validate_non_empty(array)) return Texture();
  if (p_noise && !validate_same_shape(array, *p_noise)) return Texture();

  // get the colormap and reverse if needed
  const auto colormap_colors = get_colormap_data(cmap);
  if (reverse) std::swap(vmin, vmax);

  // initialize color texture and normalization factors
  const int nc = static_cast<int>(colormap_colors.size());
  glm::vec2 normalization_factors = array.normalization_coeff(vmin, vmax);
  normalization_factors.x *= (nc - 1);
  normalization_factors.y *= (nc - 1);

  Texture color3(array.shape, 3);

  // lambda function to apply colormap
  auto apply_colormap = [&](float value) -> glm::vec3
  {
    int   q = static_cast<int>(value);
    float t = value - q;

    if (q < nc - 1)
      return (1.f - t) * colormap_colors[q] + t * colormap_colors[q + 1];
    else
      return colormap_colors[q];
  };

  // process each pixel
  for (int j = 0; j < array.shape.y; ++j)
    for (int i = 0; i < array.shape.x; ++i)
    {
      float value = array(i, j);
      if (p_noise) value += (*p_noise)(i, j);

      float vnorm = normalization_factors.x * value + normalization_factors.y;

      float normalized_value = std::clamp(vnorm,
                                          0.f,
                                          static_cast<float>(nc - 1));

      glm::vec3 color = apply_colormap(normalized_value);

      // assign color values to the tensor
      color3(i, j, 0) = color.x;
      color3(i, j, 1) = color.y;
      color3(i, j, 2) = color.z;
    }

  // apply hillshading if required
  if (hillshading) apply_hillshade(color3, array);

  return color3;
}

Texture colorize(const Array                  &array,
                 float                         vmin,
                 float                         vmax,
                 const std::vector<float>     &positions,
                 const std::vector<glm::vec3> &colormap_colors,
                 bool                          reverse,
                 const Array                  *p_noise)
{
  if (!validate_non_empty(array)) return Texture();
  if (p_noise && !validate_same_shape(array, *p_noise)) return Texture();
  if (!validate_non_empty(positions, "Colormap positions")) return Texture();
  if (!validate_non_empty(colormap_colors, "Colormap colors")) return Texture();
  if (positions.size() != colormap_colors.size())
  {
    hmap::log::warn(
        "Colormap positions size ({}) does not match colors size ({})",
        positions.size(),
        colormap_colors.size());
    return Texture();
  }

  auto cpos = positions;
  auto colors = colormap_colors;

  if (reverse)
  {
    std::reverse(colors.begin(), colors.end());
    std::reverse(cpos.begin(), cpos.end());

    for (auto &p : cpos)
      p = 1.f - p;
  }

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

  Texture out(array.shape, 3);

  for (int j = 0; j < array.shape.y; ++j)
  {
    for (int i = 0; i < array.shape.x; ++i)
    {
      float v = array(i, j) + (p_noise ? (*p_noise)(i, j) : 0.f);
      v = (v - vmin) / (vmax - vmin);
      v = std::clamp(v, 0.f, 1.f);

      out(i, j, 0) = citp_r(v);
      out(i, j, 1) = citp_g(v);
      out(i, j, 2) = citp_b(v);
    }
  }

  return out;
}

Texture colorize_bivariate(const Array                  &a1,
                           const Array                  &a2,
                           const glm::vec2               range1,
                           const glm::vec2               range2,
                           const std::vector<float>     &positions1,
                           const std::vector<float>     &positions2,
                           const std::vector<glm::vec3> &colormap_colors1,
                           const std::vector<glm::vec3> &colormap_colors2,
                           MixMethod                     method,
                           bool                          reverse1,
                           bool                          reverse2,
                           const Array                  *p_noise1,
                           const Array                  *p_noise2)
{
  if (!validate_non_empty(a1) || !validate_same_shape(a1, a2)) return Texture();
  if (p_noise1 && !validate_same_shape(a1, *p_noise1)) return Texture();
  if (p_noise2 && !validate_same_shape(a1, *p_noise2)) return Texture();
  if (!validate_non_empty(positions1, "Colormap positions 1")) return Texture();
  if (!validate_non_empty(colormap_colors1, "Colormap colors 1"))
    return Texture();
  if (positions1.size() != colormap_colors1.size())
  {
    hmap::log::warn(
        "Colormap 1 positions size ({}) does not match colors size ({})",
        positions1.size(),
        colormap_colors1.size());
    return Texture();
  }
  if (!validate_non_empty(positions2, "Colormap positions 2")) return Texture();
  if (!validate_non_empty(colormap_colors2, "Colormap colors 2"))
    return Texture();
  if (positions2.size() != colormap_colors2.size())
  {
    hmap::log::warn(
        "Colormap 2 positions size ({}) does not match colors size ({})",
        positions2.size(),
        colormap_colors2.size());
    return Texture();
  }

  const glm::ivec2 shape = a1.shape;
  Texture          color3(shape, 3);

  // --- Colormap preparation

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

  // --- Build color interpolator

  auto make_color_interpolator = [](const std::vector<float>     &positions,
                                    const std::vector<glm::vec3> &colors)
  {
    std::vector<float> r;
    std::vector<float> g;
    std::vector<float> b;

    r.reserve(colors.size());
    g.reserve(colors.size());
    b.reserve(colors.size());

    for (const auto &color : colors)
    {
      r.push_back(color[0]);
      g.push_back(color[1]);
      b.push_back(color[2]);
    }

    return std::array<Interpolator1D, 3>{
        Interpolator1D(positions, r, hmap::InterpolationMethod1D::LINEAR),
        Interpolator1D(positions, g, hmap::InterpolationMethod1D::LINEAR),
        Interpolator1D(positions, b, hmap::InterpolationMethod1D::LINEAR)};
  };

  auto [cpos1,
        colors1] = prepare_colormap(positions1, colormap_colors1, reverse1);

  auto [cpos2,
        colors2] = prepare_colormap(positions2, colormap_colors2, reverse2);

  auto citp1 = make_color_interpolator(cpos1, colors1);
  auto citp2 = make_color_interpolator(cpos2, colors2);

  // --- Color mixing

  auto mix_colors = [method](const glm::vec3 &color1, const glm::vec3 &color2)
  {
    switch (method)
    {
    case MixMethod::MM_SQRT_AVG:
    case MixMethod::MM_LINEAR: return 0.5f * (color1 + color2);

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

  // --- Apply colormaps and mix

  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      float v1 = a1(i, j);
      float v2 = a2(i, j);

      if (p_noise1) v1 += (*p_noise1)(i, j);

      if (p_noise2) v2 += (*p_noise2)(i, j);

      v1 = (v1 - range1.x) / (range1.y - range1.x);
      v2 = (v2 - range2.x) / (range2.y - range2.x);

      v1 = std::clamp(v1, 0.f, 1.f);
      v2 = std::clamp(v2, 0.f, 1.f);

      const glm::vec3 color1{citp1[0](v1), citp1[1](v1), citp1[2](v1)};

      const glm::vec3 color2{citp2[0](v2), citp2[1](v2), citp2[2](v2)};

      const glm::vec3 color = mix_colors(color1, color2);

      color3(i, j, 0) = color[0];
      color3(i, j, 1) = color[1];
      color3(i, j, 2) = color[2];
    }
  }

  return color3;
}

Texture colorize_grayscale(const Array &array)
{
  if (!validate_non_empty(array)) return Texture();
  Texture color1 = Texture(array);
  color1.remap();
  return color1;
}

Texture colorize_histogram(const Array &array)
{
  if (!validate_non_empty(array)) return Texture();
  Texture color1 = Texture(array.shape, 1);

  // normalization factors
  float a = 0.f;
  float b = 0.f;
  float vmin = array.min();
  float vmax = array.max();

  if (vmin != vmax)
  {
    a = 1.f / (vmax - vmin) * (float)(array.shape.x - 1);
    b = -vmin / (vmax - vmin) * (float)(array.shape.x - 1);
  }

  // compute histogram
  std::vector<int> hist(array.shape.x);
  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
      hist[(int)(a * array(i, j) + b)] += 1;

  int hmax = 1;
  if (!hist.empty()) hmax = *std::max_element(hist.begin(), hist.end());

  for (auto &v : hist)
    v = (int)((float)v / (float)hmax * (float)(array.shape.y - 1));

  // create histogram image
  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
      if (j < hist[i]) color1(i, j, 0) = 1.f;

  return color1;
}

Texture colorize_slope_height_heatmap(const Array &array, int cmap)
{
  if (!validate_non_empty(array)) return Texture();
  Array dz = gradient_norm(array);

  // normalization factors / 1
  float a1 = 0.f;
  float b1 = 0.f;
  float vmin1 = array.min();
  float vmax1 = array.max();

  if (vmin1 != vmax1)
  {
    a1 = 1.f / (vmax1 - vmin1) * (float)(array.shape.x - 1);
    b1 = -vmin1 / (vmax1 - vmin1) * (float)(array.shape.x - 1);
  }

  // normalization factors / 2
  float a2 = 0.f;
  float b2 = 0.f;
  float vmin2 = dz.min();
  float vmax2 = dz.max();

  if (vmin2 != vmax2)
  {
    a2 = 1.f / (vmax2 - vmin2) * (float)(array.shape.y - 1);
    b2 = -vmin2 / (vmax2 - vmin2) * (float)(array.shape.y - 1);
  }

  // compute 2D histogram
  Array sum = Array(array.shape);

  for (int j = 0; j < array.shape.y; j++)
    for (int i = 0; i < array.shape.x; i++)
    {
      int p = (int)(a1 * array(i, j) + b1);
      int q = (int)(a2 * dz(i, j) + b2);

      sum(p, q) += 1.f;
    }

  bool    hillshading = false;
  Texture col3 = colorize(sum, sum.min(), sum.max(), cmap, hillshading);

  return col3;
}

Texture colorize_vec2(const Array &array1, const Array &array2)
{
  if (!validate_non_empty(array1) || !validate_same_shape(array1, array2))
    return Texture();

  // create image
  Texture col3 = Texture(array1.shape, 3);

  // normalization factors / 1
  float a1 = 0.f;
  float b1 = 0.f;
  float vmin1 = array1.min();
  float vmax1 = array1.max();

  if (vmin1 != vmax1)
  {
    a1 = 1.f / (vmax1 - vmin1);
    b1 = -vmin1 / (vmax1 - vmin1);
  }

  // normalization factors / 2
  float a2 = 0.f;
  float b2 = 0.f;
  float vmin2 = array2.min();
  float vmax2 = array2.max();

  if (vmin2 != vmax2)
  {
    a2 = 1.f / (vmax2 - vmin2);
    b2 = -vmin2 / (vmax2 - vmin2);
  }

  for (int j = 0; j < array1.shape.y; j++)
    for (int i = 0; i < array1.shape.x; i++)
    {
      float u = a1 * array1(i, j) + b1;
      float v = a2 * array2(i, j) + b2;
      float w = u * v * (1.f - u) * (1.f - v);

      col3(i, j, 0) = u;
      col3(i, j, 1) = v;
      col3(i, j, 2) = w;
    }

  return col3;
}

Array luminance(const Texture &tex)
{
  if (!validate_non_empty(tex, 3)) return Array();
  return 0.299f * tex[0] + 0.587f * tex[1] + 0.114f * tex[2];
}

Texture mix(const Texture &tex1, const Texture &tex2, MixMethod method)
{
  if (!validate_non_empty(tex1, 4) || !validate_channels(tex1, 4))
    return Texture();
  if (!validate_non_empty(tex2, 4) || !validate_channels(tex2, 4))
    return Texture();
  if (!validate_same_shape(tex1, tex2)) return Texture();

  Texture out(tex1.shape, 4);

  const Array &a1 = tex1[3];
  const Array &a2 = tex2[3];

  Array t = a2 / (a2 + a1 * (1.f - a2));

  if (method == MM_MIXBOX)
  {
    for (int nch = 0; nch < 3; ++nch)
    {
      out[nch] = Array(tex1.shape, 0.f);
    }

    int size = tex1.shape.x * tex1.shape.y;
    for (int idx = 0; idx < size; ++idx)
    {
      float r1 = tex1[0].vector[idx];
      float g1 = tex1[1].vector[idx];
      float b1 = tex1[2].vector[idx];

      float r2 = tex2[0].vector[idx];
      float g2 = tex2[1].vector[idx];
      float b2 = tex2[2].vector[idx];

      float ratio = t.vector[idx];

      float r, g, b;
      mixbox_lerp_float(r1, g1, b1, r2, g2, b2, ratio, &r, &g, &b);

      out[0].vector[idx] = r;
      out[1].vector[idx] = g;
      out[2].vector[idx] = b;
    }
  }
  else
  {
    for (int nch = 0; nch < 3; ++nch)
    {
      if (method == MM_SQRT_AVG)
      {
        out[nch] = pow((1.f - t) * tex1[nch] * tex1[nch] +
                           t * tex2[nch] * tex2[nch],
                       0.5f);
      }
      else
      {
        out[nch] = lerp(tex1[nch], tex2[nch], t);
      }
    }
  }

  out[3] = a1 + a2 * (1.f - a1);
  return out;
}

Texture mix(const Texture &tex1,
            const Texture &tex2,
            const Array   &mask,
            MixMethod      method,
            float          gain_factor)
{
  if (!validate_non_empty(tex1, 3)) return Texture();
  if (!validate_non_empty(tex2, 3)) return Texture();
  if (!validate_same_shape(tex1, tex2)) return Texture();
  if (!validate_non_empty(mask)) return Texture();
  if (!validate_same_shape(tex1, mask)) return Texture();

  int     num_channels = std::min(tex1.num_channels(), tex2.num_channels());
  Texture out(tex1.shape, num_channels);

  Array t = mask;
  if (gain_factor != 1.f)
  {
    gain(t, gain_factor);
  }

  if (method == MM_MIXBOX)
  {
    for (int nch = 0; nch < 3; ++nch)
    {
      out[nch] = Array(tex1.shape, 0.f);
    }

    int size = tex1.shape.x * tex1.shape.y;
    for (int idx = 0; idx < size; ++idx)
    {
      float r1 = tex1[0].vector[idx];
      float g1 = tex1[1].vector[idx];
      float b1 = tex1[2].vector[idx];

      float r2 = tex2[0].vector[idx];
      float g2 = tex2[1].vector[idx];
      float b2 = tex2[2].vector[idx];

      float ratio = t.vector[idx];

      float r, g, b;
      mixbox_lerp_float(r1, g1, b1, r2, g2, b2, ratio, &r, &g, &b);

      out[0].vector[idx] = r;
      out[1].vector[idx] = g;
      out[2].vector[idx] = b;
    }
  }
  else
  {
    for (int nch = 0; nch < 3; ++nch)
    {
      if (method == MM_SQRT_AVG)
      {
        out[nch] = pow((1.f - t) * tex1[nch] * tex1[nch] +
                           t * tex2[nch] * tex2[nch],
                       0.5f);
      }
      else
      {
        out[nch] = lerp(tex1[nch], tex2[nch], t);
      }
    }
  }

  for (int nch = 3; nch < num_channels; ++nch)
  {
    out[nch] = lerp(tex1[nch], tex2[nch], t);
  }

  return out;
}

Texture mix(const std::vector<const Texture *> &texs, MixMethod method)
{
  if (!validate_non_empty(texs, "Texture list")) return Texture();
  for (const auto *t : texs)
  {
    if (!t)
    {
      hmap::log::warn("Null texture pointer encountered in mix");
      return Texture();
    }
    if (!validate_non_empty(*t, 4) || !validate_channels(*t, 4))
    {
      return Texture();
    }
  }

  Texture out = *texs.front();

  for (size_t k = 1; k < texs.size(); ++k)
  {
    out = mix(out, *(texs[k]), method);
  }

  return out;
}

Texture mix_normal_map(const Texture          &nmap_base,
                       const Texture          &nmap_detail,
                       float                   detail_scaling,
                       NormalMapBlendingMethod blending_method)
{
  if (!validate_non_empty(nmap_base, 3)) return Texture();
  if (!validate_non_empty(nmap_detail, 3)) return Texture();
  if (!validate_same_shape(nmap_base, nmap_detail)) return Texture();

  Texture out(nmap_base.shape, 4);
  // Copy alpha from base (or detail)
  if (nmap_base.num_channels() == 4)
    out[3] = nmap_base[3];
  else
    out[3] = Array(nmap_base.shape, 1.f);

  std::function<glm::vec3(glm::vec3 &, glm::vec3 &)> blending_fct;

  switch (blending_method)
  {
  case NormalMapBlendingMethod::NMAP_LINEAR:
  {
    blending_fct = [](glm::vec3 &n1, glm::vec3 &n2) { return n1 + n2; };
  }
  break;
  case NormalMapBlendingMethod::NMAP_DERIVATIVE:
  {
    blending_fct = [](glm::vec3 &n1, glm::vec3 &n2)
    {
      return glm::vec3(n1.x * n2.z + n2.x * n1.z,
                       n1.y * n2.z + n2.y * n1.z,
                       n1.z * n2.z);
    };
  }
  break;
  case NormalMapBlendingMethod::NMAP_UDN:
  {
    blending_fct = [](glm::vec3 &n1, glm::vec3 &n2)
    { return glm::vec3(n1.x + n2.x, n1.y + n2.y, n1.z); };
  }
  break;
  case NormalMapBlendingMethod::NMAP_UNITY:
  {
    blending_fct = [](glm::vec3 &n1, glm::vec3 &n2)
    {
      glm::vec3 m0 = glm::vec3(n1.z, n1.x, -n1.x);
      glm::vec3 m1 = glm::vec3(n1.x, n1.z, -n1.y);
      glm::vec3 m2 = glm::vec3(n1.x, n1.y, n1.z);

      return glm::vec3(n2.x * m0.x + n2.y * m1.x + n2.z * m2.x,
                       n2.x * m0.y + n2.y * m1.y + n2.z * m2.y,
                       n2.x * m0.z + n2.y * m1.z + n2.z * m2.z);
    };
  }
  break;
  case NormalMapBlendingMethod::NMAP_WHITEOUT:
  default:
  {
    blending_fct = [](glm::vec3 &n1, glm::vec3 &n2)
    { return glm::vec3(n1.x + n2.x, n1.y + n2.y, n1.z * n2.z); };
  }
  }

  for (int j = 0; j < nmap_base.shape.y; j++)
  {
    for (int i = 0; i < nmap_base.shape.x; i++)
    {
      glm::vec3 v111 = glm::vec3(1.f, 1.f, 1.f);
      glm::vec3 n1 = 2.f * glm::vec3(nmap_base(i, j, 0),
                                     nmap_base(i, j, 1),
                                     nmap_base(i, j, 2)) -
                     v111;
      glm::vec3 n2 = 2.f * glm::vec3(nmap_detail(i, j, 0),
                                     nmap_detail(i, j, 1),
                                     nmap_detail(i, j, 2)) -
                     v111;

      n2.x *= detail_scaling;
      n2.y *= detail_scaling;
      n2.z *= detail_scaling;

      glm::vec3 vn = blending_fct(n1, n2);
      vn = glm::normalize(vn);

      out(i, j, 0) = 0.5f * vn.x + 0.5f;
      out(i, j, 1) = 0.5f * vn.y + 0.5f;
      out(i, j, 2) = 0.5f * vn.z + 0.5f;
    }
  }

  return out;
}

} // namespace hmap
