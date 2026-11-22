/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <string>

#include <glm/glm.hpp>

#include "highmap/array.hpp"
#include "highmap/colormaps.hpp"
#include "highmap/export/image.hpp"
#include "highmap/internal/string_utils.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"
#include "highmap/math/core.hpp"
#include "highmap/transform.hpp"

namespace hmap
{

namespace
{

struct CubemapFace
{
  glm::vec3 n;
  glm::vec3 t;
  glm::vec3 b;
  int       i0;
  int       j0;
};

inline float pou_weight(float t, float overlap)
{
  float a = std::abs(t);
  float inner = 1.f - overlap;
  float outer = 1.f + overlap;

  if (a <= inner) return 1.f;
  if (a >= outer) return 0.f;

  float s = (a - inner) / (outer - inner);
  return 1.f - smoothstep3(s);
}

inline float sample_bilinear(const Array &arr, float x, float y)
{
  int   x0 = std::clamp(static_cast<int>(std::floor(x)), 0, arr.shape.x - 2);
  int   y0 = std::clamp(static_cast<int>(std::floor(y)), 0, arr.shape.y - 2);
  float dx = x - static_cast<float>(x0);
  float dy = y - static_cast<float>(y0);

  float v00 = arr(x0, y0);
  float v10 = arr(x0 + 1, y0);
  float v01 = arr(x0, y0 + 1);
  float v11 = arr(x0 + 1, y0 + 1);

  float v0 = v00 + dx * (v10 - v00);
  float v1 = v01 + dx * (v11 - v01);
  return v0 + dy * (v1 - v0);
}

} // namespace

//----------------------------------------------------------------------
// Main operator
//----------------------------------------------------------------------

void export_as_cubemap(const std::string &fname,
                       const Array       &z,
                       int                cubemap_resolution,
                       float              overlap,
                       int /* ir */,
                       Cmap   cmap,
                       bool   splitted,
                       Array *p_cubemap)
{
  if (!validate_non_empty(z)) return;
  if (!validate_cmap(cmap)) return;

  const int        res = cubemap_resolution;
  const glm::ivec2 cm_shape(res, res);
  const int        noverlap = static_cast<int>(overlap * res);

  hmap::log::trace("cubemap_resolution: {}", res);
  hmap::log::trace("noverlap: {}", noverlap);
  hmap::log::trace("overlap: {}", overlap);

  glm::ivec2 work_shape(4 * res + 2 * noverlap, 3 * res + 2 * noverlap);
  Array      zfull = z.resample_to_shape_bicubic(work_shape);

  // face order: 0=right, 1=left, 2=top, 3=bottom, 4=front, 5=back
  const CubemapFace faces[6] = {
      {{1.f, 0.f, 0.f},
       {0.f, 0.f, -1.f},
       {0.f, 1.f, 0.f},
       noverlap + 2 * res,
       noverlap + res}, // right
      {{-1.f, 0.f, 0.f},
       {0.f, 0.f, 1.f},
       {0.f, 1.f, 0.f},
       noverlap,
       noverlap + res}, // left
      {{0.f, 1.f, 0.f},
       {1.f, 0.f, 0.f},
       {0.f, 0.f, -1.f},
       noverlap + res,
       noverlap + 2 * res}, // top
      {{0.f, -1.f, 0.f},
       {1.f, 0.f, 0.f},
       {0.f, 0.f, 1.f},
       noverlap + res,
       noverlap}, // bottom
      {{0.f, 0.f, 1.f},
       {1.f, 0.f, 0.f},
       {0.f, 1.f, 0.f},
       noverlap + res,
       noverlap + res}, // front
      {{0.f, 0.f, -1.f},
       {-1.f, 0.f, 0.f},
       {0.f, 1.f, 0.f},
       noverlap + 3 * res,
       noverlap + res}, // back
  };

  Array out_faces[6] = {
      Array(cm_shape),
      Array(cm_shape),
      Array(cm_shape),
      Array(cm_shape),
      Array(cm_shape),
      Array(cm_shape),
  };

  const float inv_res_m1 = res > 1 ? 2.f / static_cast<float>(res - 1) : 0.f;
  const float half_res = 0.5f * static_cast<float>(res);
  const float tau = 1.f / (1.f + overlap);

  // precompute 1D coordinates and self-weights
  std::vector<float> coord(res);
  std::vector<float> w_self(res);
  for (int k = 0; k < res; k++)
  {
    coord[k] = -1.f + static_cast<float>(k) * inv_res_m1;
    w_self[k] = pou_weight(coord[k], overlap);
  }

  // adjacent face indices: {+t, -t, +b, -b}
  const int neighbors[6][4] = {
      {5, 4, 2, 3}, // right
      {4, 5, 2, 3}, // left
      {0, 1, 5, 4}, // top
      {0, 1, 4, 5}, // bottom
      {0, 1, 2, 3}, // front
      {1, 0, 2, 3}, // back
  };

  auto accumulate_face =
      [&](const CubemapFace &cand, const glm::vec3 &p, float &sum, float &wsum)
  {
    float pn = glm::dot(p, cand.n);
    if (pn <= 0.f) return;

    float inv_pn = 1.f / pn;
    float ug = glm::dot(p, cand.t) * inv_pn;
    float vg = glm::dot(p, cand.b) * inv_pn;

    float w = pou_weight(ug, overlap) * pou_weight(vg, overlap);
    if (w > 0.f)
    {
      float sx = static_cast<float>(cand.i0) + (ug + 1.f) * half_res;
      float sy = static_cast<float>(cand.j0) + (vg + 1.f) * half_res;
      sum += w * sample_bilinear(zfull, sx, sy);
      wsum += w;
    }
  };

  // 3D partition-of-unity blending: h(p) depends only on the 3D cube point p,
  // guaranteeing C1/C2 continuity across all 12 edges and 8 triple-corners.
  for (int f = 0; f < 6; f++)
  {
    const CubemapFace &face = faces[f];
    Array             &out = out_faces[f];

    for (int j = 0; j < res; j++)
    {
      float v = coord[j];
      float wv = w_self[j];
      bool  v_border = std::abs(v) > tau;
      int   nb_v = v > 0.f ? neighbors[f][2] : neighbors[f][3];

      for (int i = 0; i < res; i++)
      {
        float u = coord[i];
        float sx0 = static_cast<float>(face.i0) + (u + 1.f) * half_res;
        float sy0 = static_cast<float>(face.j0) + (v + 1.f) * half_res;
        float val0 = sample_bilinear(zfull, sx0, sy0);

        bool u_border = std::abs(u) > tau;
        if (!u_border && !v_border)
        {
          out(i, j) = val0;
          continue;
        }

        glm::vec3 p = face.n + u * face.t + v * face.b;
        float     wsum = w_self[i] * wv;
        float     sum = wsum * val0;

        if (u_border)
        {
          int nb_u = u > 0.f ? neighbors[f][0] : neighbors[f][1];
          accumulate_face(faces[nb_u], p, sum, wsum);
        }
        if (v_border) accumulate_face(faces[nb_v], p, sum, wsum);

        out(i, j) = sum / wsum;
      }
    }
  }

  const Array &z_right = out_faces[0];
  const Array &z_left = out_faces[1];
  const Array &z_top = out_faces[2];
  const Array &z_bottom = out_faces[3];
  const Array &z_front = out_faces[4];
  const Array &z_back = out_faces[5];

  if (p_cubemap || !splitted)
  {
    Array cross(glm::ivec2(4 * res, 3 * res));
    cross.set_slice(glm::ivec4(0, res, res, 2 * res), z_left);
    cross.set_slice(glm::ivec4(res, 2 * res, res, 2 * res), z_front);
    cross.set_slice(glm::ivec4(2 * res, 3 * res, res, 2 * res), z_right);
    cross.set_slice(glm::ivec4(3 * res, 4 * res, res, 2 * res), z_back);
    cross.set_slice(glm::ivec4(res, 2 * res, 0, res), z_bottom);
    cross.set_slice(glm::ivec4(res, 2 * res, 2 * res, 3 * res), z_top);

    if (p_cubemap) *p_cubemap = cross;
    if (!splitted) cross.to_png(fname, cmap, false);
  }

  if (splitted)
  {
    std::filesystem::path path(fname);

    z_right.to_png(add_filename_suffix(path, "_00_right").string(), cmap);
    z_left.to_png(add_filename_suffix(path, "_01_left").string(), cmap);
    z_top.to_png(add_filename_suffix(path, "_02_top").string(), cmap);
    z_bottom.to_png(add_filename_suffix(path, "_03_bottom").string(), cmap);
    z_front.to_png(add_filename_suffix(path, "_04_front").string(), cmap);
    z_back.to_png(add_filename_suffix(path, "_05_back").string(), cmap);
  }
}

} // namespace hmap
