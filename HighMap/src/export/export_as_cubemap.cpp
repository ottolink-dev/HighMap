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
                       Cmap               cmap,
                       bool               splitted,
                       Array             *p_cubemap)
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

  // 3D partition-of-unity blending: h(p) depends only on the 3D cube point p,
  // guaranteeing C1/C2 continuity across all 12 edges and 8 triple-corners.
  for (int f = 0; f < 6; f++)
  {
    const CubemapFace &face = faces[f];
    Array             &out = out_faces[f];

    for (int j = 0; j < res; j++)
    {
      float v = -1.f + static_cast<float>(j) * inv_res_m1;

      for (int i = 0; i < res; i++)
      {
        float     u = -1.f + static_cast<float>(i) * inv_res_m1;
        glm::vec3 p = face.n + u * face.t + v * face.b;

        float sum = 0.f;
        float wsum = 0.f;

        for (int g = 0; g < 6; g++)
        {
          const CubemapFace &cand = faces[g];
          float              pn = glm::dot(p, cand.n);

          if (pn <= 0.f) continue;

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
        }

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

  if (p_cubemap)
  {
    *p_cubemap = Array(glm::ivec2(4 * res, 3 * res));
    p_cubemap->set_slice(glm::ivec4(0, res, res, 2 * res), z_left);
    p_cubemap->set_slice(glm::ivec4(res, 2 * res, res, 2 * res), z_front);
    p_cubemap->set_slice(glm::ivec4(2 * res, 3 * res, res, 2 * res), z_right);
    p_cubemap->set_slice(glm::ivec4(3 * res, 4 * res, res, 2 * res), z_back);
    p_cubemap->set_slice(glm::ivec4(res, 2 * res, 0, res), z_bottom);
    p_cubemap->set_slice(glm::ivec4(res, 2 * res, 2 * res, 3 * res), z_top);
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
  else
  {
    Array cross(glm::ivec2(4 * res, 3 * res));
    cross.set_slice(glm::ivec4(0, res, res, 2 * res), z_left);
    cross.set_slice(glm::ivec4(res, 2 * res, res, 2 * res), z_front);
    cross.set_slice(glm::ivec4(2 * res, 3 * res, res, 2 * res), z_right);
    cross.set_slice(glm::ivec4(3 * res, 4 * res, res, 2 * res), z_back);
    cross.set_slice(glm::ivec4(res, 2 * res, 0, res), z_bottom);
    cross.set_slice(glm::ivec4(res, 2 * res, 2 * res, 3 * res), z_top);
    cross.to_png(fname, cmap, false);
  }
}

} // namespace hmap
