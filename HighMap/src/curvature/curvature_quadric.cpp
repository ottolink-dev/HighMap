/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/curvature.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"

namespace hmap
{

Array curvature_quadric(const Array &z, int ir, CurvatureType curvature_type)
{
  if (!validate_non_empty(z)) return Array();

  const glm::ivec2 &shape = z.shape;
  const int         nx = shape.x;
  const int         ny = shape.y;
  const int         type_id = static_cast<int>(curvature_type);

  Array out(shape);

  auto sample = [&](int x, int y) -> float
  {
    int xc = std::clamp(x, 0, nx - 1);
    int yc = std::clamp(y, 0, ny - 1);
    return z(xc, yc);
  };

  if (ir == 0)
  {
    for (int gy = 0; gy < ny; ++gy)
      for (int gx = 0; gx < nx; ++gx)
      {
        const float zc = sample(gx, gy);

        float zx = (sample(gx + 1, gy) - sample(gx - 1, gy)) / 2.f;
        float zy = (sample(gx, gy + 1) - sample(gx, gy - 1)) / 2.f;

        float zxx = sample(gx + 1, gy) + sample(gx - 1, gy) - 2.f * zc;
        float zyy = sample(gx, gy + 1) + sample(gx, gy - 1) - 2.f * zc;

        float zxy = (sample(gx + 1, gy + 1) + sample(gx - 1, gy - 1) -
                     sample(gx + 1, gy - 1) - sample(gx - 1, gy + 1)) /
                    4.f;

        float p = zx * zx + zy * zy;
        float q = p + 1.f;
        float val = 0.f;

        switch (type_id)
        {
        case 0: // MIN
        {
          float H = ((1.f + zy * zy) * zxx - 2.f * zxy * zx * zy +
                     (1.f + zx * zx) * zyy) /
                    (2.f * std::sqrt(q));
          float K = (zxx * zyy - zxy * zxy) / (q * q);
          val = H - std::sqrt(std::fmax(H * H - K, 0.f));
          break;
        }

        case 1: // MAX
        {
          float H = ((1.f + zy * zy) * zxx - 2.f * zxy * zx * zy +
                     (1.f + zx * zx) * zyy) /
                    (2.f * std::sqrt(q));
          float K = (zxx * zyy - zxy * zxy) / (q * q);
          val = H + std::sqrt(std::fmax(H * H - K, 0.f));
          break;
        }

        case 2: // MEAN
          val = ((1.f + zy * zy) * zxx - 2.f * zxy * zx * zy +
                 (1.f + zx * zx) * zyy) /
                (2.f * std::sqrt(q));
          break;

        case 3: // GAUSSIAN
          val = (zxx * zyy - zxy * zxy) / (q * q);
          break;

        case 4: // PROFILE
          val = (p <= 1e-16f)
                    ? 0.f
                    : (zxx * zx * zx + 2.f * zxy * zx * zy + zyy * zy * zy) /
                          (p * std::sqrt(q));
          break;

        case 5: // CONTOUR
          val = (p <= 1e-16f)
                    ? 0.f
                    : (zxx * zy * zy - 2.f * zxy * zx * zy + zyy * zx * zx) /
                          std::sqrt(q);
          break;

        case 6: // TANGENTIAL
          val = (p <= 1e-16f)
                    ? 0.f
                    : (zxx * zy * zy - 2.f * zxy * zx * zy + zyy * zx * zx) /
                          (p * std::sqrt(q));
          break;

        case 7: // ACCUMULATION
        {
          float H = ((1.f + zy * zy) * zxx - 2.f * zxy * zx * zy +
                     (1.f + zx * zx) * zyy) /
                    (2.f * std::sqrt(q));
          float K = (zxx * zyy - zxy * zxy) / (q * q);
          val = H * H - K * K;
          break;
        }

        case 8: // SHAPE INDEX
        {
          float H = ((1.f + zy * zy) * zxx - 2.f * zxy * zx * zy +
                     (1.f + zx * zx) * zyy) /
                    (2.f * std::sqrt(q));
          float K = (zxx * zyy - zxy * zxy) / (q * q);
          float d = std::sqrt(std::fmax(H * H - K, 0.f));
          val = 2.f / float(M_PI) * std::atan(H / (d + 1e-30f));
          val = val * 0.5f + 0.5f;
          break;
        }

        case 9: // UNSPHERICITY
        {
          float H = ((1.f + zy * zy) * zxx - 2.f * zxy * zx * zy +
                     (1.f + zx * zx) * zyy) /
                    (2.f * std::sqrt(q));
          float K = (zxx * zyy - zxy * zxy) / (q * q);
          val = std::sqrt(std::fmax(H * H - K, 0.f));
          break;
        }

        case 10: // RING
        {
          float num = ((zx * zx - zy * zy) * zxy - zx * zy * (zxx - zyy));
          float den = (p + 1e-30f) * (1.f + p);
          float tmp = num / den;
          val = tmp * tmp;
          break;
        }

        case 11: // ROTOR
        {
          float num = ((zx * zx - zy * zy) * zxy - zx * zy * (zxx - zyy));
          val = num / std::pow(p + 1e-6f, 1.5f);
          break;
        }

        default: val = 0.f; break;
        }

        out(gx, gy) = val;
      }

    extrapolate_borders(out);
  }
  else
  {
    for (int gy = 0; gy < ny; ++gy)
      for (int gx = 0; gx < nx; ++gx)
      {
        const float z0 = sample(gx, gy);

        float M00 = 0.f, M11 = 0.f, M01 = 0.f;
        float M22 = 0.f, M33 = 0.f, M44 = 0.f;

        float r0 = 0.f, r1 = 0.f, r2 = 0.f, r3 = 0.f, r4 = 0.f;

        for (int i = -ir; i <= ir; ++i)
          for (int j = -ir; j <= ir; ++j)
          {
            float x = float(i);
            float y = float(j);

            float dz = sample(gx + i, gy + j) - z0;

            float x2 = x * x;
            float y2 = y * y;

            M00 += x2 * x2;
            M11 += y2 * y2;
            M01 += x2 * y2;
            M22 += x2 * y2;
            M33 += x2;
            M44 += y2;

            r0 += dz * x2;
            r1 += dz * y2;
            r2 += dz * x * y;
            r3 += dz * x;
            r4 += dz * y;
          }

        float det = M00 * M11 - M01 * M01;

        float a = (M11 * r0 - M01 * r1) / det;
        float b = (M00 * r1 - M01 * r0) / det;
        float c = r2 / M22;
        float d = r3 / M33;
        float e = r4 / M44;

        float p = d * d + e * e;
        float val = 0.f;

        switch (type_id)
        {
        case 0: val = a + b - std::sqrt((a - b) * (a - b) + c * c); break;

        case 1: val = a + b + std::sqrt((a - b) * (a - b) + c * c); break;

        case 2: val = a + b; break;

        case 3: val = (a + b) * (a + b) - (a - b) * (a - b) + c * c; break;

        case 4:
          val = (p <= 1e-16f) ? 0.f
                              : 2.f * (a * d * d + b * e * e + c * e * d) /
                                    (p * std::sqrt(1.f + p));
          break;

        case 5:
          val = (p <= 1e-16f)
                    ? 0.f
                    : 2.f * (b * d * d + a * e * e - c * d * e) / std::sqrt(p);
          break;

        case 6:
          val = (p <= 1e-16f) ? 0.f
                              : 2.f * (b * d * d + a * e * e - c * d * e) / p;
          break;

        case 7: // ACCUMULATION
        {
          float H = a + b;
          float K = (a + b) * (a + b) - (a - b) * (a - b) + c * c;
          val = H * H - K * K;
          break;
        }

        case 8: // SHAPE INDEX
        {
          float H = a + b;
          float K = (a + b) * (a + b) - (a - b) * (a - b) + c * c;
          float disc = std::sqrt(std::fmax(H * H - K, 0.f));
          val = 2.f / float(M_PI) * std::atan(H / (disc + 1e-30f));
          val = val * 0.5f + 0.5f;
          break;
        }

        case 9: // UNSPHERICITY
        {
          float H = a + b;
          float K = (a + b) * (a + b) - (a - b) * (a - b) + c * c;
          val = std::sqrt(std::fmax(H * H - K, 0.f));
          break;
        }

        case 10: // RING
        {
          float num = ((d * d - e * e) * c - d * e * (2.f * a - 2.f * b));
          float den = (p + 1e-30f) * (1.f + p);
          float tmp = num / den;
          val = tmp * tmp;
          break;
        }

        case 11: // ROTOR
        {
          float num = ((d * d - e * e) * c - d * e * (2.f * a - 2.f * b));
          val = num / std::pow(p + 1e-6f, 1.5f);
          break;
        }

        default: val = 0.f; break;
        }

        out(gx, gy) = val;
      }
  }

  return out;
}

} // namespace hmap

namespace hmap::gpu
{

Array curvature_quadric(const Array  &z,
                        int           ir,
                        CurvatureType curvature_type,
                        bool          approx_algo)
{
  if (!validate_non_empty(z)) return Array();

  const glm::ivec2 shape = z.shape;
  int              type_id = static_cast<int>(curvature_type);

  Array out(shape);

  auto run = clwrapper::Run("curvature_quadric");

  Array zf;

  if (approx_algo)
  {
    zf = z;
    gpu::smooth_cpulse(zf, ir);
    ir = 0;
    run.bind_imagef("z", zf.vector, shape.x, shape.y);
  }
  else
  {
    run.bind_imagef("z", z.vector, shape.x, shape.y);
  }

  run.bind_imagef("out", out.vector, shape.x, shape.y, true);

  run.bind_arguments(shape.x, shape.y, ir, type_id);

  run.execute({shape.x, shape.y});

  run.read_imagef("out");

  if (ir == 0) extrapolate_borders(out);

  return out;
}

} // namespace hmap::gpu
