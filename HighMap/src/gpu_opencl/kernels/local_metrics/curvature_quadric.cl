R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel curvature_quadric(read_only image2d_t  img_in,
                              write_only image2d_t img_out,
                              const int            nx,
                              const int            ny,
                              const int            ir,
                              const int            curv_type)
{
  // Adapted from terrain-descriptors
  // Original author: oargudo
  // License: MIT (see THIRD_PARTY_LICENSES.md)

  const int2 g = {get_global_id(0), get_global_id(1)};
  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  const float dx = 1.f;
  const float dy = 1.f;
  float       val = 0.f;

  /* ============================================================
   * CASE 1: FIXED 3x3 DERIVATIVES (ir == 0)
   * ============================================================ */
  if (ir == 0)
  {
    const float zc = read_imagef(img_in, sampler, g).x;

    float zx = (read_imagef(img_in, sampler, (int2)(g.x + 1, g.y)).x -
                read_imagef(img_in, sampler, (int2)(g.x - 1, g.y)).x) /
               (2.f * dx);

    float zy = (read_imagef(img_in, sampler, (int2)(g.x, g.y + 1)).x -
                read_imagef(img_in, sampler, (int2)(g.x, g.y - 1)).x) /
               (2.f * dy);

    float zxx = (read_imagef(img_in, sampler, (int2)(g.x + 1, g.y)).x +
                 read_imagef(img_in, sampler, (int2)(g.x - 1, g.y)).x -
                 2.f * zc) /
                (dx * dx);

    float zyy = (read_imagef(img_in, sampler, (int2)(g.x, g.y + 1)).x +
                 read_imagef(img_in, sampler, (int2)(g.x, g.y - 1)).x -
                 2.f * zc) /
                (dy * dy);

    float zxy = (read_imagef(img_in, sampler, (int2)(g.x + 1, g.y + 1)).x +
                 read_imagef(img_in, sampler, (int2)(g.x - 1, g.y - 1)).x -
                 read_imagef(img_in, sampler, (int2)(g.x + 1, g.y - 1)).x -
                 read_imagef(img_in, sampler, (int2)(g.x - 1, g.y + 1)).x) /
                (4.f * dx * dy);

    float p = zx * zx + zy * zy;
    float q = p + 1.f;

    switch (curv_type)
    {
    case 0: // MIN
      val = (zc == 0.f) ? 0.f
                        : ((1.f + zy * zy) * zxx - 2.f * zxy * zx * zy +
                           (1.f + zx * zx) * zyy) /
                              (2.f * sqrt(q));
      break;

    case 1: // MAX
      val = (zc == 0.f)
                ? 0.f
                : ((1.f + zy * zy) * zxx - 2.f * zxy * zx * zy +
                   (1.f + zx * zx) * zyy) /
                      (2.f * sqrt(q)); // same base; post-processed if needed
      break;

    case 2: // MEAN
      val = ((1.f + zy * zy) * zxx - 2.f * zxy * zx * zy +
             (1.f + zx * zx) * zyy) /
            (2.f * sqrt(q));
      break;

    case 3: // GAUSSIAN
      val = (zxx * zyy - zxy * zxy) / (q * q);
      break;

    case 4: // PROFILE
      val = (p <= 1e-16f)
                ? 0.f
                : (zxx * zx * zx + 2.f * zxy * zx * zy + zyy * zy * zy) /
                      (p * sqrt(q));
      break;

    case 5: // CONTOUR
      val = (p <= 1e-16f)
                ? 0.f
                : (zxx * zy * zy - 2.f * zxy * zx * zy + zyy * zx * zx) /
                      sqrt(q);
      break;

    case 6: // TANGENTIAL
      val = (p <= 1e-16f)
                ? 0.f
                : (zxx * zy * zy - 2.f * zxy * zx * zy + zyy * zx * zx) /
                      (p * sqrt(q));
      break;

    case 7: // ACCUMULATION
    {
      float H = ((1.f + zy * zy) * zxx - 2.f * zxy * zx * zy +
                 (1.f + zx * zx) * zyy) /
                (2.f * sqrt(q));

      float K = (zxx * zyy - zxy * zxy) / (q * q);

      val = H * H - K * K;
      break;
    }

    case 8: // SHAPE INDEX
    {
      float H = ((1.f + zy * zy) * zxx - 2.f * zxy * zx * zy +
                 (1.f + zx * zx) * zyy) /
                (2.f * sqrt(q));

      float K = (zxx * zyy - zxy * zxy) / (q * q);

      float d = sqrt(fmax(H * H - K, 0.f));

      val = 2.f / M_PI * atan(H / (d + 1e-30f));
      val = val * 0.5f + 0.5f;
      break;
    }

    case 9: // UNSPHERICITY
    {
      float H = ((1.f + zy * zy) * zxx - 2.f * zxy * zx * zy +
                 (1.f + zx * zx) * zyy) /
                (2.f * sqrt(q));

      float K = (zxx * zyy - zxy * zxy) / (q * q);

      val = sqrt(fmax(H * H - K, 0.f));
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

      val = num / pow(p + 1e-6f, 1.5f);
      break;
    }

    default: val = 0.f; break;
    }
  }

  /* ============================================================
   * CASE 2: QUADRIC FIT (ir != 0)
   * ============================================================ */

  else
  {
    const float z0 = read_imagef(img_in, sampler, g).x;

    float M00 = 0.f, M11 = 0.f, M01 = 0.f;
    float M22 = 0.f, M33 = 0.f, M44 = 0.f;

    float r0 = 0.f, r1 = 0.f, r2 = 0.f, r3 = 0.f, r4 = 0.f;

    for (int i = -ir; i <= ir; ++i)
      for (int j = -ir; j <= ir; ++j)
      {
        float x = (float)i * dx;
        float y = (float)j * dy;

        int xi = g.x + i;
        int yj = g.y + j;

        float z = read_imagef(img_in, sampler, (int2)(xi, yj)).x - z0;

        float x2 = x * x;
        float y2 = y * y;

        M00 += x2 * x2;
        M11 += y2 * y2;
        M01 += x2 * y2;
        M22 += x2 * y2;
        M33 += x2;
        M44 += y2;

        r0 += z * x2;
        r1 += z * y2;
        r2 += z * x * y;
        r3 += z * x;
        r4 += z * y;
      }

    float det = M00 * M11 - M01 * M01;

    float a = (M11 * r0 - M01 * r1) / det;
    float b = (M00 * r1 - M01 * r0) / det;
    float c = r2 / M22;
    float d = r3 / M33;
    float e = r4 / M44;

    float p = d * d + e * e;

    switch (curv_type)
    {
    case 0: val = a + b - sqrt((a - b) * (a - b) + c * c); break;

    case 1: val = a + b + sqrt((a - b) * (a - b) + c * c); break;

    case 2: val = a + b; break;

    case 3: val = (a + b) * (a + b) - (a - b) * (a - b) + c * c; break;

    case 4:
      val = (p <= 1e-16f) ? 0.f
                          : 2.f * (a * d * d + b * e * e + c * e * d) /
                                (p * sqrt(1.f + p));
      break;

    case 5:
      val = (p <= 1e-16f) ? 0.f
                          : 2.f * (b * d * d + a * e * e - c * d * e) / sqrt(p);
      break;

    case 6:
      val = (p <= 1e-16f) ? 0.f : 2.f * (b * d * d + a * e * e - c * d * e) / p;
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

      float d = sqrt(fmax(H * H - K, 0.f));

      val = 2.f / M_PI * atan(H / (d + 1e-30f));
      val = val * 0.5f + 0.5f;
      break;
    }

    case 9: // UNSPHERICITY
    {
      float H = a + b;

      float K = (a + b) * (a + b) - (a - b) * (a - b) + c * c;

      val = sqrt(fmax(H * H - K, 0.f));
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

      val = num / pow(p + 1e-6f, 1.5f);
      break;
    }

    default: val = 0.f; break;
    }
  }

  write_imagef(img_out, g, val);
}
)""
