R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
inline int is_fg(float v)
{
  return v > 0.5f;
}

float sample(read_only image2d_t in, sampler_t s, int x, int y, int nx, int ny)
{
  if (x < 0 || x >= nx || y < 0 || y >= ny) return 0.f;
  return read_imagef(in, s, (int2)(x, y)).x;
}

void kernel thinning(read_only image2d_t  in,
                     write_only image2d_t out,
                     int                  nx,
                     int                  ny,
                     int                  iter)
{
  int i = get_global_id(0);
  int j = get_global_id(1);

  if (i >= nx || j >= ny) return;

  sampler_t s = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE |
                CLK_FILTER_NEAREST;

  if (i == 0 || j == 0 || i == nx - 1 || j == ny - 1)
  {
    float val = read_imagef(in, s, (int2)(i, j)).x;
    write_imagef(out, (int2)(i, j), val);
    return;
  }

  float p2 = sample(in, s, i, j + 1, nx, ny);
  float p3 = sample(in, s, i + 1, j + 1, nx, ny);
  float p4 = sample(in, s, i + 1, j, nx, ny);
  float p5 = sample(in, s, i + 1, j - 1, nx, ny);
  float p6 = sample(in, s, i, j - 1, nx, ny);
  float p7 = sample(in, s, i - 1, j - 1, nx, ny);
  float p8 = sample(in, s, i - 1, j, nx, ny);
  float p9 = sample(in, s, i - 1, j + 1, nx, ny);

  int a = (!is_fg(p8) && is_fg(p9)) + (!is_fg(p9) && is_fg(p2)) +
          (!is_fg(p2) && is_fg(p3)) + (!is_fg(p3) && is_fg(p4)) +
          (!is_fg(p4) && is_fg(p5)) + (!is_fg(p5) && is_fg(p6)) +
          (!is_fg(p6) && is_fg(p7)) + (!is_fg(p7) && is_fg(p8));

  int b = is_fg(p2) + is_fg(p3) + is_fg(p4) + is_fg(p5) + is_fg(p6) +
          is_fg(p7) + is_fg(p8) + is_fg(p9);

  int m1 = iter == 0 ? (is_fg(p8) * is_fg(p2) * is_fg(p4))
                     : (is_fg(p8) * is_fg(p2) * is_fg(p6));

  int m2 = iter == 0 ? (is_fg(p2) * is_fg(p4) * is_fg(p6))
                     : (is_fg(p8) * is_fg(p4) * is_fg(p6));

  float val = read_imagef(in, s, (int2)(i, j)).x;

  float marker = (a == 1 && (b >= 2 && b <= 6) && m1 == 0 && m2 == 0) ? 1.f
                                                                      : 0.f;

  write_imagef(out, (int2)(i, j), val * (1.f - marker));
}

void kernel relative_distance_from_skeleton(read_only image2d_t  array,
                                            read_only image2d_t  sk,
                                            read_only image2d_t  border,
                                            write_only image2d_t rdist,
                                            const int            nx,
                                            const int            ny,
                                            const int            ir)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  if (read_imagef(array, sampler, g).x == 0.f) return;

  float dmin_sk = FLT_MAX;
  float dmin_bd = FLT_MAX;

  for (int p = g.x - ir; p <= g.x + ir; p++)
  {
    const int dp2 = (g.x - p) * (g.x - p); // ← hoisted out of inner loop
    for (int q = g.y - ir; q <= g.y + ir; q++)
    {
      const float vsk = read_imagef(sk, sampler, (int2)(p, q)).x;
      const float vbd = read_imagef(border, sampler, (int2)(p, q)).x;

      // skip non-feature cells early
      if (vsk != 1.f && vbd != 1.f) continue;

      const int d2 = dp2 + (g.y - q) * (g.y - q);

      if (vsk == 1.f && d2 < dmin_sk) dmin_sk = d2;
      if (vbd == 1.f && d2 < dmin_bd) dmin_bd = d2;
    }
  }

  const float sum = dmin_sk + dmin_bd;
  write_imagef(rdist, g, sum > 0.f ? dmin_bd / sum : 0.f);
}
)""
