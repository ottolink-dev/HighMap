R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel water_depth_filter(read_only image2d_t  depth,
                               read_only image2d_t  water_mask,
                               read_only image2d_t  zt,
                               write_only image2d_t zt_out,
                               const int            nx,
                               const int            ny,
                               const int            ir)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  const int i = g.x;
  const int j = g.y;

  if (TGET(water_mask, i, j) == 0.f) return;

  float val = 0.f;
  float sum = 0.f;

  for (int r = -ir; r <= ir; r++)
    for (int s = -ir; s <= ir; s++)
    {
      if (is_within_radius_and_inside(i, j, g.x, g.y, ir, nx, ny))
      {
        float current = TGET(depth, i + r, j + s);

        // only taken into account "water" cells
        if (current > 0.f)
        {
          float r = 1.f; // max(0.f, 1.f - hypot(r, s) / ir);
          val += r * current;
          sum += r;
        }
      }
    }

  if (sum > 0.f) TSET(zt_out, i, j, val / sum);
}
)""
