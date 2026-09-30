R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel rotate_hmap(read_only image2d_t  in,
                        write_only image2d_t out,
                        const int            nx,
                        const int            ny,
                        const float          angle,
                        const int            zoom_in)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler_in = CLK_NORMALIZED_COORDS_TRUE |
                               CLK_ADDRESS_MIRRORED_REPEAT | CLK_FILTER_LINEAR;

  float2 pos = {(float)g.x / (nx - 1.f), (float)g.y / (ny - 1.f)};

  // rotate
  float2 pos_c = pos - (float2)(0.5f, 0.5f);
  float  alpha = -angle / 180.f * 3.14159f;
  float  ca = cos(alpha);
  float  sa = sin(alpha);

  float zoom = zoom_in == 1 ? 1.f / (fabs(ca) + fabs(sa)) : 1.f;

  pos.x = 0.5f + zoom * (pos_c.x * ca - pos_c.y * sa);
  pos.y = 0.5f + zoom * (pos_c.x * sa + pos_c.y * ca);

  float val = read_imagef(in, sampler_in, pos).x;

  write_imagef(out, g, val);
}
)""
