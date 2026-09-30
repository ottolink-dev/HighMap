R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel directional_blur(read_only image2d_t  in,
                             read_only image2d_t  angle,
                             write_only image2d_t out,
                             const int            nx,
                             const int            ny,
                             const float          radius,
                             const int            steps)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

  float val = 0.f;
  float sum = 0.f;
  float dr = radius / steps;

  float  alpha = M_PI / 180.f * read_imagef(angle, sampler, g).x;
  float2 dir = {cos(alpha), sin(alpha)};

  for (float r = 0.f; r <= radius; r += dr)
  {
    float  w = 1.f - smoothstep3(r / radius);
    float2 pos = (float2)(g.x + r * dir.x, g.y + r * dir.y);

    val += w * read_imagef(in, sampler, pos).x;
    sum += w;
  }

  if (sum > 1e-8f)
    val /= sum;
  else
    val = 0.f;

  write_imagef(out, g, val);
}
)""
