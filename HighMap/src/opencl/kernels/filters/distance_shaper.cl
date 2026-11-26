R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel distance_shaper_rt(read_only image2d_t  in,
                               write_only image2d_t out,
                               global const float  *thetas,
                               global const float  *slopes,
                               const int            nx,
                               const int            ny,
                               const int            directions,
                               const float          threshold,
                               const float          z_max,
                               const int            mode)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

  float val_orig = read_imagef(in, sampler, g).x;

  float       best_val = val_orig;
  float       second_val = -1e20f;
  const float max_dist = 1.41421356f * (float)max(nx, ny);

  for (int k = 0; k < directions; ++k)
  {
    float th = thetas[k];
    float s = slopes[k];

    float2 dir = (float2)(cos(th), sin(th));
    float  dir_max = val_orig;

    for (float r = 1.0f; r <= max_dist; r += 1.0f)
    {
      float2 pos = (float2)(g.x - r * dir.x, g.y - r * dir.y);

      if (pos.x < 0.0f || pos.x >= (float)nx || pos.y < 0.0f ||
          pos.y >= (float)ny)
      {
        break;
      }

      float h_up = read_imagef(in, sampler, pos).x;

      if (h_up >= threshold)
      {
        float cand = h_up - r * s;

        if (cand > dir_max)
        {
          dir_max = cand;
        }
      }

      // early out: upper bound cannot exceed dir_max
      if (z_max - r * s <= dir_max)
      {
        break;
      }
    }

    if (dir_max > best_val)
    {
      second_val = best_val;
      best_val = dir_max;
    }
    else if (dir_max > second_val)
    {
      second_val = dir_max;
    }
  }

  float res = best_val;

  // mode 0 is ASTERISK
  if (mode == 0)
  {
    if (second_val > -1e19f && best_val > threshold)
    {
      float num = second_val - threshold;
      float den = best_val - threshold;

      if (den > 1e-5f)
      {
        float ratio = clamp(num / den, 0.0f, 1.0f);
        float cut = 0.35f * ratio * ratio * (second_val - threshold);
        res = max(val_orig, best_val - cut);
      }
    }
  }

  write_imagef(out, g, res);
}
)""
