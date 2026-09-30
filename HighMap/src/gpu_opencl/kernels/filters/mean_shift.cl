R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel mean_shift(read_only image2d_t  img_in,
                       write_only image2d_t img_out,
                       const int            nx,
                       const int            ny,
                       const int            ir,
                       const float          talus,
                       const int            talus_weighted)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  float val = 0.f;
  float sum = 0.f;
  float norm = 0.f;
  int   count = 0;
  float v_ij = read_imagef(img_in, sampler, g).x;

  if (talus_weighted == 1)
  {
    // weighted by the difference
    for (int p = g.x - ir; p < g.x + ir + 1; ++p)
      for (int q = g.y - ir; q < g.y + ir + 1; ++q)
      {
        float v_pq = read_imagef(img_in, sampler, (int2)(p, q)).x;
        float dv = fabs(v_ij - v_pq);
        if (dv < talus)
        {
          float weight = 1.f - dv / talus;
          sum += v_pq * weight;
          norm += weight;
        }
      }
    val = sum / norm;
  }
  else
  {
    // basic average
    for (int p = g.x - ir; p < g.x + ir + 1; ++p)
      for (int q = g.y - ir; q < g.y + ir + 1; ++q)
      {
        float v_pq = read_imagef(img_in, sampler, (int2)(p, q)).x;
        float dv = fabs(v_ij - v_pq);
        if (dv < talus)
        {
          sum += v_pq;
          count++;
        }
      }
    val = sum / (float)count;
  }

  write_imagef(img_out, g, val);
}
)""
