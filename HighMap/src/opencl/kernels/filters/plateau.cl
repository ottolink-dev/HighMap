R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel plateau_post(global float *array,
                         global float *amin,
                         global float *amax,
                         global float *mask,
                         const int     nx,
                         const int     ny,
                         const float   factor,
                         const int     has_mask)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int index = linear_index(g.x, g.y, nx);

  float val = array[index];
  float val0 = val;
  float amin_v = amin[index];
  float amax_v = amax[index];

  // to [0, 1]
  val = (val - amin_v) / (amax_v - amin_v + 1e-6);
  val = max(0.f, val);

  // gain
  val = val < 0.5f ? 0.5f * pow_float(2.f * val, factor)
                   : 1.f - 0.5f * pow_float(2.f * (1.f - val), factor);

  // revert to orginal amplitude
  val = amin_v + (amax_v - amin_v) * val;

  // mask
  if (has_mask)
  {
    float mask_v = mask[index];
    val = (1.f - mask_v) * val0 + mask_v * val;
  }

  array[index] = val;
}
)""
