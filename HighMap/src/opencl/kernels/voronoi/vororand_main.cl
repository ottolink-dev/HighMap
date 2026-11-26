R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

void kernel vororand(global float   *output,
                     global float   *noise_x,
                     global float   *noise_y,
                     constant float *xp,
                     constant float *yp,
                     const int       nx,
                     const int       ny,
                     const int       npoints,
                     const float     k_smoothing,
                     const float     exp_sigma,
                     const int       return_type,
                     const int       has_noise_x,
                     const int       has_noise_y,
                     const float4    bbox)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  int    index = linear_index(g.x, g.y, nx);
  float  dx = has_noise_x > 0 ? noise_x[index] : 0.f;
  float  dy = has_noise_y > 0 ? noise_y[index] : 0.f;
  float2 pos = g_to_xy(g, nx, ny, 1.f, 1.f, dx, dy, bbox);

  float val = 0.f;

  switch (return_type)
  {
  case 0: // F1
  {
    float min1 = FLT_MAX;

    for (int i = 0; i < npoints; ++i)
    {
      float2 p = (float2)(xp[i], yp[i]);
      float2 diff = p - pos;
      float  dist = dot(diff, diff);
      min1 = smin(min1, dist, k_smoothing);
    }
    val = min(10.f, min1);
  }
  break;
  //
  case 1: // F2
  {
    float min1 = FLT_MAX;
    float min2 = FLT_MAX;

    for (int i = 0; i < npoints; ++i)
    {
      float2 p = (float2)(xp[i], yp[i]);
      float2 diff = p - pos;
      float  dist = dot(diff, diff);

      float new_min1 = smin(min1, dist, k_smoothing);
      float new_min2 = smin(min2, smax(min1, dist, k_smoothing), k_smoothing);
      min1 = new_min1;
      min2 = new_min2;
    }
    val = min2;
  }
  break;
  //
  case 2: // F1 * F2
  {
    float min1 = FLT_MAX;
    float min2 = FLT_MAX;

    for (int i = 0; i < npoints; ++i)
    {
      float2 p = (float2)(xp[i], yp[i]);
      float2 diff = p - pos;
      float  dist = dot(diff, diff);

      float new_min1 = smin(min1, dist, k_smoothing);
      float new_min2 = smin(min2, smax(min1, dist, k_smoothing), k_smoothing);
      min1 = new_min1;
      min2 = new_min2;
    }
    val = min1 * min2;
  }
  break;
  //
  case 3: // F1 / F2
  {
    float min1 = FLT_MAX;
    float min2 = FLT_MAX;

    for (int i = 0; i < npoints; ++i)
    {
      float2 p = (float2)(xp[i], yp[i]);
      float2 diff = p - pos;
      float  dist = dot(diff, diff);

      float new_min1 = smin(min1, dist, k_smoothing);
      float new_min2 = smin(min2, smax(min1, dist, k_smoothing), k_smoothing);
      min1 = new_min1;
      min2 = new_min2;
    }
    val = min1 / min2;
  }
  break;
  //
  case 4: // F2 - F1
  {
    float min1 = FLT_MAX;
    float min2 = FLT_MAX;

    for (int i = 0; i < npoints; ++i)
    {
      float2 p = (float2)(xp[i], yp[i]);
      float2 diff = p - pos;
      float  dist = dot(diff, diff);

      float new_min1 = smin(min1, dist, k_smoothing);
      float new_min2 = smin(min2, smax(min1, dist, k_smoothing), k_smoothing);
      min1 = new_min1;
      min2 = new_min2;
    }
    val = min2 - min1;
  }
  break;
    //
  case 5: // edge distance exp.
  {
    float  min1 = FLT_MAX;
    float2 diff_min = (float2)(0.f, 0.f);

    for (int i = 0; i < npoints; ++i)
    {
      float2 p = (float2)(xp[i], yp[i]);
      float2 diff = p - pos;
      float  dist = dot(diff, diff);

      if (dist < min1)
      {
        min1 = dist;
        diff_min = diff;
      }
    }

    val = FLT_MAX;

    for (int i = 0; i < npoints; ++i)
    {
      float2 p = (float2)(xp[i], yp[i]);
      float2 diff = p - pos;

      if (dot(diff - diff_min, diff - diff_min) > 1e-5f)
      {
        float dist = dot(0.5f * (diff_min + diff), normalize(diff - diff_min));
        val = smin(val, dist, k_smoothing);
      }
    }

    val = exp(-0.5f * val * val / (exp_sigma * exp_sigma));
  }
  break;
    //
  case 6: // edge distance
  {
    float  min1 = FLT_MAX;
    float2 diff_min = (float2)(0.f, 0.f);

    for (int i = 0; i < npoints; ++i)
    {
      float2 p = (float2)(xp[i], yp[i]);
      float2 diff = p - pos;
      float  dist = dot(diff, diff);

      if (dist < min1)
      {
        min1 = dist;
        diff_min = diff;
      }
    }

    val = FLT_MAX;

    for (int i = 0; i < npoints; ++i)
    {
      float2 p = (float2)(xp[i], yp[i]);
      float2 diff = p - pos;

      if (dot(diff - diff_min, diff - diff_min) > 1e-5f)
      {
        float dist = dot(0.5f * (diff_min + diff), normalize(diff - diff_min));
        val = smin(val, dist, k_smoothing);
      }
    }
  }
  break;
    //
  case 7: // constant
  {
    float min1 = FLT_MAX;

    for (int i = 0; i < npoints; ++i)
    {
      float2 p = (float2)(xp[i], yp[i]);
      float2 diff = p - pos;
      float  dist = dot(diff, diff);

      if (k_smoothing > 1e-6f)
      {
        float h = smoothstep(-1.f, 1.f, (min1 - dist) / k_smoothing);
        val = lerp(val, (float)i, h) -
              h * (1.f - h) * k_smoothing / (1.f + 3.f * k_smoothing);
      }
      else
      {
        if (dist < min1) val = (float)i;
      }

      min1 = min(dist, min1);
    }
  }
  break;
  //
  case 8: // constant * (F2 - F1)
  {
    float min1 = FLT_MAX;
    float min2 = FLT_MAX;

    for (int i = 0; i < npoints; ++i)
    {
      float2 p = (float2)(xp[i], yp[i]);
      float2 diff = p - pos;
      float  dist = dot(diff, diff);

      if (k_smoothing > 1e-6f)
      {
        float h = smoothstep(-1.f, 1.f, (min1 - dist) / k_smoothing);
        val = lerp(val, (float)i, h) -
              h * (1.f - h) * k_smoothing / (1.f + 3.f * k_smoothing);
      }
      else
      {
        if (dist < min1) val = (float)i;
      }

      float new_min1 = smin(min1, dist, k_smoothing);
      float new_min2 = smin(min2, smax(min1, dist, k_smoothing), k_smoothing);
      min1 = new_min1;
      min2 = new_min2;
    }

    val *= (min2 - min1);
  }
  break;
  }

  output[index] = val;
}
)""
