R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

float helper_point_to_line_squared_distance(float x0,
                                            float y0,
                                            float x1,
                                            float y1,
                                            float x2,
                                            float y2)
{
  float dx = x2 - x1;
  float dy = y2 - y1;
  float num = dx * (y1 - y0) - (x1 - x0) * dy;
  return (num * num) / (dx * dx + dy * dy);
}

float2 helper_closest_point_on_line(float2 p, float2 a, float2 b)
{
  float2 ab = b - a;
  float2 ap = p - a;
  float  t = dot(ap, ab) / dot(ab, ab);
  return a + ab * t;
}

void kernel vorolines(global float   *output,
                      global float   *noise_x,
                      global float   *noise_y,
                      constant float *xp,
                      constant float *yp,
                      constant float *xs,
                      constant float *ys,
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

    for (int k = 0; k < npoints; ++k)
    {
      float dist = helper_point_to_line_squared_distance(pos.x,
                                                         pos.y,
                                                         xp[k],
                                                         yp[k],
                                                         xs[k],
                                                         ys[k]);
      min1 = smin(min1, dist, k_smoothing);
    }
    val = min1;
  }
  break;
    //
  case 1: // F2
  {
    float min1 = FLT_MAX;
    float min2 = FLT_MAX;

    for (int k = 0; k < npoints; ++k)
    {
      float dist = helper_point_to_line_squared_distance(pos.x,
                                                         pos.y,
                                                         xp[k],
                                                         yp[k],
                                                         xs[k],
                                                         ys[k]);
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

    for (int k = 0; k < npoints; ++k)
    {
      float dist = helper_point_to_line_squared_distance(pos.x,
                                                         pos.y,
                                                         xp[k],
                                                         yp[k],
                                                         xs[k],
                                                         ys[k]);
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

    for (int k = 0; k < npoints; ++k)
    {
      float dist = helper_point_to_line_squared_distance(pos.x,
                                                         pos.y,
                                                         xp[k],
                                                         yp[k],
                                                         xs[k],
                                                         ys[k]);
      float new_min1 = smin(min1, dist, k_smoothing);
      float new_min2 = smin(min2, smax(min1, dist, k_smoothing), k_smoothing);
      min1 = new_min1;
      min2 = new_min2;
    }
    val = min1 / max(1e-2f, min2);
  }
  break;
    //
  case 4: // F2 - F1
  {
    float min1 = FLT_MAX;
    float min2 = FLT_MAX;

    for (int k = 0; k < npoints; ++k)
    {
      float dist = helper_point_to_line_squared_distance(pos.x,
                                                         pos.y,
                                                         xp[k],
                                                         yp[k],
                                                         xs[k],
                                                         ys[k]);
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
    float2 diff_min;

    for (int k = 0; k < npoints; ++k)
    {
      float dist = helper_point_to_line_squared_distance(pos.x,
                                                         pos.y,
                                                         xp[k],
                                                         yp[k],
                                                         xs[k],
                                                         ys[k]);

      if (dist < min1)
      {
        min1 = dist;
        float2 p = helper_closest_point_on_line(pos,
                                                (float2)(xp[k], yp[k]),
                                                (float2)(xs[k], ys[k]));
        diff_min = p - pos;
      }
    }

    val = FLT_MAX;

    for (int k = 0; k < npoints; ++k)
    {
      float2 p = helper_closest_point_on_line(pos,
                                              (float2)(xp[k], yp[k]),
                                              (float2)(xs[k], ys[k]));
      float2 diff = p - pos;

      if (dot(diff - diff_min, diff - diff_min) > 1e-9f)
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
    float2 diff_min;

    for (int k = 0; k < npoints; ++k)
    {
      float dist = helper_point_to_line_squared_distance(pos.x,
                                                         pos.y,
                                                         xp[k],
                                                         yp[k],
                                                         xs[k],
                                                         ys[k]);

      if (dist < min1)
      {
        min1 = dist;
        float2 p = helper_closest_point_on_line(pos,
                                                (float2)(xp[k], yp[k]),
                                                (float2)(xs[k], ys[k]));
        diff_min = p - pos;
      }
    }

    val = FLT_MAX;

    for (int k = 0; k < npoints; ++k)
    {
      float2 p = helper_closest_point_on_line(pos,
                                              (float2)(xp[k], yp[k]),
                                              (float2)(xs[k], ys[k]));
      float2 diff = p - pos;

      if (dot(diff - diff_min, diff - diff_min) > 1e-9f)
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

    for (int k = 0; k < npoints; ++k)
    {
      float dist = helper_point_to_line_squared_distance(pos.x,
                                                         pos.y,
                                                         xp[k],
                                                         yp[k],
                                                         xs[k],
                                                         ys[k]);
      float edge_value = (float)k;

      if (k_smoothing > 1e-6f)
      {
        float h = smoothstep(-1.f, 1.f, (min1 - dist) / k_smoothing);
        val = lerp(val, edge_value, h) -
              h * (1.f - h) * k_smoothing / (1.f + 3.f * k_smoothing);
      }
      else
      {
        if (dist < min1) val = edge_value;
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

    for (int k = 0; k < npoints; ++k)
    {
      float dist = helper_point_to_line_squared_distance(pos.x,
                                                         pos.y,
                                                         xp[k],
                                                         yp[k],
                                                         xs[k],
                                                         ys[k]);
      float edge_value = (float)k;
      ;

      if (k_smoothing > 1e-6f)
      {
        float h = smoothstep(-1.f, 1.f, (min1 - dist) / k_smoothing);
        val = lerp(val, edge_value, h) -
              h * (1.f - h) * k_smoothing / (1.f + 3.f * k_smoothing);
      }
      else
      {
        if (dist < min1) val = edge_value;
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
