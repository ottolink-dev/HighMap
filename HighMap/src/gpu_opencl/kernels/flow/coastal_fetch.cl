R""(
/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
void kernel coastal_fetch(read_only image2d_t  z,
                          read_only image2d_t  compute_mask,
                          write_only image2d_t out,
                          const int            nx,
                          const int            ny,
                          const int            ndirections)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

  // early skip if requested
  if (read_imagef(compute_mask, sampler, g).x <= 0.f)
  {
    write_imagef(out, g, 0.f);
    return;
  }

  // full domain diagonal in worst case scenario
  const int    nsteps = (int)(1.414f * max(nx, ny));
  const float2 pos0 = (float2)(g.x + 0.5f, g.y + 0.5f);
  const float  z0 = read_imagef(z, sampler, g).x;

  float fetch = 0.f;
  int   nhit = 0;

  for (int idir = 0; idir < ndirections; ++idir)
  {
    float  theta = 2.f * M_PI * (float)idir / (float)ndirections;
    float2 dir = (float2)(cos(theta), sin(theta));

    for (int k = 0; k < nsteps; ++k)
    {
      float2 ray = pos0 + (float)(k + 1) * dir;

      if (!is_inside((int)ray.x, (int)ray.y, nx + 1, ny + 1))
      {
        fetch += (float)nsteps; // open boundary = full fetch
        nhit++;
        break;
      }

      float zr = read_imagef(z, sampler, ray).x;

      if (zr > z0)
      {
        fetch += (float)(k + 1); // distance to blocking land
        nhit++;
        break;
      }
    }
  }

  write_imagef(out, g, (nhit > 0) ? fetch / (float)nhit : 0.f);
}

kernel void coastal_fetch_directional(read_only image2d_t  z,
                                      read_only image2d_t  compute_mask,
                                      write_only image2d_t out,
                                      const int            nx,
                                      const int            ny,
                                      const int            ndirections,
                                      const float2         wind_dir,
                                      const float          directional_exp)
{
  // wind_dir: unit vector, prevailing directio

  const int2 g = {get_global_id(0), get_global_id(1)};
  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

  // early skip if requested
  if (read_imagef(compute_mask, sampler, g).x <= 0.f)
  {
    write_imagef(out, g, 0.f);
    return;
  }

  const int    nsteps = (int)(1.414f * max(nx, ny));
  const float2 pos0 = (float2)(g.x + 0.5f, g.y + 0.5f);
  const float  z0 = read_imagef(z, sampler, g).x;

  float fetch_weighted = 0.f;
  float weight_sum = 0.f;

  for (int idir = 0; idir < ndirections; ++idir)
  {
    float  theta = 2.f * M_PI * (float)idir / (float)ndirections;
    float2 dir = (float2)(cos(theta), sin(theta));

    // dot weight: cosine of angle between this ray and prevailing direction
    // clamp to 0 so back-hemisphere directions are ignored
    float w = fmax(0.f, dot(dir, wind_dir));

    // raise to a power to sharpen the directional selectivity
    w = pow(w, directional_exp);

    float ray_fetch = (float)nsteps; // default: open, full fetch

    for (int k = 0; k < nsteps; ++k)
    {
      float2 ray = pos0 + (float)(k + 1) * dir;

      if (!is_inside((int)ray.x, (int)ray.y, nx, ny))
      {
        ray_fetch = (float)nsteps;
        break;
      }

      float zr = read_imagef(z, sampler, ray).x;

      if (zr > z0)
      {
        ray_fetch = (float)(k + 1);
        break;
      }
    }

    fetch_weighted += w * ray_fetch;
    weight_sum += w;
  }

  float result = (weight_sum > 0.f) ? fetch_weighted / weight_sum : 0.f;
  write_imagef(out, g, result);
}
)""
