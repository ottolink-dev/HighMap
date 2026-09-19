R""(
/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
__constant const int   thermal_olsen_di[8] = {-1, 0, 0, 1, -1, -1, 1, 1};
__constant const int   thermal_olsen_dj[8] = {0, 1, -1, 0, -1, 1, -1, 1};
__constant const float thermal_olsen_c[8] = {1.f,
                                             1.f,
                                             1.f,
                                             1.f,
                                             1.414213562f,
                                             1.414213562f,
                                             1.414213562f,
                                             1.414213562f};
__constant const int   thermal_olsen_rev[8] = {3, 2, 1, 0, 7, 6, 5, 4};

void kernel thermal_olsen(global const float *z_in,
                          global float       *z_out,
                          const global float *talus,
                          const int           nx,
                          const int           ny)
{
  int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;
  if (apply_boundaries_buffer_io(z_in, z_out, g.x, g.y, nx, ny, 3)) return;

  // --- thermal erosion

  int   idx = linear_index(g.x, g.y, nx);
  float delta = 0.f;
  float z_val = z_in[idx];
  float talus_val = talus[idx];

  // 1. Loss from current cell (g) to lower neighbors:
  float dsum = 0.f;
  float dmax = 0.f;
  float dz[8];

#pragma unroll
  for (int p = 0; p < 8; ++p)
  {
    int2  gp = (int2)(g.x + thermal_olsen_di[p], g.y + thermal_olsen_dj[p]);
    float diff = z_val - z_in[linear_index(gp.x, gp.y, nx)];
    dz[p] = diff;

    if (diff > talus_val * thermal_olsen_c[p])
    {
      dsum += diff;
      dmax = max(dmax, diff);
    }
  }

  if (dmax > 0.f && dsum > 0.f)
  {
#pragma unroll
    for (int p = 0; p < 8; ++p)
    {
      if (dz[p] > talus_val * thermal_olsen_c[p])
      {
        float amount = 0.25f * (dmax - talus_val * thermal_olsen_c[p]) * dz[p] /
                       dsum;
        delta -= amount;
      }
    }
  }

  // 2. Gain to current cell (g) from higher neighbors (gk):
  uint shift = hash21u((uint)g.x, (uint)g.y) & 7;

  for (int k = 0; k < 8; ++k)
  {
    int  s = (k + shift) & 7;
    int2 gk = (int2)(g.x + thermal_olsen_di[s], g.y + thermal_olsen_dj[s]);
    int  gk_idx = linear_index(gk.x, gk.y, nx);

    float talus_k = talus[gk_idx];
    float z_k = z_in[gk_idx];

    float dsum_k = 0.f;
    float dmax_k = 0.f;
    float dz_k[8];

#pragma unroll
    for (int p = 0; p < 8; ++p)
    {
      int2  gp = (int2)(gk.x + thermal_olsen_di[p], gk.y + thermal_olsen_dj[p]);
      float diff = z_k - z_in[linear_index(gp.x, gp.y, nx)];
      dz_k[p] = diff;

      if (diff > talus_k * thermal_olsen_c[p])
      {
        dsum_k += diff;
        dmax_k = max(dmax_k, diff);
      }
    }

    if (dmax_k > 0.f && dsum_k > 0.f)
    {
      int   rev_p = thermal_olsen_rev[s];
      float diff_to_c = dz_k[rev_p];
      if (diff_to_c > talus_k * thermal_olsen_c[rev_p])
      {
        float amount = 0.25f * (dmax_k - talus_k * thermal_olsen_c[rev_p]) *
                       diff_to_c / dsum_k;
        delta += amount;
      }
    }
  }

  z_out[idx] = z_val + delta;
}
)""
