R""(
/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

inline void musgrave_cell_outflow(const float w_curr,
                                  const float h_curr,
                                  const float s_curr,
                                  const float h_nbrs[8],
                                  const float c_capacity,
                                  const float c_erosion,
                                  const float c_deposition,
                                  float      *out_dw_tot,
                                  float       out_dw[8],
                                  float       out_ds[8],
                                  float      *out_s_kept,
                                  float      *out_dz)
{
  // using d8_inv_dist

  float dh_tot = 0.f;
  float dh[8];

  for (int k = 0; k < 8; ++k)
  {
    float diff = (h_curr - h_nbrs[k]) * d8_inv_dist[k];
    dh[k] = max(0.f, diff);
    dh_tot += dh[k];
  }

  if (dh_tot > 0.f)
  {
    float dw_tot = min(w_curr, dh_tot);
    *out_dw_tot = dw_tot;

    float cs_tot = c_capacity * dw_tot;

    if (s_curr > cs_tot)
    {
      float ds_dep = c_deposition * (s_curr - cs_tot);
      *out_dz = ds_dep;
      *out_s_kept = (1.f - c_deposition) * (s_curr - cs_tot);

      for (int k = 0; k < 8; ++k)
      {
        float frac = dh[k] / dh_tot;
        out_dw[k] = dw_tot * frac;
        out_ds[k] = c_capacity * out_dw[k];
      }
    }
    else
    {
      float ds_ero = c_erosion * (cs_tot - s_curr);
      *out_dz = -ds_ero;
      *out_s_kept = 0.f;

      float s_tot_moving = s_curr + ds_ero;

      for (int k = 0; k < 8; ++k)
      {
        float frac = dh[k] / dh_tot;
        out_dw[k] = dw_tot * frac;
        out_ds[k] = s_tot_moving * frac;
      }
    }
  }
  else
  {
    *out_dw_tot = 0.f;
    for (int k = 0; k < 8; ++k)
    {
      out_dw[k] = 0.f;
      out_ds[k] = 0.f;
    }

    float ds_dep = c_deposition * s_curr;
    *out_dz = ds_dep;
    *out_s_kept = (1.f - c_deposition) * s_curr;
  }
}

void kernel hydraulic_musgrave(read_only image2d_t  z,
                               read_only image2d_t  w,
                               read_only image2d_t  s,
                               read_only image2d_t  moisture_map,
                               write_only image2d_t z_new,
                               write_only image2d_t w_new,
                               write_only image2d_t s_new,
                               const int            nx,
                               const int            ny,
                               const float          c_capacity,
                               const float          c_erosion,
                               const float          c_deposition,
                               const float          water_level,
                               const float          evap_rate)
{
  const int2 g = {get_global_id(0), get_global_id(1)};

  if (g.x >= nx || g.y >= ny) return;

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  const int i = g.x;
  const int j = g.y;

  if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
  {
    TSET(z_new, i, j, TGET(z, i, j));
    TSET(w_new, i, j, TGET(w, i, j));
    TSET(s_new, i, j, TGET(s, i, j));
    return;
  }

  // using d8_di, d8_dj, d8_rev

  // effective current water and surface height with rain / evaporation
  float w_center = (1.f - evap_rate) * TGET(w, i, j) +
                   evap_rate * TGET(moisture_map, i, j) * water_level;
  float z_center = TGET(z, i, j);
  float s_center = TGET(s, i, j);
  float h_center = z_center + w_center;

  // collect heights around center
  float h_nbrs_center[8];
  for (int k = 0; k < 8; ++k)
  {
    int   ni = i + d8_di[k];
    int   nj = j + d8_dj[k];
    float w_nbr = (1.f - evap_rate) * TGET(w, ni, nj) +
                  evap_rate * TGET(moisture_map, ni, nj) * water_level;
    h_nbrs_center[k] = TGET(z, ni, nj) + w_nbr;
  }

  float dw_out_tot = 0.f;
  float dw_out[8];
  float ds_out[8];
  float s_kept = 0.f;
  float dz_self = 0.f;

  musgrave_cell_outflow(w_center,
                        h_center,
                        s_center,
                        h_nbrs_center,
                        c_capacity,
                        c_erosion,
                        c_deposition,
                        &dw_out_tot,
                        dw_out,
                        ds_out,
                        &s_kept,
                        &dz_self);

  // gather incoming flow from neighbors
  float w_in_tot = 0.f;
  float s_in_tot = 0.f;

  for (int k = 0; k < 8; ++k)
  {
    int ni = i + d8_di[k];
    int nj = j + d8_dj[k];

    if (ni > 0 && ni < nx - 1 && nj > 0 && nj < ny - 1)
    {
      float w_nbr = (1.f - evap_rate) * TGET(w, ni, nj) +
                    evap_rate * TGET(moisture_map, ni, nj) * water_level;
      float z_nbr = TGET(z, ni, nj);
      float s_nbr = TGET(s, ni, nj);
      float h_nbr = z_nbr + w_nbr;

      float h_nbrs_of_nbr[8];
      for (int m = 0; m < 8; ++m)
      {
        int   nni = ni + d8_di[m];
        int   nnj = nj + d8_dj[m];
        float w_nnbr = (1.f - evap_rate) * TGET(w, nni, nnj) +
                       evap_rate * TGET(moisture_map, nni, nnj) * water_level;
        h_nbrs_of_nbr[m] = TGET(z, nni, nnj) + w_nnbr;
      }

      float nbr_dw_tot = 0.f;
      float nbr_dw_out[8];
      float nbr_ds_out[8];
      float nbr_s_kept = 0.f;
      float nbr_dz = 0.f;

      musgrave_cell_outflow(w_nbr,
                            h_nbr,
                            s_nbr,
                            h_nbrs_of_nbr,
                            c_capacity,
                            c_erosion,
                            c_deposition,
                            &nbr_dw_tot,
                            nbr_dw_out,
                            nbr_ds_out,
                            &nbr_s_kept,
                            &nbr_dz);

      int k_to_me = d8_rev[k];
      w_in_tot += nbr_dw_out[k_to_me];
      s_in_tot += nbr_ds_out[k_to_me];
    }
  }

  TSET(z_new, i, j, z_center + dz_self);
  TSET(w_new, i, j, max(0.f, w_center - dw_out_tot + w_in_tot));
  TSET(s_new, i, j, max(0.f, s_kept + s_in_tot));
}
)""
