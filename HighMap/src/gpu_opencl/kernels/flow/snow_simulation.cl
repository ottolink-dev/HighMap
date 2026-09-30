R""(
/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

kernel void snow_simulation(read_only image2d_t  z,
                            read_only image2d_t  s,
                            read_only image2d_t  talus,
                            read_only image2d_t  melting_map,
                            write_only image2d_t s_out,
                            int                  nx,
                            int                  ny,
                            float                dt,
                            float                snow_depth,
                            float                k_snow,
                            float                k_melt,
                            float                k_visc,
                            float                k_depth_ratio,
                            float                k_depth_slope_ratio,
                            float                k_creep,
                            int                  outflow_boundaries)
{
  int i = get_global_id(0);
  int j = get_global_id(1);

  const sampler_t sampler = CLK_NORMALIZED_COORDS_FALSE |
                            CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;

  if (i >= nx || j >= ny) return;

  const float diag = 0.70710678f;
  const float slope_repose = TGET(talus, i, j);
  const float sc = TGET(s, i, j);
  const float zc = TGET(z, i, j);

  //
  const float sc0 = k_depth_ratio * snow_depth;
  const float depth_factor = 1.f +
                             k_depth_slope_ratio * smoothstep(0.f, sc0, sc);
  const float slope_repose_eff = slope_repose * depth_factor;

  float Hc = zc + sc;

  // --- outgoing flux from center

  float sl = (outflow_boundaries && i == 0) ? 0.f : TGET(s, i - 1, j);
  float sr = (outflow_boundaries && i == nx - 1) ? 0.f : TGET(s, i + 1, j);
  float sd = (outflow_boundaries && j == 0) ? 0.f : TGET(s, i, j - 1);
  float su = (outflow_boundaries && j == ny - 1) ? 0.f : TGET(s, i, j + 1);

  float sdl = (outflow_boundaries && (i == 0 || j == 0))
                  ? 0.f
                  : TGET(s, i - 1, j - 1);
  float sdr = (outflow_boundaries && (i == nx - 1 || j == 0))
                  ? 0.f
                  : TGET(s, i + 1, j - 1);
  float sul = (outflow_boundaries && (i == 0 || j == ny - 1))
                  ? 0.f
                  : TGET(s, i - 1, j + 1);
  float sur = (outflow_boundaries && (i == nx - 1 || j == ny - 1))
                  ? 0.f
                  : TGET(s, i + 1, j + 1);

  float Hl = (outflow_boundaries && i == 0) ? zc : (TGET(z, i - 1, j) + sl);
  float Hr = (outflow_boundaries && i == nx - 1) ? zc
                                                 : (TGET(z, i + 1, j) + sr);
  float Hd = (outflow_boundaries && j == 0) ? zc : (TGET(z, i, j - 1) + sd);
  float Hu = (outflow_boundaries && j == ny - 1) ? zc
                                                 : (TGET(z, i, j + 1) + su);

  float Hdl = (outflow_boundaries && (i == 0 || j == 0))
                  ? zc
                  : (TGET(z, i - 1, j - 1) + sdl);
  float Hdr = (outflow_boundaries && (i == nx - 1 || j == 0))
                  ? zc
                  : (TGET(z, i + 1, j - 1) + sdr);
  float Hul = (outflow_boundaries && (i == 0 || j == ny - 1))
                  ? zc
                  : (TGET(z, i - 1, j + 1) + sul);
  float Hur = (outflow_boundaries && (i == nx - 1 || j == ny - 1))
                  ? zc
                  : (TGET(z, i + 1, j + 1) + sur);

  float sx = 0.5f * (Hr - Hl);
  float sy = 0.5f * (Hu - Hd);
  float slope = hypot(sx, sy);

  // Avalanche flux (for slope > repose angle)
  float excess = max(0.f, slope - slope_repose_eff);
  float flux_out_avalanche = k_snow * excess * dt;

  // Flat-terrain creep flux (sub-talus settling/drift)
  float flat_factor = max(0.f, 1.f - (slope / max(slope_repose_eff, 1e-4f)));
  float flux_out_creep = k_creep * flat_factor * slope * dt;

  float flux_out = flux_out_avalanche + flux_out_creep;
  flux_out = fmin(flux_out, sc);

  // depth at which snow starts to lock for avalanche flow
  float mobility = exp(-sc / max(sc0, 1e-4f));
  flux_out = (flux_out_avalanche * mobility) + flux_out_creep;
  flux_out = fmin(flux_out, sc);

  float wl = max(0.f, Hc - Hl);
  float wr = max(0.f, Hc - Hr);
  float wd = max(0.f, Hc - Hd);
  float wu = max(0.f, Hc - Hu);
  float wdl = max(0.f, Hc - Hdl);
  float wdr = max(0.f, Hc - Hdr);
  float wul = max(0.f, Hc - Hul);
  float wur = max(0.f, Hc - Hur);

  float wsum = wl + wr + wd + wu + diag * (wdl + wdr + wul + wur);

  float out = (wsum > 1e-6f) ? flux_out : 0.f;

  // --- incoming flux from neighbors

  float in = 0.f;

  // LEFT neighbor (i-1,j) sends RIGHT to (i,j)
  if (!(outflow_boundaries && i == 0))
  {
    float Hn = Hl;
    float sn = sl;
    float slope_n = max(0.f, Hn - Hc);
    float excess_n = max(0.f, slope_n - slope_repose_eff);
    float flat_n = max(0.f, 1.f - (slope_n / max(slope_repose_eff, 1e-4f)));
    float f_av = k_snow * excess_n * dt * exp(-sn / max(sc0, 1e-4f));
    float f_cr = k_creep * flat_n * slope_n * dt;
    float factor_n = fmin(sn, f_av + f_cr);

    float sum = max(0.f, Hn - Hr) +
                max(0.f, Hn - (TGET(z, i - 1, j - 1) + TGET(s, i - 1, j - 1))) *
                    diag +
                max(0.f, Hn - (TGET(z, i - 1, j + 1) + TGET(s, i - 1, j + 1))) *
                    diag;

    if (sum > 1e-6f) in += factor_n * slope_n / sum;
  }

  // RIGHT neighbor (i+1,j) sends LEFT to (i,j)
  if (!(outflow_boundaries && i == nx - 1))
  {
    float Hn = Hr;
    float sn = sr;
    float slope_n = max(0.f, Hn - Hc);
    float excess_n = max(0.f, slope_n - slope_repose_eff);
    float flat_n = max(0.f, 1.f - (slope_n / max(slope_repose_eff, 1e-4f)));
    float f_av = k_snow * excess_n * dt * exp(-sn / max(sc0, 1e-4f));
    float f_cr = k_creep * flat_n * slope_n * dt;
    float factor_n = fmin(sn, f_av + f_cr);

    float sum = max(0.f, Hn - Hl) +
                max(0.f, Hn - (TGET(z, i + 1, j - 1) + TGET(s, i + 1, j - 1))) *
                    diag +
                max(0.f, Hn - (TGET(z, i + 1, j + 1) + TGET(s, i + 1, j + 1))) *
                    diag;

    if (sum > 1e-6f) in += factor_n * slope_n / sum;
  }

  // DOWN neighbor (i,j-1) sends UP to (i,j)
  if (!(outflow_boundaries && j == 0))
  {
    float Hn = Hd;
    float sn = sd;
    float slope_n = max(0.f, Hn - Hc);
    float excess_n = max(0.f, slope_n - slope_repose_eff);
    float flat_n = max(0.f, 1.f - (slope_n / max(slope_repose_eff, 1e-4f)));
    float f_av = k_snow * excess_n * dt * exp(-sn / max(sc0, 1e-4f));
    float f_cr = k_creep * flat_n * slope_n * dt;
    float factor_n = fmin(sn, f_av + f_cr);

    float sum = max(0.f, Hn - Hu) +
                max(0.f, Hn - (TGET(z, i - 1, j - 1) + TGET(s, i - 1, j - 1))) *
                    diag +
                max(0.f, Hn - (TGET(z, i + 1, j - 1) + TGET(s, i + 1, j - 1))) *
                    diag;

    if (sum > 1e-6f) in += factor_n * slope_n / sum;
  }

  // UP neighbor (i,j+1) sends DOWN to (i,j)
  if (!(outflow_boundaries && j == ny - 1))
  {
    float Hn = Hu;
    float sn = su;
    float slope_n = max(0.f, Hn - Hc);
    float excess_n = max(0.f, slope_n - slope_repose_eff);
    float flat_n = max(0.f, 1.f - (slope_n / max(slope_repose_eff, 1e-4f)));
    float f_av = k_snow * excess_n * dt * exp(-sn / max(sc0, 1e-4f));
    float f_cr = k_creep * flat_n * slope_n * dt;
    float factor_n = fmin(sn, f_av + f_cr);

    float sum = max(0.f, Hn - Hd) +
                max(0.f, Hn - (TGET(z, i - 1, j + 1) + TGET(s, i - 1, j + 1))) *
                    diag +
                max(0.f, Hn - (TGET(z, i + 1, j + 1) + TGET(s, i + 1, j + 1))) *
                    diag;

    if (sum > 1e-6f) in += factor_n * slope_n / sum;
  }

  float sc_new = sc - out + in;

  // --- surface viscosity (smooths total surface H = z + s to fill hollows)

  float H_avg = (Hl + Hr + Hd + Hu + diag * (Hdl + Hdr + Hul + Hur)) /
                (4.f + 4.f * diag);

  float s_target = max(0.f, H_avg - zc);

  sc_new = mix(sc_new, s_target, k_visc);

  // --- melting

  float mt = TGET(melting_map, i, j); // in (0, 1]
  float melt = k_melt * mt * dt;
  melt = fmin(melt, sc_new);
  sc_new -= melt;

  // --- update output

  TSET(s_out, i, j, max(0.f, sc_new));
}
)""
