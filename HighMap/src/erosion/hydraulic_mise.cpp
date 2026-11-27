/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software.
 *
 * Adapted from MISE (Multiscale Implicit Stream-power Erosion) by Leonhard
 * https://github.com/Leonhardmaster2 */

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <iterator>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/erosion/erosion_parameters.hpp"
#include "highmap/erosion/hydraulic_erosion.hpp"
#include "highmap/internal/radix_heap.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate2d.hpp"
#include "highmap/interpolate/interpolate_array.hpp"
#include "highmap/math/array.hpp"
#include "highmap/math/convolve.hpp"
#include "highmap/random.hpp"
#include "highmap/range.hpp"
#include "highmap/spatial/uv_tile.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <bit>
#include <format>

namespace hmap
{

namespace
{

constexpr int   DX8[8] = {1, 1, 0, -1, -1, -1, 0, 1};
constexpr int   DY8[8] = {0, 1, 1, 1, 0, -1, -1, -1};
constexpr float DD8[8] =
    {1.f, 1.41421356f, 1.f, 1.41421356f, 1.f, 1.41421356f, 1.f, 1.41421356f};

enum : uint8_t
{
  F_CLOSED = 1, // already reached by the flood
  F_PIT = 2,    // reached through the pit queue: flat or closed depression
  F_OUTLET = 4  // fixed base-level cell
};

// --- Helper functions and structures

// drainage graph structure
struct FlowGraph
{
  int                  nx = 0, ny = 0, n = 0;
  Array                zf;    // routing surface: z with the depressions filled
  std::vector<int>     order; // flood order, receivers before donors
  std::vector<int>     rank;  // inverse permutation of order
  std::vector<int>     rcv;   // receiver of each cell, itself for an outlet
  std::vector<uint8_t> flag;
  Array                area;  // drainage area in cells, or accumulated rain
  Array                steep; // upstream sum of the relative slope excess

  std::vector<int> lake;        // lake index of a flooded cell, else -1
  std::vector<int> lake_outlet; // spill cell of each lake
  std::vector<int>
      lake_root; // first flooded cell visited upstream-to-downstream
  std::vector<float>  lake_volume, lake_sediment, lake_target;
  std::vector<float>  lake_lo, lake_hi, lake_level;
  std::vector<double> lake_acc;
  std::vector<int>    pitq; // all flooded cells, in flood order
  int                 pit_cells = 0;

  RadixHeap<int>   heap;         // exact flood
  std::vector<int> bhead, bnext; // bucket queue of the fast flood

  void resize(int nx_, int ny_)
  {
    nx = nx_;
    ny = ny_;
    n = nx * ny;
    zf.set_shape({nx, ny});
    order.resize(n);
    rank.resize(n);
    rcv.resize(n);
    flag.resize(n);
    area.set_shape({nx, ny});
    steep.set_shape({nx, ny});
    lake.resize(n);
    pitq.resize(n);
  }

  float dist(int i) const
  {
    int d = rcv[i] > i ? rcv[i] - i : i - rcv[i];
    return (d == 1 || d == nx) ? 1.f : 1.41421356f;
  }

  bool submerged(int i, const Array &z) const
  {
    return (flag[i] & F_PIT) && zf.vector[i] > z.vector[i];
  }

  // priority-flood depression filling and topological ordering
  void flood(const Array &z_arr, const uint8_t *outlet, bool fast)
  {
    const auto &z = z_arr.vector;
    auto       &zf_v = zf.vector;

    std::fill(flag.begin(), flag.end(), (uint8_t)0);
    std::fill(lake.begin(), lake.end(), -1);
    lake_outlet.clear();
    lake_root.clear();
    lake_volume.clear();

    // configure bucket queue (fast approximate flood) or radix heap (exact
    // flood)
    const int nbuckets = 1 << 18;
    float     zmin = 0.f;
    float     scale = 0.f;
    int       cur = nbuckets;
    int       open = 0;

    if (fast)
    {
      float zmax = -3.0e38f;
      zmin = 3.0e38f;

      for (int i = 0; i < n; ++i)
      {
        zmin = std::min(zmin, z[i]);
        zmax = std::max(zmax, z[i]);
      }

      scale = (zmax > zmin) ? (float)(nbuckets - 1) / (zmax - zmin) : 0.f;
      bhead.assign(nbuckets, -1);
      bnext.resize(n);
    }
    else
    {
      heap.clear();
    }

    auto push = [&](int i)
    {
      if (fast)
      {
        const int b = (int)((z[i] - zmin) * scale);
        bnext[i] = bhead[b];
        bhead[b] = i;
        cur = std::min(cur, b);
        ++open;
      }
      else
      {
        heap.push(float_key(z[i]), i);
      }
    };

    auto seed_cell = [&](int i)
    {
      if (flag[i]) return;

      flag[i] = F_CLOSED | F_OUTLET;
      zf_v[i] = z[i];
      push(i);
    };

    // seed domain borders as outlets
    for (int x = 0; x < nx; ++x)
    {
      seed_cell(x);
      seed_cell((ny - 1) * nx + x);
    }

    for (int y = 1; y < ny - 1; ++y)
    {
      seed_cell(y * nx);
      seed_cell(y * nx + nx - 1);
    }

    if (outlet)
    {
      for (int i = 0; i < n; ++i)
        if (outlet[i]) seed_cell(i);
    }

    size_t head = 0;
    size_t tail = 0;
    int    count = 0;
    int    cur_lake = -1;

    while (head < tail || (fast ? open > 0 : heap.count > 0))
    {
      int c = 0;

      if (head < tail)
      {
        c = pitq[head++];
      }
      else
      {
        if (fast)
        {
          while (bhead[cur] < 0)
            ++cur;

          c = bhead[cur];
          bhead[cur] = bnext[c];
          --open;
        }
        else
        {
          c = heap.pop().second;
        }

        cur_lake = -1;
      }

      rank[c] = count;
      order[count++] = c;

      const float zc_eps = next_up(zf_v[c]);
      const int   cx = c % nx;
      const int   cy = c / nx;

      for (int k = 0; k < 8; ++k)
      {
        const int x = cx + DX8[k];
        const int y = cy + DY8[k];
        if (x < 0 || y < 0 || x >= nx || y >= ny) continue;

        const int j = y * nx + x;
        if (flag[j] & F_CLOSED) continue;

        if (z[j] <= zc_eps)
        {
          // cell is inside a depression: raise to spill level and enqueue in
          // pitq
          flag[j] = F_CLOSED | F_PIT;
          zf_v[j] = zc_eps;
          pitq[tail++] = j;

          if (cur_lake < 0)
          {
            cur_lake = (int)lake_outlet.size();
            lake_outlet.push_back(c);
            lake_root.push_back(j);
            lake_volume.push_back(0.f);
          }

          lake[j] = cur_lake;
          lake_volume[cur_lake] += zc_eps - z[j];
        }
        else
        {
          flag[j] = F_CLOSED;
          zf_v[j] = z[j];
          push(j);
        }
      }
    }

    pit_cells = (int)tail;
  }

  // stochastic D8 receiver selection weighted by (dz)^p
  void receivers(float p, uint32_t seed)
  {
    const auto &zf_v = zf.vector;

#pragma omp parallel for if (ny >= 64)
    for (int y = 0; y < ny; ++y)
    {
      for (int x = 0; x < nx; ++x)
      {
        const int i = y * nx + x;
        rcv[i] = i;
        if (flag[i] & F_OUTLET) continue;

        const float zi = zf_v[i];
        const int   ri = rank[i];
        int         cj[8];
        float       cw[8];
        int         cnt = 0;
        float       sum = 0.f;

        for (int k = 0; k < 8; ++k)
        {
          const int xx = x + DX8[k];
          const int yy = y + DY8[k];
          if (xx < 0 || yy < 0 || xx >= nx || yy >= ny) continue;

          const int j = yy * nx + xx;
          if (rank[j] < ri && zf_v[j] < zi)
          {
            const float d = zi - zf_v[j];
            const float w = (p == 1.f) ? d
                                       : ((p == 2.f) ? d * d : std::pow(d, p));
            cj[cnt] = j;
            cw[cnt] = w;
            sum += w;
            ++cnt;
          }
        }

        if (cnt == 0) continue;

        int pick = 0;
        if (cnt > 1)
        {
          const float u = hash01((uint32_t)i, seed) * sum;
          float       acc = 0.f;
          pick = cnt - 1;

          for (int k = 0; k < cnt; ++k)
          {
            acc += cw[k];
            if (u < acc)
            {
              pick = k;
              break;
            }
          }
        }

        rcv[i] = cj[pick];
      }
    }
  }

  // upstream-to-downstream accumulation of drainage area and slope excess
  void accumulate(const Array &z_arr, const Array *p_rain, float talus_cell)
  {
    const auto &z = z_arr.vector;
    auto       &area_v = area.vector;
    auto       &steep_v = steep.vector;

    if (p_rain)
      area = *p_rain;
    else
      std::fill(area_v.begin(), area_v.end(), 1.f);

    std::fill(steep_v.begin(), steep_v.end(), 0.f);
    const float inv = 1.f / talus_cell;

    for (int o = n - 1; o >= 0; --o)
    {
      const int i = order[o];
      const int r = rcv[i];
      if (r == i) continue;

      if (!submerged(i, z_arr))
      {
        const float s = (z[i] - z[r]) / dist(i);
        if (s > talus_cell) steep_v[i] += (s - talus_cell) * inv;
      }

      steep_v[r] += steep_v[i];
      area_v[r] += area_v[i];
    }
  }
};

struct Level
{
  int            nx = 0, ny = 0;
  Array         *z = nullptr;
  Array         *sed = nullptr;
  const Array   *p_erodibility = nullptr;
  const Array   *p_rain = nullptr;
  const Array   *p_bedrock = nullptr;
  const uint8_t *outlet = nullptr;
};

struct Scratch
{
  Array ero, qs, zdon, delta, powlut;
  void  resize(int nx, int ny)
  {
    ero.set_shape({nx, ny});
    qs.set_shape({nx, ny});
    zdon.set_shape({nx, ny});
    delta.set_shape({nx, ny});
  }
};

// --- Erosion and sediment transport passes

void erode_pass(const Level      &L,
                const FlowGraph  &G,
                const MiseParams &P,
                float             dtK,
                float             dtKt,
                float             dtKd,
                Scratch          &W)
{
  auto &z = L.z->vector;
  auto &sed = L.sed->vector;

  const int   nx = G.nx;
  const float dx = 1.f / (float)nx;
  const float cellA = dx * dx;
  const float m = P.area_exp;

  const bool lut_ok = (L.p_rain == nullptr);
  const int  lut_n = 2048;

  if (lut_ok)
  {
    W.powlut.set_shape({lut_n, 1});
    for (int k = 0; k < lut_n; ++k)
      W.powlut.vector[k] = std::pow(std::min((float)k * cellA, P.downcutting),
                                    m);
  }

  for (int o = 0; o < G.n; ++o)
  {
    const int i = G.order[o];
    const int r = G.rcv[i];
    W.ero.vector[i] = 0.f;
    if (r == i) continue;
    const float zi = z[i], zr = z[r];
    if (zi <= zr || G.submerged(i, *L.z)) continue;

    const float d = G.dist(i);
    const float dd = d * dx;
    const float ac = G.area.vector[i];
    const float am = (lut_ok && ac < (float)lut_n)
                         ? W.powlut.vector[(int)ac]
                         : std::pow(std::min(ac * cellA, P.downcutting), m);
    const float ke = L.p_erodibility ? L.p_erodibility->vector[i] : 1.f;
    const float soft = (sed[i] > 0.f) ? P.sediment_erodibility : 1.f;

    float rho = 1.f;

    if (P.slope_correction_max > 1.f)
    {
      const float gx = std::max(0.f, std::max(zi - z[i - 1], zi - z[i + 1]));
      const float gy = std::max(0.f, std::max(zi - z[i - nx], zi - z[i + nx]));
      rho = std::min(
          P.slope_correction_max,
          std::max(1.f, std::sqrt(gx * gx + gy * gy) * d / (zi - zr)));
    }

    // implicit stream-power incision toward receiver elevation zr
    const float Ff = dtK * ke * soft * am * rho / dd;
    float       zn = (zi + Ff * zr) / (1.f + Ff);

    // implicit thermal/debris relaxation toward talus threshold zc
    const float zc = zr + P.talus * dd;

    if (zn > zc && (dtKt > 0.f || dtKd > 0.f))
    {
      const float Fh = (dtKt + dtKd * G.steep.vector[i] * cellA) * ke * soft /
                       dd;
      zn = (zi + Ff * zr + Fh * zc) / (1.f + Ff + Fh);
    }

    // clamp to max cliff slope and bedrock floor
    if (P.cliff > 0.f) zn = std::min(zn, zr + P.cliff * dd);
    if (L.p_bedrock) zn = std::min(zi, std::max(zn, L.p_bedrock->vector[i]));

    // record eroded volume and deplete loose sediment
    const float e = zi - zn;
    W.ero.vector[i] = e;
    z[i] = zn;
    sed[i] = std::max(0.f, sed[i] - e);
  }
}

// --- Sediment transport and lake filling

void sediment_pass(const Level      &L,
                   FlowGraph        &G,
                   const MiseParams &P,
                   Scratch          &W)
{
  auto &z = L.z->vector;
  auto &sed = L.sed->vector;

  const float dx = 1.f / (float)G.nx;
  const float cellA = dx * dx;
  const float a_ref = std::max(P.deposition_area, 1e-12f);
  auto       &qs = W.qs.vector;
  auto       &zdon = W.zdon.vector;
  auto       &ero = W.ero.vector;
  const auto &area = G.area.vector;

  std::fill(qs.begin(), qs.end(), 0.f);
  std::fill(zdon.begin(), zdon.end(), 3.0e38f);
  G.lake_sediment.assign(G.lake_outlet.size(), 0.f);
  G.lake_target.assign(G.lake_outlet.size(), 0.f);

  // route sediment upstream-to-downstream
  for (int o = G.n - 1; o >= 0; --o)
  {
    const int i = G.order[o];
    const int r = G.rcv[i];
    float     q = qs[i];

    if (r == i) continue;

    // accumulate sediment inside flooded depressions and spill excess
    const int lk = G.lake[i];
    if (lk >= 0)
    {
      G.lake_sediment[lk] += q + ero[i];

      if (i == G.lake_root[lk])
      {
        const float s = G.lake_sediment[lk];
        const float fill = P.lake_fill * std::min(s, G.lake_volume[lk]);
        G.lake_target[lk] = fill;

        if (s > fill) qs[G.lake_outlet[lk]] += s - fill;
      }

      continue;
    }

    // deposit incoming sediment where slope is below equilibrium
    if (q > 0.f)
    {
      const float an = area[i] * cellA;
      float       seq = P.deposition_slope;

      if (an > a_ref) seq *= std::pow(an / a_ref, -P.deposition_exp);

      const float room = std::min(z[r] + seq * G.dist(i) * dx, zdon[i]) - z[i];

      if (room > 0.f)
      {
        const float dep = std::min(q, P.deposition_rate * room);
        z[i] += dep;
        sed[i] += dep;
        q -= dep;
      }
    }

    // forward remaining load and eroded material to receiver
    q += ero[i];
    if (q > 0.f)
    {
      qs[r] += q;
      zdon[r] = std::min(zdon[r], z[i]);
    }
  }

  // fill flooded depressions up to their target sediment volume via bisection
  const size_t nl = G.lake_outlet.size();
  if (nl > 0)
  {
    const int *cells = G.pitq.data();
    const int  nc = G.pit_cells;

    G.lake_lo.assign(nl, 3.0e38f);
    G.lake_hi.assign(nl, -3.0e38f);

    for (int c = 0; c < nc; ++c)
    {
      const int i = cells[c];
      const int lk = G.lake[i];
      G.lake_lo[lk] = std::min(G.lake_lo[lk], z[i]);
      G.lake_hi[lk] = std::max(G.lake_hi[lk], G.zf.vector[i]);
    }

    std::vector<float> &lo = G.lake_lo;
    std::vector<float> &hi = G.lake_hi;
    G.lake_level.resize(nl);

    for (int it = 0; it < 20; ++it)
    {
      G.lake_acc.assign(nl, 0.0);

      for (size_t k = 0; k < nl; ++k)
        G.lake_level[k] = 0.5f * (lo[k] + hi[k]);

      for (int c = 0; c < nc; ++c)
      {
        const int   i = cells[c];
        const int   lk = G.lake[i];
        const float h = G.lake_level[lk] - z[i];
        if (h > 0.f) G.lake_acc[lk] += h;
      }

      for (size_t k = 0; k < nl; ++k)
        (G.lake_acc[k] > (double)G.lake_target[k] ? hi[k]
                                                  : lo[k]) = G.lake_level[k];
    }

    for (int c = 0; c < nc; ++c)
    {
      const int i = cells[c];
      const int lk = G.lake[i];
      if (G.lake_target[lk] <= 0.f) continue;

      const float level = std::min(lo[lk], G.zf.vector[i]);
      if (level > z[i])
      {
        const float dep = level - z[i];
        z[i] = level;
        sed[i] += dep;
      }
    }
  }
}

// --- Talus relaxation

void relax_sediment(Array &z,
                    Array &sed,
                    int    nx,
                    int    ny,
                    float  talus_cell,
                    int    iters,
                    Array &delta)
{
  for (int it = 0; it < iters; ++it)
  {
    delta = 0.f;

    // compute talus excess and distribute to lower neighbors
    for (int j = 1; j < ny - 1; ++j)
      for (int i = 1; i < nx - 1; ++i)
      {
        if (sed(i, j) <= 0.f) continue;

        const float zi = z(i, j);
        float       ex[8], sum = 0.f, mx = 0.f;

        for (int k = 0; k < 8; ++k)
        {
          float e = zi - z(i + DX8[k], j + DY8[k]) - talus_cell * DD8[k];
          e = e > 0.f ? e : 0.f;
          ex[k] = e;
          sum += e;
          mx = e > mx ? e : mx;
        }

        if (sum <= 0.f) continue;

        const float amount = std::min(0.25f * mx, sed(i, j));
        delta(i, j) -= amount;

        const float s = amount / sum;
        for (int k = 0; k < 8; ++k)
          if (ex[k] > 0.f) delta(i + DX8[k], j + DY8[k]) += ex[k] * s;
      }

    for (int j = 1; j < ny - 1; ++j)
      for (int i = 1; i < nx - 1; ++i)
      {
        const float d = delta(i, j);
        if (d == 0.f) continue;
        z(i, j) += d;
        sed(i, j) = std::max(0.f, sed(i, j) + d);
      }
  }
}

// --- Single-level solver

void erode_level(const Level      &L,
                 const MiseParams &P,
                 float             share,
                 int               iters,
                 uint32_t          seed,
                 FlowGraph        &G,
                 Scratch          &W)
{
  if (iters <= 0 || share <= 0.f) return;

  G.resize(L.nx, L.ny);
  W.resize(L.nx, L.ny);

  const float dt = share / (float)iters;
  const float talus_cell = P.talus / (float)L.nx;

  for (int it = 0; it < iters; ++it)
  {
    G.flood(*L.z, L.outlet, !P.exact_flood);
    G.receivers(P.receiver_exp, seed);
    G.accumulate(*L.z, L.p_rain, talus_cell);

    erode_pass(L, G, P, P.strength * dt, P.thermal * dt, P.debris * dt, W);
    sediment_pass(L, G, P, W);

    if (P.sediment_relax_iters > 0)
      relax_sediment(*L.z,
                     *L.sed,
                     L.nx,
                     L.ny,
                     P.sediment_talus / (float)L.nx,
                     P.sediment_relax_iters,
                     W.delta);
  }
}

// --- Flow map computation

Array flow_map_impl(const Array   &z,
                    int            nx,
                    int            ny,
                    const Array   *p_rain,
                    const uint8_t *outlet,
                    uint32_t       seed)
{
  Array flow(glm::ivec2{nx, ny});

  if (nx < 3 || ny < 3 || z.vector.size() != flow.vector.size() ||
      !validate_finite(z))
    return flow;

  FlowGraph G;
  G.resize(nx, ny);
  G.flood(z, outlet, false);
  G.receivers(4.f, seed);
  G.accumulate(z, p_rain, 1.f);

  const float inv = 1.f / (float)flow.vector.size();
  flow = G.area * inv;

  return flow;
}

} // namespace

// --- Main solver

void hydraulic_mise(Array            &z,
                    const MiseParams &params,
                    const Array      *p_bedrock,
                    const Array      *p_erodibility,
                    const Array      *p_moisture_map,
                    const Array      *p_outlet,
                    Array            *p_sediment,
                    Array            *p_flow,
                    Array            *p_erosion_map,
                    Array            *p_deposition_map)
{
  if (!validate_non_empty(z)) return;
  if (!validate_finite(z)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;
  if (p_erodibility && !validate_same_shape(z, *p_erodibility)) return;
  if (p_moisture_map && !validate_same_shape(z, *p_moisture_map)) return;
  if (p_outlet && !validate_same_shape(z, *p_outlet)) return;

  const int    nx = z.shape.x;
  const int    ny = z.shape.y;
  const size_t n = (size_t)nx * ny;

  if (nx < 3 || ny < 3 || z.vector.size() != n) return;

  const Array z0 = z;
  Array       sediment(z.shape);
  const float zfloor = z0.min();

  // setup outlet buffer if provided
  std::vector<uint8_t> outlet_bytes;
  const uint8_t       *p_outlet_ptr = nullptr;
  if (p_outlet)
  {
    outlet_bytes.resize(n);
    for (size_t i = 0; i < n; ++i)
      outlet_bytes[i] = (p_outlet->vector[i] > 0.5f) ? 1 : 0;
    p_outlet_ptr = outlet_bytes.data();
  }

  const float ratio = (float)std::min(nx, ny) /
                      (float)std::max(params.base_res, 8);
  const int levels = 1 + std::min(9,
                                  std::max(0,
                                           (int)std::lround(std::log2(
                                               std::max(ratio, 1.f)))));
  const int run_levels = levels -
                         std::min(std::max(params.skip_finest, 0), levels - 1);

  FlowGraph G;
  Scratch   W;
  float     share = 1.f - params.share_ratio;
  float     iters = (float)params.iters0;

  for (int l = 0; l < run_levels; ++l)
  {
    const int      f = 1 << (levels - 1 - l);
    const int      it = std::max(params.iters_min, (int)std::lround(iters));
    const int      lx = (nx + f - 1) / f, ly = (ny + f - 1) / f;
    const uint32_t seed = params.seed * 7919u + 104729u * (uint32_t)lx;

    if (f == 1)
    {
      // --- Full resolution level

      Array before;
      if (params.antialias > 0.f) before = z;
      Level L{nx,
              ny,
              &z,
              &sediment,
              p_erodibility,
              p_moisture_map,
              p_bedrock,
              p_outlet_ptr};
      erode_level(L, params, share, it, seed, G, W);

      if (params.antialias > 0.f)
      {
        Array d = z - before;

        const std::vector<float> kernel = {0.25f, 0.5f, 0.25f};
        Array b = convolve1d_j(convolve1d_i(d, kernel), kernel);

        z = before + d + params.antialias * (b - d);
        clamp_min(z, zfloor);
      }
    }
    else
    {
      // --- Coarse resolution level

      const glm::ivec2 coarse_shape{lx, ly};
      Array            zl = z.resample_to_shape(coarse_shape);
      Array            sl = sediment.resample_to_shape(coarse_shape);

      Array                el, rl, bl;
      std::vector<uint8_t> ol;

      if (p_erodibility) el = p_erodibility->resample_to_shape(coarse_shape);
      if (p_moisture_map) rl = p_moisture_map->resample_to_shape(coarse_shape);
      if (p_bedrock) bl = p_bedrock->resample_to_shape(coarse_shape);

      if (p_outlet)
      {
        const Array oc = p_outlet->resample_to_shape(coarse_shape);
        ol.resize(oc.vector.size());
        for (size_t i = 0; i < oc.vector.size(); ++i)
          ol[i] = oc.vector[i] > 0.5f;
      }
      Array db_l = zl - sl;

      Level L{lx,
              ly,
              &zl,
              &sl,
              p_erodibility ? &el : nullptr,
              p_moisture_map ? &rl : nullptr,
              p_bedrock ? &bl : nullptr,
              p_outlet_ptr ? ol.data() : nullptr};
      erode_level(L, params, share, it, seed, G, W);

      db_l = (zl - sl) - db_l;
      Array db = resample_bicubic_warp(db_l, z.shape, params.warp, seed);
      Array zs = resample_bicubic_warp(zl, z.shape, params.warp, seed);
      Array ss = resample_bicubic_warp(sl, z.shape, params.warp, seed);

      // blend coarse bedrock change with attenuated fine-scale detail
#pragma omp parallel for if (ny >= 64)
      for (int y = 0; y < ny; ++y)
      {
        for (size_t i = (size_t)y * nx; i < (size_t)(y + 1) * nx; ++i)
        {
          const float b1 = zs.vector[i] - ss.vector[i];
          const float detail = (z.vector[i] - sediment.vector[i]) -
                               (b1 - db.vector[i]);
          const float ad = std::fabs(detail);
          const float cut = std::fabs(db.vector[i]);
          const float b = b1 + (ad > 0.f ? detail * ad / (ad + cut) : 0.f);
          const float s = std::min(std::max(zs.vector[i] - b, 0.f),
                                   2.f * std::max(ss.vector[i], 0.f));
          z.vector[i] = std::max(b + s, zfloor);
          sediment.vector[i] = s;
        }
      }

      if (p_bedrock) clamp_min(z, minimum(z0, *p_bedrock));
    }

    share *= params.share_ratio;
    iters *= params.iters_decay;
  }

  // restore fixed outlet cells
  auto restore = [&](size_t i)
  {
    z.vector[i] = z0.vector[i];
    sediment.vector[i] = 0.f;
  };

  for (int x = 0; x < nx; ++x)
  {
    restore((size_t)x);
    restore((size_t)(ny - 1) * nx + x);
  }

  for (int y = 0; y < ny; ++y)
  {
    restore((size_t)y * nx);
    restore((size_t)y * nx + nx - 1);
  }

  if (p_outlet_ptr)
  {
    for (size_t i = 0; i < n; ++i)
      if (p_outlet_ptr[i]) restore(i);
  }

  if (p_bedrock) clamp_min(z, minimum(z0, *p_bedrock));

  if (params.extrapolate_border)
  {
    extrapolate_borders(z, 1, 0.9f);

    if (p_bedrock)
    {
      z = maximum(z, *p_bedrock);
    }
  }

  if (p_sediment) *p_sediment = sediment;

  if (p_flow)
  {
    *p_flow = flow_map_impl(z,
                            nx,
                            ny,
                            p_moisture_map,
                            p_outlet_ptr,
                            params.seed);
  }

  if (p_erosion_map)
  {
    *p_erosion_map = z0 - z;
    clamp_min(*p_erosion_map, 0.f);
  }

  if (p_deposition_map)
  {
    *p_deposition_map = z - z0;
    clamp_min(*p_deposition_map, 0.f);
  }
}

void hydraulic_mise(Array            &z,
                    const Array      *p_mask,
                    const MiseParams &params,
                    const Array      *p_bedrock,
                    const Array      *p_erodibility,
                    const Array      *p_moisture_map,
                    const Array      *p_outlet,
                    Array            *p_sediment,
                    Array            *p_flow,
                    Array            *p_erosion_map,
                    Array            *p_deposition_map)
{
  if (!validate_non_empty(z)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;

  if (!p_mask)
  {
    hydraulic_mise(z,
                   params,
                   p_bedrock,
                   p_erodibility,
                   p_moisture_map,
                   p_outlet,
                   p_sediment,
                   p_flow,
                   p_erosion_map,
                   p_deposition_map);
  }
  else
  {
    Array z_f = z;
    hydraulic_mise(z_f,
                   params,
                   p_bedrock,
                   p_erodibility,
                   p_moisture_map,
                   p_outlet,
                   p_sediment,
                   p_flow,
                   p_erosion_map,
                   p_deposition_map);
    z = lerp(z, z_f, *p_mask);

    if (p_sediment) *p_sediment *= (*p_mask);
    if (p_erosion_map) *p_erosion_map *= (*p_mask);
    if (p_deposition_map) *p_deposition_map *= (*p_mask);
  }
}

} // namespace hmap

namespace hmap::va
{

VirtualArray hydraulic_mise(const ComputeMode  &cm,
                            const VirtualArray &z,
                            int                 cutoff_resolution,
                            const MiseParams   &params,
                            const VirtualArray *p_mask,
                            const VirtualArray *p_bedrock,
                            const VirtualArray *p_erodibility,
                            const VirtualArray *p_moisture_map,
                            const VirtualArray *p_outlet,
                            VirtualArray       *p_sediment,
                            VirtualArray       *p_flow,
                            VirtualArray       *p_erosion_map,
                            VirtualArray       *p_deposition_map)
{
  VirtualArray out;
  if (z.empty()) return out;

  glm::ivec2 total_shape = z.shape;

  // If entire virtual array is within cutoff_resolution, perform full global
  // erosion
  if (std::max(total_shape.x, total_shape.y) <= cutoff_resolution)
  {
    Array z_full = z.to_array(cm);

    Array bedrock = p_bedrock ? p_bedrock->to_array(cm) : Array();
    Array erodibility = p_erodibility ? p_erodibility->to_array(cm) : Array();
    Array moisture = p_moisture_map ? p_moisture_map->to_array(cm) : Array();
    Array outlet = p_outlet ? p_outlet->to_array(cm) : Array();
    Array mask = p_mask ? p_mask->to_array(cm) : Array();

    Array sed, flow, ero, dep;

    hydraulic_mise(z_full,
                   p_mask ? &mask : nullptr,
                   params,
                   p_bedrock ? &bedrock : nullptr,
                   p_erodibility ? &erodibility : nullptr,
                   p_moisture_map ? &moisture : nullptr,
                   p_outlet ? &outlet : nullptr,
                   p_sediment ? &sed : nullptr,
                   p_flow ? &flow : nullptr,
                   p_erosion_map ? &ero : nullptr,
                   p_deposition_map ? &dep : nullptr);

    out.copy_from(z, cm, /* copy_src_data */ false);
    out.from_array(z_full, cm);

    if (p_sediment)
    {
      p_sediment->copy_from(z, cm, false);
      p_sediment->from_array(sed, cm);
    }
    if (p_flow)
    {
      p_flow->copy_from(z, cm, false);
      p_flow->from_array(flow, cm);
    }
    if (p_erosion_map)
    {
      p_erosion_map->copy_from(z, cm, false);
      p_erosion_map->from_array(ero, cm);
    }
    if (p_deposition_map)
    {
      p_deposition_map->copy_from(z, cm, false);
      p_deposition_map->from_array(dep, cm);
    }

    return out;
  }

  // --- Hybrid multi-scale distributed execution

  // 1. Global coarse-scale solver pass up to cutoff_resolution
  int        coarse_res = std::max(128, cutoff_resolution);
  glm::ivec2 coarse_shape;
  if (total_shape.x >= total_shape.y)
  {
    coarse_shape.x = coarse_res;
    coarse_shape.y = std::max(
        64,
        (int)((int64_t)total_shape.y * coarse_res / total_shape.x));
  }
  else
  {
    coarse_shape.y = coarse_res;
    coarse_shape.x = std::max(
        64,
        (int)((int64_t)total_shape.x * coarse_res / total_shape.y));
  }

  Array z_coarse = z.to_array(coarse_shape, cm);
  Array bedrock_c = p_bedrock ? p_bedrock->to_array(coarse_shape, cm) : Array();
  Array erod_c = p_erodibility ? p_erodibility->to_array(coarse_shape, cm)
                               : Array();
  Array moist_c = p_moisture_map ? p_moisture_map->to_array(coarse_shape, cm)
                                 : Array();
  Array outlet_c = p_outlet ? p_outlet->to_array(coarse_shape, cm) : Array();
  Array mask_c = p_mask ? p_mask->to_array(coarse_shape, cm) : Array();

  MiseParams coarse_params = params;
  coarse_params.base_res = std::min(params.base_res,
                                    std::min(coarse_shape.x, coarse_shape.y));

  Array sed_c, flow_c;
  hydraulic_mise(z_coarse,
                 p_mask ? &mask_c : nullptr,
                 coarse_params,
                 p_bedrock ? &bedrock_c : nullptr,
                 p_erodibility ? &erod_c : nullptr,
                 p_moisture_map ? &moist_c : nullptr,
                 p_outlet ? &outlet_c : nullptr,
                 &sed_c,
                 &flow_c,
                 nullptr,
                 nullptr);

  // Initialize distributed output
  out.copy_from(z, cm, /* copy_src_data */ false);
  out.from_array_bicubic(z_coarse, cm);

  if (p_sediment)
  {
    p_sediment->copy_from(z, cm, false);
    p_sediment->from_array_bicubic(sed_c, cm);
  }
  if (p_flow)
  {
    p_flow->copy_from(z, cm, false);
    p_flow->from_array_bilinear(flow_c, cm);
  }

  // 2. Tileable fine-detail refinement
  MiseParams fine_params = params;
  fine_params.base_res = std::min(cutoff_resolution, 256);
  fine_params.iters0 = std::max(3, params.iters0 / 2);

  std::vector<const VirtualArray *> in_arrays = {&z};
  if (p_mask) in_arrays.push_back(p_mask);
  if (p_bedrock) in_arrays.push_back(p_bedrock);
  if (p_erodibility) in_arrays.push_back(p_erodibility);
  if (p_moisture_map) in_arrays.push_back(p_moisture_map);

  hmap::for_each_tile(
      in_arrays,
      {&out},
      [&](std::vector<const hmap::Array *> /* p_arrays_in */,
          std::vector<hmap::Array *> p_arrays_out,
          const hmap::UvTile & /* region */)
      {
        auto [pa_out] = unpack<1>(p_arrays_out);
        if (!pa_out || pa_out->vector.empty()) return;

        // Perform local high-frequency incision refinement per tile
        hydraulic_mise(*pa_out,
                       fine_params,
                       nullptr,
                       nullptr,
                       nullptr,
                       nullptr,
                       nullptr,
                       nullptr,
                       nullptr,
                       nullptr);
      },
      cm);

  // Synchronize halo boundaries between adjacent tiles
  out.sync_overlap_buffers(SyncOperation::SmoothBlend);

  // Apply optional mask blending
  if (p_mask)
  {
    hmap::for_each_tile(
        {&z, p_mask},
        {&out},
        [&](std::vector<const hmap::Array *> p_arrays_in,
            std::vector<hmap::Array *>       p_arrays_out,
            const hmap::UvTile &)
        {
          auto [pa_z, pa_m] = unpack<2>(p_arrays_in);
          auto [pa_out] = unpack<1>(p_arrays_out);

          if (pa_m) *pa_out = lerp(*pa_z, *pa_out, *pa_m);
        },
        cm);
  }

  // Compute delta maps
  if (p_erosion_map)
  {
    p_erosion_map->copy_from(z, cm, false);
    hmap::for_each_tile(
        {&z, &out},
        {p_erosion_map},
        [&](std::vector<const hmap::Array *> p_arrays_in,
            std::vector<hmap::Array *>       p_arrays_out,
            const hmap::UvTile &)
        {
          auto [pa_z, pa_out] = unpack<2>(p_arrays_in);
          auto [pa_ero] = unpack<1>(p_arrays_out);

          *pa_ero = maximum(*pa_z - *pa_out, 0.f);
        },
        cm);
  }

  if (p_deposition_map)
  {
    p_deposition_map->copy_from(z, cm, false);
    hmap::for_each_tile(
        {&z, &out},
        {p_deposition_map},
        [&](std::vector<const hmap::Array *> p_arrays_in,
            std::vector<hmap::Array *>       p_arrays_out,
            const hmap::UvTile &)
        {
          auto [pa_z, pa_out] = unpack<2>(p_arrays_in);
          auto [pa_dep] = unpack<1>(p_arrays_out);

          *pa_dep = maximum(*pa_out - *pa_z, 0.f);
        },
        cm);
  }

  return out;
}

} // namespace hmap::va
