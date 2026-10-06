/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software.
 *
 * Adapted from MISE (Multiscale Implicit Stream-power Erosion) by Leonhard
 * https://github.com/Leonhardmaster2 */

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <thread>
#include <utility>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/blending.hpp"
#include "highmap/erosion.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/interpolate/interpolate2d.hpp"
#include "highmap/math/array.hpp"
#include "highmap/random.hpp"
#include "highmap/range.hpp"

#include <bit>

namespace hmap
{

namespace
{

using Field = std::vector<float>;
using Clock = std::chrono::steady_clock;

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

// test all floats for non-nan and non-infinite values
bool all_finite(const float *v, size_t n)
{
  for (size_t i = 0; i < n; ++i)
  {
    uint32_t u = 0;
    std::memcpy(&u, v + i, 4);
    if ((u & 0x7f800000u) == 0x7f800000u) return false;
  }
  return true;
}

// area-weighted box downsampling
template <class F> void parallel_rows(int rows, int threads, F &&fn)
{
  if (threads <= 1 || rows < 64)
  {
    fn(0, rows);
    return;
  }
  std::vector<std::thread> pool;
  pool.reserve(threads);
  for (int t = 0; t < threads; ++t)
  {
    int a = (int)((int64_t)rows * t / threads);
    int b = (int)((int64_t)rows * (t + 1) / threads);
    pool.emplace_back([=, &fn] { fn(a, b); });
  }
  for (auto &th : pool)
    th.join();
}

Field downsample_area(const float *src,
                      int          sw,
                      int          sh,
                      int          dw,
                      int          dh,
                      int          threads)
{
  Field       tmp((size_t)dw * sh);
  Field       dst((size_t)dw * dh);
  const float rx = (float)sw / (float)dw;
  const float ry = (float)sh / (float)dh;
  parallel_rows(sh,
                threads,
                [&](int y0, int y1)
                {
                  for (int y = y0; y < y1; ++y)
                    for (int x = 0; x < dw; ++x)
                    {
                      float a = x * rx, b = (x + 1) * rx, s = 0.f;
                      for (int i = (int)a; i < sw && (float)i < b; ++i)
                        s += src[(size_t)y * sw + i] *
                             (std::min(b, (float)(i + 1)) -
                              std::max(a, (float)i));
                      tmp[(size_t)y * dw + x] = s / rx;
                    }
                });
  parallel_rows(dh,
                threads,
                [&](int y0, int y1)
                {
                  for (int y = y0; y < y1; ++y)
                  {
                    float a = y * ry, b = (y + 1) * ry;
                    for (int x = 0; x < dw; ++x)
                    {
                      float s = 0.f;
                      for (int j = (int)a; j < sh && (float)j < b; ++j)
                        s += tmp[(size_t)j * dw + x] *
                             (std::min(b, (float)(j + 1)) -
                              std::max(a, (float)j));
                      dst[(size_t)y * dw + x] = s / ry;
                    }
                  }
                });
  return dst;
}

void extrapolate_border_impl(std::vector<float> &z,
                             int                 nx,
                             int                 ny,
                             float               zmin,
                             float               zmax,
                             const uint8_t      *outlet,
                             const float        *bedrock,
                             std::vector<float> *delta)
{
  if (nx < 3 || ny < 3) return;
  auto set = [&](int x, int y, int x1, int y1, int x2, int y2)
  {
    const size_t i = (size_t)y * nx + x;
    if (outlet && outlet[i]) return;
    float v = std::min(
        zmax,
        std::max(zmin,
                 2.f * z[(size_t)y1 * nx + x1] - z[(size_t)y2 * nx + x2]));
    if (bedrock) v = std::max(v, std::min(bedrock[i], z[i]));
    if (delta && delta->size() == z.size()) (*delta)[i] += v - z[i];
    z[i] = v;
  };
  for (int y = 0; y < ny; ++y)
  {
    set(0, y, 1, y, 2, y);
    set(nx - 1, y, nx - 2, y, nx - 3, y);
  }
  for (int x = 0; x < nx; ++x)
  {
    set(x, 0, x, 1, x, 2);
    set(x, ny - 1, x, ny - 2, x, ny - 3);
  }
}

uint32_t float_key(float f)
{
  uint32_t u = 0;
  std::memcpy(&u, &f, 4);
  return (u & 0x80000000u) ? ~u : (u | 0x80000000u);
}

float key_float(uint32_t k)
{
  uint32_t u = (k & 0x80000000u) ? (k & 0x7fffffffu) : ~k;
  float    f = 0.f;
  std::memcpy(&f, &u, 4);
  return f;
}

float next_up(float f)
{
  return key_float(float_key(f) + 1u);
}

int thread_count(int requested)
{
  int t = requested > 0 ? requested : (int)std::thread::hardware_concurrency();
  return std::max(1, std::min(t, 16));
}

// monotone priority queue (radix heap)
struct RadixHeap
{
  std::vector<std::pair<uint32_t, int>> bucket[33];
  uint32_t                              last = 0;
  size_t                                count = 0;

  void clear()
  {
    for (auto &b : bucket)
      b.clear();
    last = 0;
    count = 0;
  }
  static int slot(uint32_t x)
  {
    return (int)std::bit_width(x);
  }
  void push(uint32_t key, int idx)
  {
    bucket[slot(key ^ last)].emplace_back(key, idx);
    ++count;
  }
  std::pair<uint32_t, int> pop()
  {
    if (bucket[0].empty())
    {
      int i = 1;
      while (bucket[i].empty())
        ++i;
      uint32_t mn = 0xffffffffu;
      for (auto &e : bucket[i])
        mn = std::min(mn, e.first);
      last = mn;
      for (auto &e : bucket[i])
        bucket[slot(e.first ^ last)].push_back(e);
      bucket[i].clear();
    }
    auto e = bucket[0].back();
    bucket[0].pop_back();
    --count;
    return e;
  }
};

// drainage graph structure
struct FlowGraph
{
  int                  nx = 0, ny = 0, n = 0;
  Field                zf;    // routing surface: z with the depressions filled
  std::vector<int>     order; // flood order, receivers before donors
  std::vector<int>     rank;  // inverse permutation of order
  std::vector<int>     rcv;   // receiver of each cell, itself for an outlet
  std::vector<uint8_t> flag;
  Field                area;  // drainage area in cells, or accumulated rain
  Field                steep; // upstream sum of the relative slope excess

  std::vector<int> lake;        // lake index of a flooded cell, else -1
  std::vector<int> lake_outlet; // spill cell of each lake
  std::vector<int>
      lake_root; // first flooded cell visited upstream-to-downstream
  std::vector<float>  lake_volume, lake_sediment, lake_target;
  std::vector<float>  lake_lo, lake_hi, lake_level;
  std::vector<double> lake_acc;
  std::vector<int>    pitq; // all flooded cells, in flood order
  int                 pit_cells = 0;

  RadixHeap        heap;         // exact flood
  std::vector<int> bhead, bnext; // bucket queue of the fast flood

  void resize(int nx_, int ny_)
  {
    nx = nx_;
    ny = ny_;
    n = nx * ny;
    zf.resize(n);
    order.resize(n);
    rank.resize(n);
    rcv.resize(n);
    flag.resize(n);
    area.resize(n);
    steep.resize(n);
    lake.resize(n);
    pitq.resize(n);
  }

  float dist(int i) const
  {
    int d = rcv[i] > i ? rcv[i] - i : i - rcv[i];
    return (d == 1 || d == nx) ? 1.f : 1.41421356f;
  }

  bool submerged(int i, const Field &z) const
  {
    return (flag[i] & F_PIT) && zf[i] > z[i];
  }

  void flood(const Field &z, const uint8_t *outlet, bool fast)
  {
    std::fill(flag.begin(), flag.end(), (uint8_t)0);
    std::fill(lake.begin(), lake.end(), -1);
    lake_outlet.clear();
    lake_root.clear();
    lake_volume.clear();

    const int nbuckets = 1 << 18;
    float     zmin = 0.f, scale = 0.f;
    int       cur = nbuckets, open = 0;
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
      heap.clear();

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
        heap.push(float_key(z[i]), i);
    };
    auto seed_cell = [&](int i)
    {
      if (flag[i]) return;
      flag[i] = F_CLOSED | F_OUTLET;
      zf[i] = z[i];
      push(i);
    };
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
      for (int i = 0; i < n; ++i)
        if (outlet[i]) seed_cell(i);

    size_t head = 0, tail = 0;
    int    count = 0;
    int    cur_lake = -1;
    while (head < tail || (fast ? open > 0 : heap.count > 0))
    {
      int c = 0;
      if (head < tail)
        c = pitq[head++];
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
          c = heap.pop().second;
        cur_lake = -1;
      }
      rank[c] = count;
      order[count++] = c;

      const float zc_eps = next_up(zf[c]);
      const int   cx = c % nx, cy = c / nx;
      for (int k = 0; k < 8; ++k)
      {
        const int x = cx + DX8[k], y = cy + DY8[k];
        if (x < 0 || y < 0 || x >= nx || y >= ny) continue;
        const int j = y * nx + x;
        if (flag[j] & F_CLOSED) continue;
        if (z[j] <= zc_eps)
        {
          flag[j] = F_CLOSED | F_PIT;
          zf[j] = zc_eps;
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
          zf[j] = z[j];
          push(j);
        }
      }
    }
    pit_cells = (int)tail;
  }

  void receivers(float p, uint32_t seed, int threads)
  {
    parallel_rows(ny,
                  threads,
                  [&](int y0, int y1)
                  {
                    for (int y = y0; y < y1; ++y)
                      for (int x = 0; x < nx; ++x)
                      {
                        const int i = y * nx + x;
                        rcv[i] = i;
                        if (flag[i] & F_OUTLET) continue;
                        const float zi = zf[i];
                        const int   ri = rank[i];
                        int         cj[8];
                        float       cw[8];
                        int         cnt = 0;
                        float       sum = 0.f;
                        for (int k = 0; k < 8; ++k)
                        {
                          const int xx = x + DX8[k], yy = y + DY8[k];
                          if (xx < 0 || yy < 0 || xx >= nx || yy >= ny)
                            continue;
                          const int j = yy * nx + xx;
                          if (rank[j] < ri && zf[j] < zi)
                          {
                            const float d = zi - zf[j];
                            const float w = (p == 1.f)
                                                ? d
                                                : ((p == 2.f) ? d * d
                                                              : std::pow(d, p));
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
                  });
  }

  void accumulate(const Field &z, const float *rain, float talus_cell)
  {
    if (rain)
      std::copy(rain, rain + n, area.begin());
    else
      std::fill(area.begin(), area.end(), 1.f);
    std::fill(steep.begin(), steep.end(), 0.f);
    const float inv = 1.f / talus_cell;
    for (int o = n - 1; o >= 0; --o)
    {
      const int i = order[o];
      const int r = rcv[i];
      if (r == i) continue;
      if (!submerged(i, z))
      {
        const float s = (z[i] - z[r]) / dist(i);
        if (s > talus_cell) steep[i] += (s - talus_cell) * inv;
      }
      steep[r] += steep[i];
      area[r] += area[i];
    }
  }
};

struct Level
{
  int            nx = 0, ny = 0;
  Field         *z = nullptr;
  Field         *sed = nullptr;
  const float   *erodibility = nullptr;
  const float   *rain = nullptr;
  const float   *bedrock = nullptr;
  const uint8_t *outlet = nullptr;
};

struct Scratch
{
  Field ero, qs, zdon, delta, powlut;
  void  resize(size_t n)
  {
    ero.resize(n);
    qs.resize(n);
    zdon.resize(n);
    delta.resize(n);
  }
};

void upsample3(const Field &a_src,
               const Field &b_src,
               const Field &c_src,
               int          sw,
               int          sh,
               Field       &a_dst,
               Field       &b_dst,
               Field       &c_dst,
               int          dw,
               int          dh,
               float        warp,
               uint32_t     seed,
               int          threads)
{
  a_dst.resize((size_t)dw * dh);
  b_dst.resize((size_t)dw * dh);
  c_dst.resize((size_t)dw * dh);
  const float rx = (float)sw / (float)dw, ry = (float)sh / (float)dh;
  parallel_rows(
      dh,
      threads,
      [&](int y0, int y1)
      {
        for (int y = y0; y < y1; ++y)
          for (int x = 0; x < dw; ++x)
          {
            float fx = (x + 0.5f) * rx - 0.5f;
            float fy = (y + 0.5f) * ry - 0.5f;
            if (warp > 0.f)
            {
              float u = fx * 0.5f, v = fy * 0.5f;
              fx += warp * vnoise(u + 13.7f, v + 3.1f, seed);
              fy += warp * vnoise(u + 71.3f, v + 47.9f, seed + 101u);
            }
            int   xi = (int)std::floor(fx), yi = (int)std::floor(fy);
            float wx[4], wy[4];
            bspline_weights(fx - (float)xi, wx);
            bspline_weights(fy - (float)yi, wy);
            float sa = 0.f, sb = 0.f, sc = 0.f;
            for (int j = 0; j < 4; ++j)
            {
              const size_t row = (size_t)std::min(sh - 1,
                                                  std::max(0, yi - 1 + j)) *
                                 sw;
              float ra = 0.f, rb = 0.f, rc = 0.f;
              for (int i = 0; i < 4; ++i)
              {
                const size_t k = row +
                                 (size_t)std::min(sw - 1,
                                                  std::max(0, xi - 1 + i));
                ra += wx[i] * a_src[k];
                rb += wx[i] * b_src[k];
                rc += wx[i] * c_src[k];
              }
              sa += wy[j] * ra;
              sb += wy[j] * rb;
              sc += wy[j] * rc;
            }
            const size_t o = (size_t)y * dw + x;
            a_dst[o] = sa;
            b_dst[o] = sb;
            c_dst[o] = sc;
          }
      });
}

void erode_pass(const Level      &L,
                const FlowGraph  &G,
                const MiseParams &P,
                float             dtK,
                float             dtKt,
                float             dtKd,
                Scratch          &W)
{
  Field      &z = *L.z;
  Field      &sed = *L.sed;
  const int   nx = G.nx;
  const float dx = 1.f / (float)nx;
  const float cellA = dx * dx;
  const float m = P.area_exp;

  const bool lut_ok = (L.rain == nullptr);
  const int  lut_n = 2048;
  if (lut_ok)
  {
    W.powlut.resize(lut_n);
    for (int k = 0; k < lut_n; ++k)
      W.powlut[k] = std::pow(std::min((float)k * cellA, P.downcutting), m);
  }

  for (int o = 0; o < G.n; ++o)
  {
    const int i = G.order[o];
    const int r = G.rcv[i];
    W.ero[i] = 0.f;
    if (r == i) continue;
    const float zi = z[i], zr = z[r];
    if (zi <= zr || G.submerged(i, z)) continue;

    const float d = G.dist(i);
    const float dd = d * dx;
    const float ac = G.area[i];
    const float am = (lut_ok && ac < (float)lut_n)
                         ? W.powlut[(int)ac]
                         : std::pow(std::min(ac * cellA, P.downcutting), m);
    const float ke = L.erodibility ? L.erodibility[i] : 1.f;
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

    // fluvial incision
    const float Ff = dtK * ke * soft * am * rho / dd;
    float       zn = (zi + Ff * zr) / (1.f + Ff);

    // hillslope relaxation
    const float zc = zr + P.talus * dd;
    if (zn > zc && (dtKt > 0.f || dtKd > 0.f))
    {
      const float Fh = (dtKt + dtKd * G.steep[i] * cellA) * ke * soft / dd;
      zn = (zi + Ff * zr + Fh * zc) / (1.f + Ff + Fh);
    }
    if (P.cliff > 0.f) zn = std::min(zn, zr + P.cliff * dd);
    if (L.bedrock) zn = std::min(zi, std::max(zn, L.bedrock[i]));

    const float e = zi - zn;
    W.ero[i] = e;
    z[i] = zn;
    sed[i] = std::max(0.f, sed[i] - e);
  }
}

void sediment_pass(const Level      &L,
                   FlowGraph        &G,
                   const MiseParams &P,
                   Scratch          &W)
{
  Field      &z = *L.z;
  Field      &sed = *L.sed;
  const float dx = 1.f / (float)G.nx;
  const float cellA = dx * dx;
  const float a_ref = std::max(P.deposition_area, 1e-12f);
  std::fill(W.qs.begin(), W.qs.end(), 0.f);
  std::fill(W.zdon.begin(), W.zdon.end(), 3.0e38f);
  G.lake_sediment.assign(G.lake_outlet.size(), 0.f);
  G.lake_target.assign(G.lake_outlet.size(), 0.f);

  for (int o = G.n - 1; o >= 0; --o)
  {
    const int i = G.order[o];
    const int r = G.rcv[i];
    float     q = W.qs[i];
    if (r == i) continue;
    const int lk = G.lake[i];
    if (lk >= 0)
    {
      G.lake_sediment[lk] += q + W.ero[i];
      if (i == G.lake_root[lk])
      {
        const float s = G.lake_sediment[lk];
        const float fill = P.lake_fill * std::min(s, G.lake_volume[lk]);
        G.lake_target[lk] = fill;
        if (s > fill) W.qs[G.lake_outlet[lk]] += s - fill;
      }
      continue;
    }
    if (q > 0.f)
    {
      const float an = G.area[i] * cellA;
      float       seq = P.deposition_slope;
      if (an > a_ref) seq *= std::pow(an / a_ref, -P.deposition_exp);
      const float room = std::min(z[r] + seq * G.dist(i) * dx, W.zdon[i]) -
                         z[i];
      if (room > 0.f)
      {
        const float dep = std::min(q, P.deposition_rate * room);
        z[i] += dep;
        sed[i] += dep;
        q -= dep;
      }
    }
    q += W.ero[i];
    if (q > 0.f)
    {
      W.qs[r] += q;
      W.zdon[r] = std::min(W.zdon[r], z[i]);
    }
  }

  const size_t nl = G.lake_outlet.size();
  if (nl > 0)
  {
    const int *cells = G.pitq.data();
    const int  nc = G.pit_cells;
    G.lake_lo.assign(nl, 3.0e38f);
    G.lake_hi.assign(nl, -3.0e38f);
    for (int c = 0; c < nc; ++c)
    {
      const int i = cells[c], lk = G.lake[i];
      G.lake_lo[lk] = std::min(G.lake_lo[lk], z[i]);
      G.lake_hi[lk] = std::max(G.lake_hi[lk], G.zf[i]);
    }
    std::vector<float> &lo = G.lake_lo, &hi = G.lake_hi;
    G.lake_level.resize(nl);
    for (int it = 0; it < 20; ++it)
    {
      G.lake_acc.assign(nl, 0.0);
      for (size_t k = 0; k < nl; ++k)
        G.lake_level[k] = 0.5f * (lo[k] + hi[k]);
      for (int c = 0; c < nc; ++c)
      {
        const int   i = cells[c], lk = G.lake[i];
        const float h = G.lake_level[lk] - z[i];
        if (h > 0.f) G.lake_acc[lk] += h;
      }
      for (size_t k = 0; k < nl; ++k)
        (G.lake_acc[k] > (double)G.lake_target[k] ? hi[k]
                                                  : lo[k]) = G.lake_level[k];
    }
    for (int c = 0; c < nc; ++c)
    {
      const int i = cells[c], lk = G.lake[i];
      if (G.lake_target[lk] <= 0.f) continue;
      const float level = std::min(lo[lk], G.zf[i]);
      if (level > z[i])
      {
        const float dep = level - z[i];
        z[i] = level;
        sed[i] += dep;
      }
    }
  }
}

void relax_sediment(Field &z,
                    Field &sed,
                    int    nx,
                    int    ny,
                    float  talus_cell,
                    int    iters,
                    Field &delta)
{
  for (int it = 0; it < iters; ++it)
  {
    std::fill(delta.begin(), delta.end(), 0.f);
    for (int y = 1; y < ny - 1; ++y)
      for (int x = 1; x < nx - 1; ++x)
      {
        const int i = y * nx + x;
        if (sed[i] <= 0.f) continue;
        const float zi = z[i];
        float       ex[8], sum = 0.f, mx = 0.f;
        for (int k = 0; k < 8; ++k)
        {
          float e = zi - z[i + DY8[k] * nx + DX8[k]] - talus_cell * DD8[k];
          e = e > 0.f ? e : 0.f;
          ex[k] = e;
          sum += e;
          mx = e > mx ? e : mx;
        }
        if (sum <= 0.f) continue;
        const float amount = std::min(0.25f * mx, sed[i]);
        delta[i] -= amount;
        const float s = amount / sum;
        for (int k = 0; k < 8; ++k)
          if (ex[k] > 0.f) delta[i + DY8[k] * nx + DX8[k]] += ex[k] * s;
      }
    for (int y = 1; y < ny - 1; ++y)
      for (int x = 1; x < nx - 1; ++x)
      {
        const int   i = y * nx + x;
        const float d = delta[i];
        if (d == 0.f) continue;
        z[i] += d;
        sed[i] = std::max(0.f, sed[i] + d);
      }
  }
}

void erode_level(const Level      &L,
                 const MiseParams &P,
                 float             share,
                 int               iters,
                 uint32_t          seed,
                 int               threads,
                 FlowGraph        &G,
                 Scratch          &W)
{
  if (iters <= 0 || share <= 0.f) return;
  G.resize(L.nx, L.ny);
  W.resize((size_t)G.n);
  Field      &z = *L.z;
  const float dt = share / (float)iters;
  const float talus_cell = P.talus / (float)L.nx;
  for (int it = 0; it < iters; ++it)
  {
    G.flood(z, L.outlet, !P.exact_flood);
    G.receivers(P.receiver_exp, seed, threads);
    G.accumulate(z, L.rain, talus_cell);
    erode_pass(L, G, P, P.strength * dt, P.thermal * dt, P.debris * dt, W);
    sediment_pass(L, G, P, W);
    if (P.sediment_relax_iters > 0)
      relax_sediment(z,
                     *L.sed,
                     L.nx,
                     L.ny,
                     P.sediment_talus / (float)L.nx,
                     P.sediment_relax_iters,
                     W.delta);
  }
}

Field flow_map_impl(const Field   &z,
                    int            nx,
                    int            ny,
                    const float   *rain,
                    const uint8_t *outlet,
                    uint32_t       seed,
                    int            threads)
{
  Field flow((size_t)nx * ny, 0.f);
  if (nx < 3 || ny < 3 || z.size() != flow.size() ||
      !all_finite(z.data(), z.size()))
    return flow;
  FlowGraph G;
  G.resize(nx, ny);
  G.flood(z, outlet, false);
  G.receivers(4.f, seed, thread_count(threads));
  G.accumulate(z, rain, 1.f);
  const float inv = 1.f / (float)flow.size();
  for (size_t i = 0; i < flow.size(); ++i)
    flow[i] = G.area[i] * inv;
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
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;
  if (p_erodibility && !validate_same_shape(z, *p_erodibility)) return;
  if (p_moisture_map && !validate_same_shape(z, *p_moisture_map)) return;
  if (p_outlet && !validate_same_shape(z, *p_outlet)) return;

  const int    nx = z.shape.x;
  const int    ny = z.shape.y;
  const size_t n = (size_t)nx * ny;

  if (nx < 3 || ny < 3 || z.vector.size() != n ||
      !all_finite(z.vector.data(), n))
    return;

  const Array z0 = z;
  Field       sediment(n, 0.f);
  const float zfloor = *std::min_element(z0.vector.begin(), z0.vector.end());
  const int   threads = thread_count(params.threads);

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

  const float *p_erodibility_ptr = p_erodibility ? p_erodibility->vector.data()
                                                 : nullptr;
  const float *p_rain_ptr = p_moisture_map ? p_moisture_map->vector.data()
                                           : nullptr;
  const float *p_bedrock_ptr = p_bedrock ? p_bedrock->vector.data() : nullptr;

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

      Field before;
      if (params.antialias > 0.f) before = z.vector;
      Level L{nx,
              ny,
              &z.vector,
              &sediment,
              p_erodibility_ptr,
              p_rain_ptr,
              p_bedrock_ptr,
              p_outlet_ptr};
      erode_level(L, params, share, it, seed, threads, G, W);

      if (params.antialias > 0.f)
      {
        Field d(n), tmp(n);
        for (size_t i = 0; i < n; ++i)
          d[i] = z.vector[i] - before[i];
        parallel_rows(
            ny,
            threads,
            [&](int y0, int y1)
            {
              for (int y = y0; y < y1; ++y)
                for (int x = 0; x < nx; ++x)
                {
                  const int xm = std::max(0, x - 1),
                            xp = std::min(nx - 1, x + 1);
                  tmp[(size_t)y * nx + x] = 0.25f * (d[(size_t)y * nx + xm] +
                                                     d[(size_t)y * nx + xp]) +
                                            0.5f * d[(size_t)y * nx + x];
                }
            });
        parallel_rows(
            ny,
            threads,
            [&](int y0, int y1)
            {
              for (int y = y0; y < y1; ++y)
              {
                const int ym = std::max(0, y - 1), yp = std::min(ny - 1, y + 1);
                for (int x = 0; x < nx; ++x)
                {
                  const size_t i = (size_t)y * nx + x;
                  const float  b = 0.25f * (tmp[(size_t)ym * nx + x] +
                                           tmp[(size_t)yp * nx + x]) +
                                  0.5f * tmp[i];
                  z.vector[i] = std::max(zfloor,
                                         before[i] + d[i] +
                                             params.antialias * (b - d[i]));
                }
              }
            });
      }
    }
    else
    {
      // --- Coarse resolution level

      Field zl = downsample_area(z.vector.data(), nx, ny, lx, ly, threads);
      Field sl = downsample_area(sediment.data(), nx, ny, lx, ly, threads);
      Field el, rl, bl;
      std::vector<uint8_t> ol;
      if (p_erodibility_ptr)
        el = downsample_area(p_erodibility_ptr, nx, ny, lx, ly, threads);
      if (p_rain_ptr) rl = downsample_area(p_rain_ptr, nx, ny, lx, ly, threads);
      if (p_bedrock_ptr)
        bl = downsample_area(p_bedrock_ptr, nx, ny, lx, ly, threads);
      if (p_outlet_ptr)
      {
        Field of(n);
        for (size_t i = 0; i < n; ++i)
          of[i] = p_outlet_ptr[i] ? 1.f : 0.f;
        const Field oc = downsample_area(of.data(), nx, ny, lx, ly, threads);
        ol.resize(oc.size());
        for (size_t i = 0; i < oc.size(); ++i)
          ol[i] = oc[i] > 0.5f;
      }
      Field db_l(zl.size());
      for (size_t i = 0; i < zl.size(); ++i)
        db_l[i] = zl[i] - sl[i];

      Level L{lx,
              ly,
              &zl,
              &sl,
              p_erodibility_ptr ? el.data() : nullptr,
              p_rain_ptr ? rl.data() : nullptr,
              p_bedrock_ptr ? bl.data() : nullptr,
              p_outlet_ptr ? ol.data() : nullptr};
      erode_level(L, params, share, it, seed, threads, G, W);

      for (size_t i = 0; i < zl.size(); ++i)
        db_l[i] = (zl[i] - sl[i]) - db_l[i];
      Field db, zs, ss;
      upsample3(db_l,
                zl,
                sl,
                lx,
                ly,
                db,
                zs,
                ss,
                nx,
                ny,
                params.warp,
                seed,
                threads);
      parallel_rows(
          ny,
          threads,
          [&](int y0, int y1)
          {
            for (size_t i = (size_t)y0 * nx; i < (size_t)y1 * nx; ++i)
            {
              const float b1 = zs[i] - ss[i];
              const float detail = (z.vector[i] - sediment[i]) - (b1 - db[i]);
              const float ad = std::fabs(detail), cut = std::fabs(db[i]);
              const float b = b1 + (ad > 0.f ? detail * ad / (ad + cut) : 0.f);
              const float s = std::min(std::max(zs[i] - b, 0.f),
                                       2.f * std::max(ss[i], 0.f));
              z.vector[i] = std::max(b + s, zfloor);
              sediment[i] = s;
            }
          });
      if (p_bedrock_ptr)
        for (size_t i = 0; i < n; ++i)
          z.vector[i] = std::max(z.vector[i],
                                 std::min(z0.vector[i], p_bedrock_ptr[i]));
    }
    share *= params.share_ratio;
    iters *= params.iters_decay;
  }

  // restore fixed outlet cells
  auto restore = [&](size_t i)
  {
    z.vector[i] = z0.vector[i];
    sediment[i] = 0.f;
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
    for (size_t i = 0; i < n; ++i)
      if (p_outlet_ptr[i]) restore(i);
  if (p_bedrock_ptr)
    for (size_t i = 0; i < n; ++i)
      z.vector[i] = std::max(z.vector[i],
                             std::min(z0.vector[i], p_bedrock_ptr[i]));

  if (params.extrapolate_border)
    extrapolate_border_impl(
        z.vector,
        nx,
        ny,
        zfloor,
        *std::max_element(z0.vector.begin(), z0.vector.end()),
        p_outlet_ptr,
        p_bedrock_ptr,
        nullptr);

  if (p_sediment)
  {
    *p_sediment = Array(z.shape);
    p_sediment->vector = sediment;
  }

  if (p_flow)
  {
    *p_flow = Array(z.shape);
    p_flow->vector = flow_map_impl(z.vector,
                                   nx,
                                   ny,
                                   p_rain_ptr,
                                   p_outlet_ptr,
                                   params.seed,
                                   params.threads);
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
