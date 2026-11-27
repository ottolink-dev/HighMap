/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software.
 *
 * Adapted from MISE (Multiscale Implicit Stream-power Erosion) by Leonhard
 * https://github.com/Leonhardmaster2 */

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <vector>

#include "highmap/hydrology/hydrology.hpp"
#include "highmap/math/array.hpp"
#include "highmap/random.hpp"

namespace hmap
{

namespace
{

constexpr int DX8[8] = {1, 1, 0, -1, -1, -1, 0, 1};
constexpr int DY8[8] = {0, 1, 1, 1, 0, -1, -1, -1};

} // namespace

void FlowGraph::resize(int nx_, int ny_)
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

float FlowGraph::dist(int i) const
{
  int d = rcv[i] > i ? rcv[i] - i : i - rcv[i];
  return (d == 1 || d == nx) ? 1.f : 1.41421356f;
}

bool FlowGraph::submerged(int i, const Array &z) const
{
  return (flag[i] & F_PIT) && zf.vector[i] > z.vector[i];
}

void FlowGraph::flood(const Array   &z_arr,
                      const uint8_t *outlet,
                      bool           fast,
                      bool           seed_borders)
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

  if (seed_borders)
  {
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

void FlowGraph::receivers(float p, uint32_t seed)
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

void FlowGraph::accumulate(const Array &z_arr,
                           const Array *p_rain,
                           float        talus_cell)
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

} // namespace hmap
