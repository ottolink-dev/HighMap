/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <queue>
#include <vector>

#include "highmap/algebra.hpp"
#include "highmap/array.hpp"
#include "highmap/erosion.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/array.hpp"
#include "highmap/primitives/functions.hpp"
#include "highmap/range.hpp"
#include "highmap/virtual_array/tile_region.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

namespace hmap
{

namespace
{

// 8-neighbor Moore pattern offsets and distances (normalized for unit spacing)
// clang-format off
constexpr int   DX[8] = {1, 1, 0, -1, -1, -1, 0, 1};
constexpr int   DY[8] = {0, -1, -1, -1, 0, 1, 1, 1};
constexpr float DIST[8] = {1.f, (float)M_SQRT2, 1.f, (float)M_SQRT2, 1.f, (float)M_SQRT2, 1.f, (float)M_SQRT2};
// clang-format on

// Solves tridiagonal linear system using Thomas' algorithm (TDMA)
void solve_tridiagonal(const std::vector<float> &lower,
                       const std::vector<float> &diag,
                       const std::vector<float> &upper,
                       const std::vector<float> &rhs,
                       std::vector<float>       &result,
                       int                       n)
{
  if (n <= 0) return;
  std::vector<float> gam(n, 0.f);

  float bet = diag[0];
  if (std::abs(bet) < 1e-12f)
  {
    result = rhs;
    return;
  }
  result[0] = rhs[0] / bet;

  for (int i = 1; i < n; ++i)
  {
    gam[i] = upper[i - 1] / bet;
    bet = diag[i] - lower[i] * gam[i];
    if (std::abs(bet) < 1e-12f)
    {
      result = rhs;
      return;
    }
    result[i] = (rhs[i] - lower[i] * result[i - 1]) / bet;
  }

  for (int i = n - 2; i >= 0; --i)
  {
    result[i] -= gam[i + 1] * result[i + 1];
  }
}

// 2D Hillslope diffusion via Alternating Direction Implicit (ADI) scheme
void hillslope_diffusion_adi(Array &z, float k_diff, float dt)
{
  if (k_diff <= 0.f || dt <= 0.f) return;

  const int nx = z.shape.x;
  const int ny = z.shape.y;
  if (nx < 3 || ny < 3) return;

  const float fx = 0.5f * k_diff;
  const float fy = 0.5f * k_diff;

  Array z_tmp = z;

  // Pass 1: solve implicitly along X (rows)
  {
    std::vector<float> lower(nx, -fx * dt);
    std::vector<float> diag(nx, 1.f + 2.f * fx * dt);
    std::vector<float> upper(nx, -fx * dt);
    std::vector<float> rhs(nx, 0.f);
    std::vector<float> row_out(nx, 0.f);

    lower[0] = 0.f;
    diag[0] = 1.f;
    upper[0] = 0.f;
    lower[nx - 1] = 0.f;
    diag[nx - 1] = 1.f;
    upper[nx - 1] = 0.f;

    for (int j = 1; j < ny - 1; ++j)
    {
      rhs[0] = z(0, j);
      rhs[nx - 1] = z(nx - 1, j);

      for (int i = 1; i < nx - 1; ++i)
      {
        rhs[i] = (1.f - 2.f * fy * dt) * z(i, j) + fy * dt * z(i, j - 1) +
                 fy * dt * z(i, j + 1);
      }

      solve_tridiagonal(lower, diag, upper, rhs, row_out, nx);

      for (int i = 0; i < nx; ++i)
      {
        z_tmp(i, j) = row_out[i];
      }
    }
  }

  // Pass 2: solve implicitly along Y (cols)
  {
    std::vector<float> lower(ny, -fy * dt);
    std::vector<float> diag(ny, 1.f + 2.f * fy * dt);
    std::vector<float> upper(ny, -fy * dt);
    std::vector<float> rhs(ny, 0.f);
    std::vector<float> col_out(ny, 0.f);

    lower[0] = 0.f;
    diag[0] = 1.f;
    upper[0] = 0.f;
    lower[ny - 1] = 0.f;
    diag[ny - 1] = 1.f;
    upper[ny - 1] = 0.f;

    for (int i = 1; i < nx - 1; ++i)
    {
      rhs[0] = z_tmp(i, 0);
      rhs[ny - 1] = z_tmp(i, ny - 1);

      for (int j = 1; j < ny - 1; ++j)
      {
        rhs[j] = (1.f - 2.f * fx * dt) * z_tmp(i, j) +
                 fx * dt * z_tmp(i - 1, j) + fx * dt * z_tmp(i + 1, j);
      }

      solve_tridiagonal(lower, diag, upper, rhs, col_out, ny);

      for (int j = 0; j < ny; ++j)
      {
        z(i, j) = col_out[j];
      }
    }
  }
}

} // namespace

// --- FastScape landscape evolution erosion algorithm

void hydraulic_fastscape(Array       &z,
                         int          iterations,
                         float        dt,
                         float        k_erosion,
                         float        m_exp,
                         float        n_exp,
                         float        k_diff,
                         float        uplift_rate,
                         bool         multiple_flow,
                         float        flow_partition_exp,
                         float        tolerance,
                         const Array *p_bedrock,
                         const Array *p_moisture_map,
                         Array       *p_erosion_map,
                         Array       *p_flow_map)
{
  if (!validate_non_empty(z)) return;
  if (p_bedrock && !validate_same_shape(z, *p_bedrock)) return;
  if (p_moisture_map && !validate_same_shape(z, *p_moisture_map)) return;

  const int nx = z.shape.x;
  const int ny = z.shape.y;
  const int n_nodes = nx * ny;

  Array z_initial;
  if (p_erosion_map) z_initial = z;

  Array flow_acc(z.shape, 0.f);

  const bool is_linear = (std::abs(n_exp - 1.f) < 1e-4f);

  struct QNode
  {
    int   i, j;
    float elev;
  };

  struct QNodeCmp
  {
    bool operator()(const QNode &a, const QNode &b) const
    {
      return a.elev > b.elev; // min-heap
    }
  };

  auto flat_idx = [nx](int i, int j) -> int { return j * nx + i; };

  // Graph structure containers
  std::vector<int>   single_receiver(n_nodes, 0);
  std::vector<float> single_dist(n_nodes, 0.f);

  // MFD containers (up to 8 receivers per cell)
  std::vector<std::array<int, 8>>   mfd_receivers(n_nodes);
  std::vector<std::array<float, 8>> mfd_weights(n_nodes);
  std::vector<std::array<float, 8>> mfd_distances(n_nodes);
  std::vector<uint8_t>              mfd_rec_count(n_nodes, 0);

  std::vector<int>   ordered_nodes(n_nodes, 0);
  std::vector<float> drainage_area(n_nodes, 0.f);

  for (int iter = 0; iter < iterations; ++iter)
  {
    // --- 1. Uplift
    if (uplift_rate != 0.f)
    {
      for (int j = 1; j < ny - 1; ++j)
        for (int i = 1; i < nx - 1; ++i)
          z(i, j) += uplift_rate * dt;
    }

    // --- 2. Flow routing and sink resolution via Priority-Flood
    std::priority_queue<QNode, std::vector<QNode>, QNodeCmp> pq;
    std::vector<uint8_t> visited(n_nodes, 0);

    int order_count = 0;

    // Boundary cells are fixed base-level outlets
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
      {
        if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
        {
          int id = flat_idx(i, j);
          single_receiver[id] = id;
          single_dist[id] = 0.f;
          visited[id] = 1;
          pq.push({i, j, z(i, j)});
        }
      }

    while (!pq.empty())
    {
      QNode curr = pq.top();
      pq.pop();

      int curr_id = flat_idx(curr.i, curr.j);
      ordered_nodes[order_count++] = curr_id;

      for (int k = 0; k < 8; ++k)
      {
        int ni = curr.i + DX[k];
        int nj = curr.j + DY[k];

        if (ni < 0 || ni >= nx || nj < 0 || nj >= ny) continue;

        int nid = flat_idx(ni, nj);
        if (!visited[nid])
        {
          visited[nid] = 1;
          single_receiver[nid] = curr_id;
          single_dist[nid] = DIST[k];

          float flood_elev = std::max(z(ni, nj), curr.elev);
          pq.push({ni, nj, flood_elev});
        }
      }
    }

    // --- 3. Build multiple flow routing if enabled
    if (multiple_flow)
    {
      for (int j = 0; j < ny; ++j)
      {
        for (int i = 0; i < nx; ++i)
        {
          int id = flat_idx(i, j);

          if (i == 0 || i == nx - 1 || j == 0 || j == ny - 1)
          {
            mfd_rec_count[id] = 1;
            mfd_receivers[id][0] = id;
            mfd_weights[id][0] = 0.f;
            mfd_distances[id][0] = 0.f;
            continue;
          }

          float   z_cur = z(i, j);
          float   weights_sum = 0.f;
          uint8_t count = 0;

          for (int k = 0; k < 8; ++k)
          {
            int ni = i + DX[k];
            int nj = j + DY[k];
            if (ni < 0 || ni >= nx || nj < 0 || nj >= ny) continue;

            float z_nbr = z(ni, nj);
            if (z_cur > z_nbr)
            {
              float dist = DIST[k];
              float slope = (z_cur - z_nbr) / dist;
              float w = std::pow(slope, flow_partition_exp);

              mfd_receivers[id][count] = flat_idx(ni, nj);
              mfd_distances[id][count] = dist;
              mfd_weights[id][count] = w;
              weights_sum += w;
              count++;
            }
          }

          if (count == 0 || weights_sum <= 0.f)
          {
            // Fallback to priority flood resolved receiver
            mfd_rec_count[id] = 1;
            mfd_receivers[id][0] = single_receiver[id];
            mfd_weights[id][0] = 1.f;
            mfd_distances[id][0] = single_dist[id];
          }
          else
          {
            mfd_rec_count[id] = count;
            for (uint8_t c = 0; c < count; ++c)
            {
              mfd_weights[id][c] /= weights_sum;
            }
          }
        }
      }
    }

    // --- 4. Drainage area accumulation (top-down / reverse topological order)
    std::fill(drainage_area.begin(), drainage_area.end(), 1.f);
    for (int idx = n_nodes - 1; idx >= 0; --idx)
    {
      int   u = ordered_nodes[idx];
      float precip = 1.f;
      if (p_moisture_map)
      {
        int ui = u % nx;
        int uj = u / nx;
        precip = (*p_moisture_map)(ui, uj);
      }

      if (multiple_flow)
      {
        uint8_t nrec = mfd_rec_count[u];
        for (uint8_t r = 0; r < nrec; ++r)
        {
          int rec = mfd_receivers[u][r];
          if (rec != u)
          {
            drainage_area[rec] += drainage_area[u] * mfd_weights[u][r] * precip;
          }
        }
      }
      else
      {
        int r = single_receiver[u];
        if (r != u)
        {
          drainage_area[r] += drainage_area[u] * precip;
        }
      }
    }

    // --- 5. Stream Power Law implicit solve (bottom-up / topological order)
    for (int idx = 0; idx < n_nodes; ++idx)
    {
      int u = ordered_nodes[idx];
      int ui = u % nx;
      int uj = u / nx;

      if (ui == 0 || ui == nx - 1 || uj == 0 || uj == ny - 1)
      {
        // Base-level boundary node: no erosion
        continue;
      }

      float z_u = z(ui, uj);
      float k_val = k_erosion;
      if (p_moisture_map) k_val *= (*p_moisture_map)(ui, uj);
      float area = drainage_area[u];

      if (multiple_flow && is_linear)
      {
        // Multi-direction flow SPL (linear case n = 1)
        uint8_t nrec = mfd_rec_count[u];
        float   num = z_u;
        float   den = 1.f;
        float   min_rec_elev = std::numeric_limits<float>::max();

        for (uint8_t r = 0; r < nrec; ++r)
        {
          int   rec = mfd_receivers[u][r];
          int   ri = rec % nx;
          int   rj = rec / nx;
          float z_r = z(ri, rj);
          min_rec_elev = std::min(min_rec_elev, z_r);

          if (rec == u) continue;

          float w = mfd_weights[u][r];
          float dist = mfd_distances[u][r];
          float factor = (k_val * dt * std::pow(area * w, m_exp)) / dist;

          num += factor * z_r;
          den += factor;
        }

        float z_u_next = num / den;
        if (min_rec_elev < std::numeric_limits<float>::max() &&
            z_u_next < min_rec_elev)
        {
          z_u_next = min_rec_elev + 1e-6f;
        }
        z(ui, uj) = z_u_next;
      }
      else
      {
        // Single flow routing or non-linear solver
        int r = single_receiver[u];
        if (r == u) continue;

        int   ri = r % nx;
        int   rj = r / nx;
        float z_r_next = z(ri, rj);
        float length = single_dist[u];

        if (z_u <= z_r_next) continue;

        float factor = k_val * dt * std::pow(area, m_exp);
        float z_u_next = z_u;

        if (is_linear)
        {
          factor /= length;
          float num = z_u + factor * z_r_next;
          float den = 1.f + factor;
          z_u_next = num / den;
        }
        else
        {
          factor /= std::pow(length, n_exp);
          double delta_0 = (double)z_u - (double)z_r_next;
          double delta_k = delta_0;

          // Newton-Raphson iteration
          for (int iter_nr = 0; iter_nr < 50; ++iter_nr)
          {
            double factor_delta_exp = (double)factor *
                                      std::pow(delta_k, (double)n_exp);
            double func = delta_k + factor_delta_exp - delta_0;

            if (std::abs(func) <= (double)tolerance) break;

            double func_deriv = 1.0 +
                                (double)n_exp * factor_delta_exp / delta_k;
            delta_k -= func / func_deriv;

            if (delta_k <= 0.0)
            {
              delta_k = 0.0;
              break;
            }
          }
          z_u_next = (float)((double)z_r_next + delta_k);
        }

        if (z_u_next < z_r_next)
        {
          z_u_next = z_r_next + 1e-6f;
        }

        z(ui, uj) = z_u_next;
      }
    }

    // --- 6. Hillslope linear diffusion
    if (k_diff > 0.f)
    {
      hillslope_diffusion_adi(z, k_diff, dt);
    }

    // --- 7. Bedrock constraint
    if (p_bedrock)
    {
      z = maximum(*p_bedrock, z);
    }
  }

  // --- Splatmap outputs
  if (p_flow_map)
  {
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i)
        flow_acc(i, j) = drainage_area[flat_idx(i, j)];
    *p_flow_map = flow_acc;
  }

  if (p_erosion_map)
  {
    *p_erosion_map = z_initial - z;
    clamp_min(*p_erosion_map, 0.f);
  }
}

void hydraulic_fastscape(Array       &z,
                         const Array *p_mask,
                         int          iterations,
                         float        dt,
                         float        k_erosion,
                         float        m_exp,
                         float        n_exp,
                         float        k_diff,
                         float        uplift_rate,
                         bool         multiple_flow,
                         float        flow_partition_exp,
                         float        tolerance,
                         const Array *p_bedrock,
                         const Array *p_moisture_map,
                         Array       *p_erosion_map,
                         Array       *p_flow_map)
{
  if (!validate_non_empty(z)) return;
  if (p_mask && !validate_same_shape(z, *p_mask)) return;

  if (!p_mask)
  {
    hydraulic_fastscape(z,
                        iterations,
                        dt,
                        k_erosion,
                        m_exp,
                        n_exp,
                        k_diff,
                        uplift_rate,
                        multiple_flow,
                        flow_partition_exp,
                        tolerance,
                        p_bedrock,
                        p_moisture_map,
                        p_erosion_map,
                        p_flow_map);
  }
  else
  {
    Array z_eroded = z;
    hydraulic_fastscape(z_eroded,
                        iterations,
                        dt,
                        k_erosion,
                        m_exp,
                        n_exp,
                        k_diff,
                        uplift_rate,
                        multiple_flow,
                        flow_partition_exp,
                        tolerance,
                        p_bedrock,
                        p_moisture_map,
                        p_erosion_map,
                        p_flow_map);
    z = lerp(z, z_eroded, *p_mask);
  }
}

} // namespace hmap

namespace hmap::va
{

VirtualArray hydraulic_fastscape(const ComputeMode  &cm,
                                 const VirtualArray &z,
                                 int                 iterations,
                                 float               dt,
                                 float               k_erosion,
                                 float               m_exp,
                                 float               n_exp,
                                 float               k_diff,
                                 float               uplift_rate,
                                 bool                multiple_flow,
                                 float               flow_partition_exp,
                                 float               tolerance,
                                 const VirtualArray *p_bedrock,
                                 const VirtualArray *p_moisture_map,
                                 VirtualArray       *p_erosion_map,
                                 VirtualArray       *p_flow_map,
                                 const VirtualArray *p_mask)
{
  VirtualArray out;
  if (z.empty()) return out;

  out.copy_from(z, cm);

  hmap::for_each_tile(
      {&z, p_mask, p_bedrock, p_moisture_map},
      {&out, p_erosion_map, p_flow_map},
      [&](std::vector<const hmap::Array *> p_arrays_in,
          std::vector<hmap::Array *>       p_arrays_out,
          const hmap::TileRegion &)
      {
        auto [pa_z, pa_mask, pa_bedrock, pa_moisture] = unpack<4>(p_arrays_in);
        auto [pa_out, pa_erosion, pa_flow] = unpack<3>(p_arrays_out);

        *pa_out = *pa_z;
        hmap::hydraulic_fastscape(*pa_out,
                                  pa_mask,
                                  iterations,
                                  dt,
                                  k_erosion,
                                  m_exp,
                                  n_exp,
                                  k_diff,
                                  uplift_rate,
                                  multiple_flow,
                                  flow_partition_exp,
                                  tolerance,
                                  pa_bedrock,
                                  pa_moisture,
                                  pa_erosion,
                                  pa_flow);
      },
      cm);

  return out;
}

} // namespace hmap::va
