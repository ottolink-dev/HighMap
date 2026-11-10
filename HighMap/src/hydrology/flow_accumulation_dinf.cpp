/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/boundary.hpp"
#include "highmap/functions.hpp"
#include "highmap/geometry/grids.hpp"
#include "highmap/gradient.hpp"
#include "highmap/hydrology/hydrology.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/math/core.hpp"
#include "highmap/primitives/coherent_noise.hpp"
#include "highmap/primitives/functions.hpp"
#include "highmap/range.hpp"
#include "highmap/virtual_array/communicator.hpp"
#include "highmap/virtual_array/tile_region.hpp"

namespace hmap
{

Array flow_accumulation_dinf(const Array      &z,
                             float             talus_ref,
                             int               iterations,
                             const TileRegion &region,
                             ICommunicator    &comm)
{
  if (!validate_non_empty(z)) return Array();

  const glm::ivec2 &shape = z.shape;
  const auto       &neighbors = neighborhood::DINF_8;
  const size_t      nb = neighbors.size();
  const int         nx = shape.x;
  const int         ny = shape.y;
  const int         ncells = nx * ny;

  // --- Fast path for single tile / non-distributed execution

  if (comm.size() <= 1)
  {
    Array  facc = constant(shape, 1.f);
    float *facc_ptr = facc.vector.data();

    std::vector<float> dinf = flow_direction_dinf_flat(z, talus_ref);

    // build nidp — parallelized with atomic increments
    std::vector<int> nidp(ncells, 0);

#pragma omp parallel for schedule(static)
    for (int j = 1; j < ny - 1; ++j)
      for (int i = 1; i < nx - 1; ++i)
      {
        const int base = (j * nx + i) * nb;
        for (size_t k = 0; k < nb; ++k)
        {
          if (dinf[base + k] > 0.f)
          {
            const int nidx = (j + neighbors[k].offset.y) * nx +
                             (i + neighbors[k].offset.x);
#pragma omp atomic
            nidp[nidx]++;
          }
        }
      }

    // seed queue with source cells (nidp == 0)
    std::vector<int> queue;
    queue.reserve(ncells / 4);

    for (int j = 1; j < ny - 1; ++j)
      for (int i = 1; i < nx - 1; ++i)
        if (nidp[j * nx + i] == 0) queue.push_back(j * nx + i);

    // topological accumulation, inherently sequential
    while (!queue.empty())
    {
      const int idx = queue.back();
      queue.pop_back();

      const int   i = idx % nx;
      const int   j = idx / nx;
      const float acc = facc_ptr[idx];
      const int   base = idx * nb;

      for (size_t k = 0; k < nb; ++k)
      {
        const float wgt = dinf[base + k];
        if (wgt == 0.f) continue;

        const int nidx = (j + neighbors[k].offset.y) * nx +
                         (i + neighbors[k].offset.x);
        facc_ptr[nidx] += acc * wgt;

        if (--nidp[nidx] == 0) queue.push_back(nidx);
      }
    }

    fill_borders(facc);
    return facc;
  }

  // Multi-tile distributed execution with communicator synchronization

  // Identify internal core bounds vs halo
  const int h_w = region.halo.x;
  const int h_e = region.halo.y;
  const int h_s = region.halo.z;
  const int h_n = region.halo.w;

  // Initialize accumulation: 1.0 on core cells, 0.0 on incoming halo cells
  Array facc(shape, 0.f);
  for (int j = h_s; j < ny - h_n; ++j)
    for (int i = h_w; i < nx - h_e; ++i)
      facc(i, j) = 1.f;

  float *facc_ptr = facc.vector.data();

  std::vector<float> dinf = flow_direction_dinf_flat(z, talus_ref);

  std::vector<int> nidp_base(ncells, 0);

#pragma omp parallel for schedule(static)
  for (int j = 1; j < ny - 1; ++j)
    for (int i = 1; i < nx - 1; ++i)
    {
      const int base = (j * nx + i) * nb;
      for (size_t k = 0; k < nb; ++k)
      {
        if (dinf[base + k] > 0.f)
        {
          const int nidx = (j + neighbors[k].offset.y) * nx +
                           (i + neighbors[k].offset.x);
#pragma omp atomic
          nidp_base[nidx]++;
        }
      }
    }

  // Determine number of sweeping iterations
  int n_sweeps = iterations;
  if (n_sweeps <= 0)
  {
    int side = std::max(1, int(std::ceil(std::sqrt(comm.size()))));
    n_sweeps = std::max(2, 2 * side);
  }

  // Delta to propagate downstream in this tile.
  // Initially, every cell's local 1.f is the delta to accumulate.
  Array  delta_acc = facc;
  float *delta_ptr = delta_acc.vector.data();

  for (int sweep = 0; sweep < n_sweeps; ++sweep)
  {
    // Local topological accumulation of delta_acc into facc
    std::vector<int> nidp = nidp_base;
    std::vector<int> queue;
    queue.reserve(ncells / 4);

    for (int j = 1; j < ny - 1; ++j)
      for (int i = 1; i < nx - 1; ++i)
        if (nidp[j * nx + i] == 0) queue.push_back(j * nx + i);

    while (!queue.empty())
    {
      const int idx = queue.back();
      queue.pop_back();

      const int   i = idx % nx;
      const int   j = idx / nx;
      const float d_val = delta_ptr[idx];
      const int   base = idx * nb;

      for (size_t k = 0; k < nb; ++k)
      {
        const float wgt = dinf[base + k];
        if (wgt == 0.f) continue;

        const int nidx = (j + neighbors[k].offset.y) * nx +
                         (i + neighbors[k].offset.x);
        if (d_val > 0.f)
        {
          const float added = d_val * wgt;
          facc_ptr[nidx] += added;
          delta_ptr[nidx] += added;
        }

        if (--nidp[nidx] == 0) queue.push_back(nidx);
      }
    }

    if (sweep == n_sweeps - 1) break;

    // Synchronize full facc across overlap boundaries using Max
    Array prev_facc = facc;
    comm.sync_overlap_buffers(facc, region, SyncOperation::Max);

    // Any positive increase on incoming halo cells becomes the new delta to
    // propagate
    delta_acc = Array(shape, 0.f);
    delta_ptr = delta_acc.vector.data();

    // Western halo cells (incoming from West)
    if (h_w > 0)
    {
      for (int j = 0; j < ny; ++j)
        for (int i = 0; i < h_w; ++i)
        {
          float diff = facc(i, j) - prev_facc(i, j);
          if (diff > 1e-4f) delta_acc(i, j) = diff;
        }
    }

    // Eastern halo cells (incoming from East)
    if (h_e > 0)
    {
      for (int j = 0; j < ny; ++j)
        for (int i = nx - h_e; i < nx; ++i)
        {
          float diff = facc(i, j) - prev_facc(i, j);
          if (diff > 1e-4f) delta_acc(i, j) = diff;
        }
    }

    // Southern halo cells (incoming from South)
    if (h_s > 0)
    {
      for (int j = 0; j < h_s; ++j)
        for (int i = 0; i < nx; ++i)
        {
          float diff = facc(i, j) - prev_facc(i, j);
          if (diff > 1e-4f) delta_acc(i, j) = diff;
        }
    }

    // Northern halo cells (incoming from North)
    if (h_n > 0)
    {
      for (int j = ny - h_n; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
        {
          float diff = facc(i, j) - prev_facc(i, j);
          if (diff > 1e-4f) delta_acc(i, j) = diff;
        }
    }
  }

  // Smooth blend along the overlap seams for clean tile stitching
  comm.sync_overlap_buffers(facc, region, SyncOperation::SmoothBlend);

  fill_borders(facc);
  return facc;
}

Array flow_accumulation_dinf_perturbed(const Array      &z,
                                       float             talus_ref,
                                       int               nsamples,
                                       glm::vec2         kw,
                                       std::uint32_t     seed,
                                       float             amp,
                                       const Array      *p_perturb_scaling,
                                       glm::vec4         bbox,
                                       int               iterations,
                                       const TileRegion &region,
                                       ICommunicator    &comm)
{
  if (!validate_non_empty(z)) return Array();
  if (p_perturb_scaling && !validate_same_shape(z, *p_perturb_scaling))
    return Array();

  const glm::ivec2 shape = z.shape;
  Array            facc(shape);

  // If a valid region bbox is provided, prefer it for coherent noise across
  // tiles
  glm::vec4 noise_bbox = (region.shape.x > 0 && region.shape.y > 0)
                             ? region.bbox
                             : bbox;

  for (int n = 0; n < nsamples; ++n)
  {
    Array zp = z;
    Array dz = noise(hmap::NoiseType::PERLIN,
                     shape,
                     kw,
                     seed + n,
                     nullptr,
                     nullptr,
                     noise_bbox);

    if (p_perturb_scaling) dz *= *p_perturb_scaling;
    zp += amp * dz;

    facc += flow_accumulation_dinf(zp, talus_ref, iterations, region, comm);
  }

  return facc / float(nsamples);
}

std::vector<Array> flow_direction_dinf(const Array &z, float talus_ref)
{
  if (!validate_non_empty(z)) return {};

  const auto  &neighbors = neighborhood::DINF_8;
  const size_t nb = neighbors.size();

  // the flow-partition exponent is defined locally based on the local
  // talus in [1, 10] (Qin et al 2007)
  Array p = Array(z.shape);
  {
    Array talus = gradient_talus(z) / talus_ref;
    clamp(talus, 0.f, 1.f);
    p = 10.f * talus + 1.f;
    clamp(p, 1.f, 3.f);
  }

  // memory consuming... every 8 direction needs a full array
  std::vector<Array> dinf(nb, {z.shape});

  for (int j = 1; j < z.shape.y - 1; j++)
    for (int i = 1; i < z.shape.x - 1; i++)
    {
      for (size_t k = 0; k < nb; k++)
      {
        float dz = z(i, j) -
                   z(i + neighbors[k].offset.x, j + neighbors[k].offset.y);
        if (dz > 0.f)
          dinf[k](i, j) = std::pow(dz * neighbors[k].inv_distance, p(i, j)) *
                          neighbors[k].weight;
      }

      // normalize
      float sum = 0.f;
      for (size_t k = 0; k < nb; k++)
        sum += dinf[k](i, j);

      if (sum > 0.f)
        for (size_t k = 0; k < nb; k++)
          dinf[k](i, j) /= sum;
    }

  return dinf;
}

std::vector<float> flow_direction_dinf_flat(const Array &z, float talus_ref)
{
  if (!validate_non_empty(z)) return {};

  const glm::ivec2 &shape = z.shape;
  const auto       &neighbors = neighborhood::DINF_8;
  const size_t      nb = neighbors.size();

  // precompute log(inv_distance) once — eliminates one log() per neighbor per
  // cell
  float log_c[8];
  for (size_t k = 0; k < nb; ++k)
    log_c[k] = neighbors[k].inv_distance > 0.f
                   ? std::log(neighbors[k].inv_distance)
                   : -std::numeric_limits<float>::infinity();

  std::vector<float> dinf(shape.x * shape.y * nb, 0.f);

  Array talus = gradient_talus(z);
  talus /= talus_ref;
  clamp_max(talus, 1.f);

#pragma omp parallel for schedule(static)
  for (int j = 1; j < shape.y - 1; ++j)
    for (int i = 1; i < shape.x - 1; ++i)
    {
      const float pij = std::clamp(10.f * talus(i, j) + 1.f, 1.f, 3.f);
      const float zij = z(i, j);

      float sum = 0.f;
      float tmp[8];

      for (size_t k = 0; k < nb; ++k)
      {
        float dz = zij -
                   z(i + neighbors[k].offset.x, j + neighbors[k].offset.y);
        if (dz > 0.f)
        {
          const float v = fast_exp(pij * (fast_log(dz) + log_c[k])) *
                          neighbors[k].weight;
          tmp[k] = v;
          sum += v;
        }
        else
          tmp[k] = 0.f;
      }

      if (sum > 0.f)
      {
        const float inv = 1.f / sum;
        const int   idx = (j * shape.x + i) * nb;
        for (size_t k = 0; k < nb; ++k)
          dinf[idx + k] = tmp[k] * inv;
      }
    }

  return dinf;
}

Array flow_direction_dinf_angle(const Array &z, float talus_ref)
{
  if (!validate_non_empty(z)) return Array();

  const auto  &neighbors = neighborhood::DINF_8;
  const size_t nb = neighbors.size();

  // flow-partition exponent (Qin et al. 2007)
  Array p = Array(z.shape);
  {
    Array talus = gradient_talus(z) / talus_ref;
    clamp_max(talus, 1.f);
    p = 10.f * talus + 1.f;
    clamp(p, 1.f, 3.f);
  }

  Array angle(z.shape, 0.f); // dominant flow direction in radians

  std::vector<float> cell_angles;

  for (size_t k = 0; k < nb; k++)
  {
    float alpha = std::atan2((float)neighbors[k].offset.y,
                             (float)neighbors[k].offset.x);
    cell_angles.push_back(alpha);
  }

  for (int j = 1; j < z.shape.y - 1; j++)
    for (int i = 1; i < z.shape.x - 1; i++)
    {
      float sum = 0.f;

      for (size_t k = 0; k < nb; k++)
      {
        float dz = z(i, j) -
                   z(i + neighbors[k].offset.x, j + neighbors[k].offset.y);
        if (dz > 0.f)
        {
          float w = std::pow(dz * neighbors[k].inv_distance, p(i, j)) *
                    neighbors[k].weight;
          angle(i, j) += w * cell_angles[k];
          sum += w;
        }
      }

      if (sum) angle(i, j) /= sum;
    }

  return angle;
}

} // namespace hmap
