/* Copyright (c) 2026 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <vector>

#include "highmap/array.hpp"
#include "highmap/authoring.hpp"
#include "highmap/internal/validation.hpp"
#include "highmap/logger.hpp"

namespace hmap
{

namespace
{

// Truncated Gaussian kernel of radius ir, equal to 1 at the center and to 0 at
// distance ir (and beyond), stored on a (2 ir + 1)^2 square window
struct PushKernel
{
  int                ir = 0;
  int                n = 0;
  std::vector<float> values;

  float at(int di, int dj) const
  {
    return this->values[(dj + this->ir) * this->n + (di + this->ir)];
  }
};

PushKernel helper_make_kernel(int ir)
{
  PushKernel kernel;
  kernel.ir = ir;
  kernel.n = 2 * ir + 1;
  kernel.values.assign(static_cast<size_t>(kernel.n) * kernel.n, 0.f);

  float sigma = static_cast<float>(ir) / 3.f;
  float inv_2s2 = 1.f / (2.f * sigma * sigma);
  float r2 = static_cast<float>(ir * ir);
  float e_r = std::exp(-r2 * inv_2s2);

  for (int dj = -ir; dj <= ir; ++dj)
    for (int di = -ir; di <= ir; ++di)
    {
      float d2 = static_cast<float>(di * di + dj * dj);
      if (d2 >= r2) continue;

      kernel.values[(dj + ir) * kernel.n + (di + ir)] = (std::exp(-d2 *
                                                                  inv_2s2) -
                                                         e_r) /
                                                        (1.f - e_r);
    }

  return kernel;
}

// kernel window clipped to the array, [i0, i1[ x [j0, j1[
struct Window
{
  int i0, i1, j0, j1;
};

Window helper_window(glm::ivec2 shape, int ic, int jc, int ir)
{
  return {std::max(0, ic - ir),
          std::min(shape.x, ic + ir + 1),
          std::max(0, jc - ir),
          std::min(shape.y, jc + ir + 1)};
}

void helper_apply_push(Array            &z,
                       const PushKernel &kernel,
                       int               ic,
                       int               jc,
                       float             amplitude)
{
  Window w = helper_window(z.shape, ic, jc, kernel.ir);

  for (int j = w.j0; j < w.j1; ++j)
    for (int i = w.i0; i < w.i1; ++i)
    {
      float g = kernel.at(i - ic, j - jc);
      if (g > 0.f) z(i, j) += amplitude * g;
    }
}

// lightweight view of a constraint
struct Term
{
  const float              *target;
  const float              *weight;
  DeformationConstraintType type;
  float                     scale;
};

bool helper_is_active(DeformationConstraintType type, float d)
{
  switch (type)
  {
  case ABOVE:
    return d < 0.f;
  case BELOW:
    return d > 0.f;
  case MATCH:
  default:
    return true;
  }
}

float helper_penalty(const std::vector<Term> &terms, size_t k, float z)
{
  float p = 0.f;
  for (const Term &t : terms)
  {
    float w = t.weight[k];
    if (w == 0.f) continue;

    float d = z - t.target[k];
    if (helper_is_active(t.type, d)) p += t.scale * w * d * d;
  }
  return p;
}

struct Candidate
{
  int   i = 0;
  int   j = 0;
  int   ik = 0; // kernel index
  float amplitude = 0.f;
  float delta = 0.f; // fitness change, negative means improvement
};

// Evaluate the best push of a given kernel at (ic, jc). The amplitude is the
// least-squares optimum of the local quadratic model of the penalty (exact
// for MATCH terms, using the currently active vertices for one-sided terms),
// clamped by the slope limit; the fitness change is then computed exactly.
bool helper_evaluate(const Array              &z,
                     const std::vector<float> &penalty,
                     const std::vector<Term>  &terms,
                     const PushKernel         &kernel,
                     int                       ic,
                     int                       jc,
                     float                     amplitude_max,
                     Candidate                &out)
{
  Window w = helper_window(z.shape, ic, jc, kernel.ir);

  float s1 = 0.f;
  float s2 = 0.f;

  for (int j = w.j0; j < w.j1; ++j)
    for (int i = w.i0; i < w.i1; ++i)
    {
      float g = kernel.at(i - ic, j - jc);
      if (g <= 0.f) continue;

      size_t k = static_cast<size_t>(j) * z.shape.x + i;
      float  zv = z.vector[k];

      for (const Term &t : terms)
      {
        float wt = t.weight[k];
        if (wt == 0.f) continue;

        float d = zv - t.target[k];
        if (!helper_is_active(t.type, d)) continue;

        float sw = t.scale * wt;
        s1 += sw * g * d;
        s2 += sw * g * g;
      }
    }

  if (s2 <= 0.f) return false;

  float amplitude = std::clamp(-s1 / s2, -amplitude_max, amplitude_max);
  if (std::abs(amplitude) <= std::numeric_limits<float>::epsilon())
    return false;

  float delta = 0.f;

  for (int j = w.j0; j < w.j1; ++j)
    for (int i = w.i0; i < w.i1; ++i)
    {
      float g = kernel.at(i - ic, j - jc);
      if (g <= 0.f) continue;

      size_t k = static_cast<size_t>(j) * z.shape.x + i;
      delta += helper_penalty(terms, k, z.vector[k] + amplitude * g) -
               penalty[k];
    }

  if (delta >= 0.f) return false;

  out.i = ic;
  out.j = jc;
  out.amplitude = amplitude;
  out.delta = delta;
  return true;
}

} // namespace

void apply_gaussian_pushes(Array                           &array,
                           const std::vector<GaussianPush> &pushes)
{
  if (!validate_non_empty(array)) return;

  std::map<int, PushKernel> kernels;

  for (const GaussianPush &p : pushes)
  {
    if (p.ir < 1) continue;

    auto it = kernels.find(p.ir);
    if (it == kernels.end())
      it = kernels.emplace(p.ir, helper_make_kernel(p.ir)).first;

    helper_apply_push(array, it->second, p.i, p.j, p.amplitude);
  }
}

Array sls_deformation(const Array                              &z,
                      const std::vector<DeformationConstraint> &constraints,
                      std::uint32_t                             seed,
                      int                                       iterations,
                      int                                       ir_min,
                      int                                       ir_max,
                      int                                       n_radii,
                      float                                     talus_max,
                      int                                       n_candidates,
                      float                                     p_best,
                      float                                     top_fraction,
                      float                                     tolerance,
                      std::vector<GaussianPush>                *p_pushes)
{
  Array z_out = z;

  if (p_pushes) p_pushes->clear();

  if (!validate_non_empty(z)) return z_out;
  if (constraints.empty() || iterations <= 0) return z_out;

  for (const DeformationConstraint &c : constraints)
  {
    if (!validate_same_shape(z, c.target)) return z_out;
    if (!validate_same_shape(z, c.weight)) return z_out;
  }

  // --- parameters

  ir_min = std::max(1, ir_min);
  ir_max = std::max(ir_min, ir_max);
  n_radii = std::max(1, n_radii);
  n_candidates = std::max(1, n_candidates);
  p_best = std::clamp(p_best, 0.f, 1.f);
  top_fraction = std::clamp(top_fraction, 0.f, 1.f);

  std::vector<Term> terms;
  terms.reserve(constraints.size());
  for (const DeformationConstraint &c : constraints)
    terms.push_back(
        {c.target.vector.data(), c.weight.vector.data(), c.type, c.scale});

  // geometrically spaced radii
  std::vector<PushKernel> kernels;
  for (int k = 0; k < n_radii; ++k)
  {
    float t = n_radii > 1 ? static_cast<float>(k) / (n_radii - 1) : 0.f;
    int   ir = static_cast<int>(
        std::round(ir_min * std::pow(static_cast<float>(ir_max) / ir_min, t)));
    if (kernels.empty() || kernels.back().ir != ir)
      kernels.push_back(helper_make_kernel(ir));
  }

  // default slope limit: cross the whole elevation range with the largest
  // radius
  if (talus_max <= 0.f)
  {
    float vmin = z.min();
    float vmax = z.max();
    for (const DeformationConstraint &c : constraints)
    {
      vmin = std::min(vmin, c.target.min());
      vmax = std::max(vmax, c.target.max());
    }
    talus_max = (vmax - vmin) / static_cast<float>(ir_max);
    if (talus_max <= 0.f) return z_out;
  }

  // --- penalty map and initial fitness

  size_t             n = z.vector.size();
  std::vector<float> penalty(n);

  for (size_t k = 0; k < n; ++k)
    penalty[k] = helper_penalty(terms, k, z_out.vector[k]);

  double fitness = std::accumulate(penalty.begin(), penalty.end(), 0.0);
  if (fitness <= 0.0) return z_out;

  double fitness_stop = static_cast<double>(tolerance) * fitness;

  // --- search

  std::mt19937                          gen(seed);
  std::uniform_real_distribution<float> dis(0.f, 1.f);

  int n_uniform = n_candidates / 2;
  int n_adaptive = n_candidates - n_uniform;
  int nx = std::max(
      1,
      static_cast<int>(std::round(
          std::sqrt(static_cast<float>(n_uniform) * z.shape.x / z.shape.y))));
  int ny = std::max(1, n_uniform / nx);

  const int stall_max = 20;
  int       stall = 0;

  std::vector<glm::ivec2> cand_ij;
  std::vector<Candidate>  cand;
  std::vector<char>       cand_valid;
  std::vector<int>        improving;

  for (int it = 0; it < iterations; ++it)
  {
    if (fitness <= fitness_stop) break;

    // candidate vertices: jittered uniform grid + penalty-weighted samples
    cand_ij.clear();

    if (n_uniform > 0)
      for (int gj = 0; gj < ny; ++gj)
        for (int gi = 0; gi < nx; ++gi)
        {
          int i = std::min(z.shape.x - 1,
                           static_cast<int>((gi + dis(gen)) * z.shape.x / nx));
          int j = std::min(z.shape.y - 1,
                           static_cast<int>((gj + dis(gen)) * z.shape.y / ny));
          cand_ij.push_back({i, j});
        }

    if (n_adaptive > 0)
    {
      std::discrete_distribution<int> dd(penalty.begin(), penalty.end());
      for (int s = 0; s < n_adaptive; ++s)
      {
        int k = dd(gen);
        cand_ij.push_back({k % z.shape.x, k / z.shape.x});
      }
    }

    // evaluate the best push per candidate over all radii
    int nc = static_cast<int>(cand_ij.size());
    cand.assign(nc, Candidate());
    cand_valid.assign(nc, 0);

#pragma omp parallel for schedule(dynamic)
    for (int c = 0; c < nc; ++c)
    {
      Candidate best;
      bool      found = false;

      for (int ik = 0; ik < static_cast<int>(kernels.size()); ++ik)
      {
        Candidate tmp;
        float     amplitude_max = talus_max * kernels[ik].ir;

        if (helper_evaluate(z_out,
                            penalty,
                            terms,
                            kernels[ik],
                            cand_ij[c].x,
                            cand_ij[c].y,
                            amplitude_max,
                            tmp) &&
            tmp.delta < best.delta)
        {
          tmp.ik = ik;
          best = tmp;
          found = true;
        }
      }

      cand[c] = best;
      cand_valid[c] = found ? 1 : 0;
    }

    improving.clear();
    for (int c = 0; c < nc; ++c)
      if (cand_valid[c]) improving.push_back(c);

    if (improving.empty())
    {
      if (++stall >= stall_max) break;
      continue;
    }
    stall = 0;

    // stochastic local search: best choice with probability p_best,
    // otherwise one of the top_fraction best choices
    std::sort(improving.begin(),
              improving.end(),
              [&cand](int a, int b) { return cand[a].delta < cand[b].delta; });

    int pick = 0;
    if (dis(gen) > p_best)
    {
      int n_top = std::max(
          1,
          static_cast<int>(std::ceil(top_fraction * improving.size())));
      pick = std::min(n_top - 1, static_cast<int>(dis(gen) * n_top));
    }

    const Candidate  &sel = cand[improving[pick]];
    const PushKernel &kernel = kernels[sel.ik];

    helper_apply_push(z_out, kernel, sel.i, sel.j, sel.amplitude);

    Window w = helper_window(z.shape, sel.i, sel.j, kernel.ir);
    for (int j = w.j0; j < w.j1; ++j)
      for (int i = w.i0; i < w.i1; ++i)
      {
        size_t k = static_cast<size_t>(j) * z.shape.x + i;
        penalty[k] = helper_penalty(terms, k, z_out.vector[k]);
      }

    fitness = std::accumulate(penalty.begin(), penalty.end(), 0.0);

    if (p_pushes) p_pushes->push_back({sel.i, sel.j, kernel.ir, sel.amplitude});
  }

  return z_out;
}

} // namespace hmap
