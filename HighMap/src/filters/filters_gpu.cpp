/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include <cmath>
#include <vector>

#include "cl_wrapper/run.hpp"

#include "highmap/array.hpp"
#include "highmap/curvature.hpp"
#include "highmap/filters.hpp"
#include "highmap/gradient.hpp"
#include "highmap/kernels.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/math/array.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/operator.hpp"
#include "highmap/range.hpp"

namespace hmap::gpu
{

void adaptive_relief(Array &array, float strength, float clamp_ratio)
{
  const glm::ivec2 &shape = array.shape;

  auto run = clwrapper::Run("adaptive_relief");

  run.bind_imagef("in", array.vector, shape.x, shape.y);
  run.bind_imagef("out", array.vector, shape.x, shape.y,
                  true); // out
  run.bind_arguments(shape.x, shape.y, strength, clamp_ratio);

  run.execute({shape.x, shape.y});
  run.read_imagef("out");
}

void adaptive_relief(Array       &array,
                     const Array *p_mask,
                     float        strength,
                     float        clamp_ratio)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a)
                  { gpu::adaptive_relief(a, strength, clamp_ratio); });
}

void expand(Array &array, int ir, int iterations)
{
  Array kernel = cubic_pulse({2 * ir + 1, 2 * ir + 1});
  gpu::expand(array, kernel, iterations);
}

void expand(Array &array, int ir, const Array *p_mask, int iterations)
{
  Array kernel = cubic_pulse({2 * ir + 1, 2 * ir + 1});

  if (!p_mask)
  {
    gpu::expand(array, kernel, iterations);
  }
  else
  {
    gpu::expand(array, kernel, p_mask, iterations);
  }
}

void expand(Array &array, const Array &kernel, int iterations)
{
  auto run = clwrapper::Run("expand");

  run.bind_imagef("z", array.vector, array.shape.x, array.shape.y);
  run.bind_imagef("weights", kernel.vector, kernel.shape.x, kernel.shape.y);
  run.bind_imagef("out",
                  array.vector,
                  array.shape.x,
                  array.shape.y,
                  true); // out
  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kernel.shape.x,
                     kernel.shape.y);

  for (int it = 0; it < iterations; ++it)
  {
    run.write_imagef("z");
    run.execute({array.shape.x, array.shape.y});
    run.read_imagef("out");
  }
}

void expand(Array       &array,
            const Array &kernel,
            const Array *p_mask,
            int          iterations)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { gpu::expand(a, kernel, iterations); });
}

void gamma_correction_local(Array &array, float gamma, int ir, float k)
{
  Array amin = gpu::local_min(array, ir);
  Array amax = gpu::local_max(array, ir);

  gpu::smooth_cpulse(amin, ir);
  gpu::smooth_cpulse(amax, ir);

  if (k != 0) // with smoothing
  {
    for (int j = 0; j < array.shape.y; j++)
      for (int i = 0; i < array.shape.x; i++)
      {
        float v = std::abs(array(i, j) - amin(i, j)) /
                  (amax(i, j) - amin(i, j) + 1e-30);
        v = std::sqrt(v * v + k);
        array(i, j) = std::pow(v, gamma) * (amax(i, j) - amin(i, j)) +
                      amin(i, j);
      }
  }
  else // without smoothing
  {
    for (int j = 0; j < array.shape.y; j++)
      for (int i = 0; i < array.shape.x; i++)
      {
        float v = std::abs(array(i, j) - amin(i, j)) /
                  (amax(i, j) - amin(i, j) + 1e-30);
        array(i, j) = std::pow(v, gamma) * (amax(i, j) - amin(i, j)) +
                      amin(i, j);
      }
  }
}

void gamma_correction_local(Array       &array,
                            float        gamma,
                            int          ir,
                            const Array *p_mask,
                            float        k)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a)
                  { gpu::gamma_correction_local(a, gamma, ir, k); });
}

void laplace(Array &array, float sigma, int iterations)
{
  auto run = clwrapper::Run("laplace");

  run.bind_buffer<float>("array", array.vector);
  run.bind_arguments(array.shape.x, array.shape.y, sigma);

  run.write_buffer("array");

  for (int it = 0; it < iterations; it++)
    run.execute({array.shape.x, array.shape.y});

  run.read_buffer("array");
}

void laplace(Array &array, const Array *p_mask, float sigma, int iterations)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { gpu::laplace(a, sigma, iterations); });
}

Array mean_shift(const Array &array,
                 int          ir,
                 float        talus,
                 int          iterations,
                 bool         talus_weighted)
{
  const glm::ivec2 shape = array.shape;
  Array            array_next = Array(shape);
  Array            array_prev = array;

  auto run = clwrapper::Run("mean_shift");

  run.bind_imagef("in", array_prev.vector, shape.x, shape.y);
  run.bind_imagef("out", array_next.vector, shape.x, shape.y, true);
  run.bind_arguments(shape.x, shape.y, ir, talus, talus_weighted ? 1 : 0);

  for (int it = 0; it < iterations; it++)
  {
    run.execute({shape.x, shape.y});
    run.read_imagef("out");

    if (iterations > 1)
    {
      array_prev = array_next;
      run.write_imagef("in");
    }
  }

  run.read_imagef("out");

  return array_next;
}

Array mean_shift(const Array &array,
                 int          ir,
                 float        talus,
                 const Array *p_mask,
                 int          iterations,
                 bool         talus_weighted)
{
  return transform_with_mask(
      array,
      p_mask,
      [&](const Array &a)
      { return gpu::mean_shift(a, ir, talus, iterations, talus_weighted); });
}

void median_3x3(Array &array)
{
  auto run = clwrapper::Run("median_3x3");

  run.bind_imagef("in", array.vector, array.shape.x, array.shape.y);
  run.bind_imagef("out", array.vector, array.shape.x, array.shape.y, true);
  run.bind_arguments(array.shape.x, array.shape.y);

  run.execute({array.shape.x, array.shape.y});

  run.read_imagef("out");
}

void median_3x3(Array &array, const Array *p_mask)
{
  apply_with_mask(array, p_mask, [&](Array &a) { gpu::median_3x3(a); });
}

Array median_pseudo(const Array &array, int ir)
{
  return (gpu::local_min(array, ir) + gpu::local_max(array, ir) +
          gpu::local_mean(array, ir)) /
         3.f;
}

void normal_displacement(Array &array, float amount, int ir, bool reverse)
{
  auto run = clwrapper::Run("normal_displacement");

  Array array_f = array;
  if (ir > 0) gpu::smooth_cpulse(array_f, ir);

  if (reverse) amount *= -1.f;

  run.bind_imagef("array", array.vector, array.shape.x, array.shape.y);
  run.bind_imagef("array_f", array_f.vector, array.shape.x, array.shape.y);
  run.bind_imagef("out", array.vector, array.shape.x, array.shape.y, true);
  run.bind_arguments(array.shape.x, array.shape.y, amount);

  run.execute({array.shape.x, array.shape.y});

  run.read_imagef("out");
}

void normal_displacement(Array       &array,
                         const Array *p_mask,
                         float        amount,
                         int          ir,
                         bool         reverse)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a)
                  { gpu::normal_displacement(a, amount, ir, reverse); });
}

void plateau(Array &array, const Array *p_mask, int ir, float factor)
{
  Array amin = gpu::local_min(array, ir);
  Array amax = gpu::local_max(array, ir);

  gpu::smooth_cpulse(amin, ir);
  gpu::smooth_cpulse(amax, ir);

  // last part
  auto run = clwrapper::Run("plateau_post");

  run.bind_buffer<float>("array", array.vector);
  run.bind_buffer<float>("amin", amin.vector);
  run.bind_buffer<float>("amax", amax.vector);

  helper_bind_optional_buffer(run, "mask", p_mask);

  run.bind_arguments(array.shape.x, array.shape.y, factor, p_mask ? 1 : 0);

  run.write_buffer("array");
  run.write_buffer("amin");
  run.write_buffer("amax");

  run.execute({array.shape.x, array.shape.y});

  run.read_buffer("array");
}

void plateau(Array &array, int ir, float factor)
{
  gpu::plateau(array, nullptr, ir, factor);
}

Array project_talus_along_direction(const Array &array,
                                    const Array &talus,
                                    int          direction,
                                    float        vmin)
{
  const glm::ivec2 &shape = array.shape;

  // no negative values, raises issue with atomic max in OpenCL
  const float offset = array.min();
  Array       out = array + offset;

  // D8 directions (clockwise, starting from +X)
  constexpr int d8_offsets[8][2] = {
      {1, 0},   // 0
      {1, -1},  // 1
      {0, -1},  // 2
      {-1, -1}, // 3
      {-1, 0},  // 4
      {-1, 1},  // 5
      {0, 1},   // 6
      {1, 1}    // 7
  };

  direction &= 7; // safety clamp

  int di = d8_offsets[direction][0];
  int dj = d8_offsets[direction][1];

  // apply
  auto run = clwrapper::Run("project_talus_along_direction");

  run.bind_buffer<float>("array", out.vector);
  run.bind_buffer<float>("talus", talus.vector);
  run.bind_buffer<float>("out", out.vector);

  run.bind_arguments(shape.x, shape.y, di, dj, vmin + offset);

  run.write_buffer("array");
  run.write_buffer("talus");
  run.write_buffer("out");

  run.execute({shape.x, shape.y});

  run.read_buffer("out");

  return out - offset;
}

Array project_talus_along_direction(const Array &array,
                                    const Array &talus,
                                    const Array *p_mask,
                                    int          direction,
                                    float        vmin)
{
  return transform_with_mask(
      array,
      p_mask,
      [&](const Array &a)
      { return project_talus_along_direction(a, talus, direction, vmin); });
}

void shrink(Array &array, int ir, int iterations)
{
  Array kernel = cubic_pulse({2 * ir + 1, 2 * ir + 1});
  gpu::shrink(array, kernel, iterations);
}

void shrink(Array &array, int ir, const Array *p_mask, int iterations)
{
  Array kernel = cubic_pulse({2 * ir + 1, 2 * ir + 1});

  if (!p_mask)
  {
    gpu::shrink(array, kernel, iterations);
  }
  else
  {
    gpu::shrink(array, kernel, p_mask, iterations);
  }
}

void shrink(Array &array, const Array &kernel, int iterations)
{

  auto run = clwrapper::Run("expand");

  run.bind_imagef("z", array.vector, array.shape.x, array.shape.y);
  run.bind_imagef("weights", kernel.vector, kernel.shape.x, kernel.shape.y);
  run.bind_imagef("out",
                  array.vector,
                  array.shape.x,
                  array.shape.y,
                  true); // out
  run.bind_arguments(array.shape.x,
                     array.shape.y,
                     kernel.shape.x,
                     kernel.shape.y);

  for (int it = 0; it < iterations; ++it)
  {
    float amax = array.max();
    array *= -1.f; // array <- amax - array;
    array += amax;

    run.write_imagef("z");
    run.execute({array.shape.x, array.shape.y});
    run.read_imagef("out");

    array *= -1.f; // array <- amax - array;
    array += amax;
  }
}

void shrink(Array       &array,
            const Array &kernel,
            const Array *p_mask,
            int          iterations)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { gpu::shrink(a, kernel, iterations); });
}

void smooth_cpulse(Array &array, int ir)
{
  // define kernel
  const int          nk = 2 * ir + 1;
  std::vector<float> k1d(nk);

  float sum = 0.f;
  float x0 = (float)ir;
  for (int i = 0; i < nk; i++)
  {
    float x = std::abs((float)i - x0) / (float)ir;
    k1d[i] = std::exp(-0.5f * x * x * 9.f); // σ ≈ ir/3
    sum += k1d[i];
  }

  // normalize
  for (int i = 0; i < nk; i++)
  {
    k1d[i] /= sum;
  }

  auto run = clwrapper::Run("smooth_cpulse");

  run.bind_imagef("in", array.vector, array.shape.x, array.shape.y);
  run.bind_imagef("weights", k1d, nk, 1);
  run.bind_imagef("out",
                  array.vector,
                  array.shape.x,
                  array.shape.y,
                  true); // out
  run.bind_arguments(array.shape.x, array.shape.y, ir);

  int pass_nb;

  pass_nb = 0; // x
  run.set_argument(6, pass_nb);
  run.execute({array.shape.x, array.shape.y});

  run.read_imagef("out");
  run.write_imagef("in");

  pass_nb = 1; // y
  run.set_argument(6, pass_nb);
  run.execute({array.shape.x, array.shape.y});

  run.read_imagef("out");
}

void smooth_cpulse(Array &array, int ir, const Array *p_mask)
{
  apply_with_mask(array, p_mask, [&](Array &a) { gpu::smooth_cpulse(a, ir); });
}

void smooth_cpulse_edge_removing(Array &array,
                                 float  talus,
                                 float  talus_width,
                                 int    ir)
{
  Array c = hmap::gradient_norm(array); // CPU version
  c = sigmoid(c, talus_width, 0.f /* vmin */, 1.f /* vmax */, talus);
  gpu::expand(c, ir);
  gpu::smooth_cpulse(array, ir, &c);
}

void smooth_fill(Array &array, int ir, float k, Array *p_deposition_map)
{
  Array array_bckp = array;

  gpu::smooth_cpulse(array, ir);
  array = gpu::maximum_smooth(array, array_bckp, k);

  if (p_deposition_map) *p_deposition_map = maximum(array - array_bckp, 0.f);
}

void smooth_fill(Array       &array,
                 int          ir,
                 const Array *p_mask,
                 float        k,
                 Array       *p_deposition_map)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a)
                  { gpu::smooth_fill(a, ir, k, p_deposition_map); });
}

void smooth_fill_holes(Array &array, int ir)
{
  Array array_smooth = array;
  gpu::smooth_cpulse(array_smooth, ir);

  // mask based on concave regions
  Array mask = gpu::curvature_quadric(array_smooth, 0, CurvatureType::CT_MEAN);
  clamp_min(mask, 0.f);
  make_binary(mask);

  int ic = (int)((float)ir / 2.f);
  if (ic > 1) gpu::smooth_cpulse(mask, ic);

  array = lerp(array, array_smooth, mask);
}

void smooth_fill_holes(Array &array, int ir, const Array *p_mask)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { gpu::smooth_fill_holes(a, ir); });
}

void smooth_fill_smear_peaks(Array &array, int ir)
{
  Array array_smooth = array;
  gpu::smooth_cpulse(array_smooth, ir);

  // mask based on concave regions
  Array mask = -gpu::curvature_quadric(array_smooth, 0, CurvatureType::CT_MEAN);
  clamp_min(mask, 0.f);
  make_binary(mask);

  int ic = (int)((float)ir / 2.f);
  if (ic > 0) gpu::smooth_cpulse(mask, ic);

  array = lerp(array, array_smooth, mask);
}

void smooth_fill_smear_peaks(Array &array, int ir, const Array *p_mask)
{
  apply_with_mask(array,
                  p_mask,
                  [&](Array &a) { gpu::smooth_fill_smear_peaks(a, ir); });
}

} // namespace hmap::gpu
