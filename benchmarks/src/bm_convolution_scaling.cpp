#include "highmap/array.hpp"
#include "highmap/convolve.hpp"
#include "highmap/filters.hpp"
#include "highmap/kernels.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/primitives.hpp"

#include <benchmark/benchmark.h>

using namespace hmap;

// ------------------------------------------------------------
// Args helper
// ------------------------------------------------------------

static void conv_scaling_args(benchmark::internal::Benchmark *b)
{
  std::vector<int> sizes = {128, 256, 512};
  std::vector<int> radii = {2, 8, 16, 32};

  for (int n : sizes)
    for (int r : radii)
    {
      if (r < n / 2) b->Args({n, r});
    }
}

// ------------------------------------------------------------
// Full 2D Convolution (CPU)
// ------------------------------------------------------------

static void BM_convolve2d_CPU(benchmark::State &state)
{
  const int n = state.range(0);
  const int r = state.range(1);
  const int ksize = 2 * r + 1;

  Array input = white(glm::vec2(n, n), 0.f, 1.f, 42);
  Array kernel = cubic_pulse(glm::ivec2(ksize, ksize));

  for (auto _ : state)
  {
    Array out = convolve2d(input, kernel);
    benchmark::DoNotOptimize(out);
  }

  state.SetItemsProcessed(int64_t(state.iterations()) * n * n);
}

// ------------------------------------------------------------
// SVD 2D Convolution (CPU, rank 3)
// ------------------------------------------------------------

static void BM_convolve2d_svd_CPU(benchmark::State &state)
{
  const int n = state.range(0);
  const int r = state.range(1);
  const int ksize = 2 * r + 1;

  Array input = white(glm::vec2(n, n), 0.f, 1.f, 42);
  Array kernel = cubic_pulse(glm::ivec2(ksize, ksize));

  for (auto _ : state)
  {
    Array out = convolve2d_svd(input, kernel, 3);
    benchmark::DoNotOptimize(out);
  }

  state.SetItemsProcessed(int64_t(state.iterations()) * n * n);
}

// ------------------------------------------------------------
// Bilateral Filter (OpenCL 2D brute-force)
// ------------------------------------------------------------

static void BM_bilateral_filter_GPU(benchmark::State &state)
{
  const int n = state.range(0);
  const int r = state.range(1);

  Array input = white(glm::vec2(n, n), 0.f, 1.f, 42);

  for (auto _ : state)
  {
    Array out = gpu::bilateral_filter(input, r, 1.0f);
    benchmark::DoNotOptimize(out);
  }

  state.SetItemsProcessed(int64_t(state.iterations()) * n * n);
}

// ------------------------------------------------------------
// Local Max (OpenCL 2D disk window)
// ------------------------------------------------------------

static void BM_local_max_GPU(benchmark::State &state)
{
  const int n = state.range(0);
  const int r = state.range(1);

  Array input = white(glm::vec2(n, n), 0.f, 1.f, 42);

  for (auto _ : state)
  {
    Array out = gpu::local_max(input, r);
    benchmark::DoNotOptimize(out);
  }

  state.SetItemsProcessed(int64_t(state.iterations()) * n * n);
}

static void BM_local_max_octagon_GPU(benchmark::State &state)
{
  const int n = state.range(0);
  const int r = state.range(1);

  Array input = white(glm::vec2(n, n), 0.f, 1.f, 42);

  for (auto _ : state)
  {
    Array out = gpu::local_max_octagon(input, r);
    benchmark::DoNotOptimize(out);
  }

  state.SetItemsProcessed(int64_t(state.iterations()) * n * n);
}

// ------------------------------------------------------------
// Local Max Square (OpenCL 2-pass separable)
// ------------------------------------------------------------

static void BM_local_max_square_GPU(benchmark::State &state)
{
  const int n = state.range(0);
  const int r = state.range(1);

  Array input = white(glm::vec2(n, n), 0.f, 1.f, 42);

  for (auto _ : state)
  {
    Array out = gpu::local_max_square(input, r);
    benchmark::DoNotOptimize(out);
  }

  state.SetItemsProcessed(int64_t(state.iterations()) * n * n);
}

// ------------------------------------------------------------
// Local Max (CPU 1D monotonic deque separable)
// ------------------------------------------------------------

static void BM_local_max_CPU(benchmark::State &state)
{
  const int n = state.range(0);
  const int r = state.range(1);

  Array input = white(glm::vec2(n, n), 0.f, 1.f, 42);

  for (auto _ : state)
  {
    Array out = local_max(input, r);
    benchmark::DoNotOptimize(out);
  }

  state.SetItemsProcessed(int64_t(state.iterations()) * n * n);
}

// ------------------------------------------------------------
// Ruggedness (OpenCL 2D disk window)
// ------------------------------------------------------------

static void BM_ruggedness_GPU(benchmark::State &state)
{
  const int n = state.range(0);
  const int r = state.range(1);

  Array input = white(glm::vec2(n, n), 0.f, 1.f, 42);

  for (auto _ : state)
  {
    Array out = gpu::ruggedness(input, r);
    benchmark::DoNotOptimize(out);
  }

  state.SetItemsProcessed(int64_t(state.iterations()) * n * n);
}

// ------------------------------------------------------------
// Sparse Max Convolution (OpenCL scatter/atomic)
// ------------------------------------------------------------

static void BM_sparse_max_conv_GPU(benchmark::State &state)
{
  const int n = state.range(0);
  const int r = state.range(1);
  const int ksize = 2 * r + 1;

  Array input = white_sparse(glm::ivec2(n, n), 0.f, 1.f, 0.05f, 42);
  Array kernel = cubic_pulse(glm::ivec2(ksize, ksize));

  for (auto _ : state)
  {
    Array out = gpu::sparse_max_convolution(input, kernel, false);
    benchmark::DoNotOptimize(out);
  }

  state.SetItemsProcessed(int64_t(state.iterations()) * n * n);
}

// Only run CPU 2D convolve for moderate sizes to avoid timeout
static void conv2d_cpu_args(benchmark::internal::Benchmark *b)
{
  b->Args({128, 2});
  b->Args({128, 8});
  b->Args({128, 16});
  b->Args({256, 2});
  b->Args({256, 8});
}

BENCHMARK(BM_convolve2d_CPU)->Apply(conv2d_cpu_args);
BENCHMARK(BM_convolve2d_svd_CPU)->Apply(conv_scaling_args);
BENCHMARK(BM_bilateral_filter_GPU)->Apply(conv_scaling_args);
BENCHMARK(BM_local_max_GPU)->Apply(conv_scaling_args);
BENCHMARK(BM_local_max_octagon_GPU)->Apply(conv_scaling_args);
BENCHMARK(BM_local_max_square_GPU)->Apply(conv_scaling_args);
BENCHMARK(BM_local_max_CPU)->Apply(conv_scaling_args);
BENCHMARK(BM_ruggedness_GPU)->Apply(conv_scaling_args);
BENCHMARK(BM_sparse_max_conv_GPU)->Apply(conv_scaling_args);
