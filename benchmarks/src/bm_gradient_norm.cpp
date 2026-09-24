#include "highmap/array.hpp"
#include "highmap/gradient.hpp"
#include "highmap/primitives.hpp"

#include <benchmark/benchmark.h>

using namespace hmap;

static void BM_gradient_norm_CPU(benchmark::State &state)
{
  const int n = state.range(0); // mesh size

  Array input = white(glm::vec2(n, n), 0.f, 1.f, 42);

  for (auto _ : state)
  {
    Array out = input;
    hmap::gradient_norm(out);
    benchmark::DoNotOptimize(out);
  }

  state.SetItemsProcessed(int64_t(state.iterations()) * n * n);
}

static void gradient_norm_args(benchmark::internal::Benchmark *b)
{
  std::vector<int> sizes = {128, 256, 512, 1024};

  for (int n : sizes)
    b->Args({n});
}

BENCHMARK(BM_gradient_norm_CPU)->Apply(gradient_norm_args);
