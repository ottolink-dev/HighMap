#include "highmap.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/dbg/timer.hpp"

#include <gtest/gtest.h>

using namespace hmap;

Array helper_generate_array()
{
  const glm::ivec2 shape = {512, 512};
  const glm::vec2  kw = {4.f, 4.f};
  const int        seed = 0;

  Array z = noise_fbm(NoiseType::PERLIN, shape, kw, seed);
  return z;
}

TEST(GpuCpu, ClosingByReconstruction)
{
  Timer::Clear();

  const float tol = 1e0f;
  const int   ir = 64;

  Array z = helper_generate_array();

  {
    Timer::Start("CPU");
    Array zc = closing_by_reconstruction(z, ir);
    Timer::Stop("CPU");

    Timer::Start("GPU");
    Array zg = gpu::closing_by_reconstruction(z, ir);
    Timer::Stop("GPU");

    bool ret = assert_almost_equal(zc, zg, tol);
    EXPECT_EQ(ret, true);
  }

  {
    Timer::Start("CPU");
    Array zc = opening_by_reconstruction(z, ir);
    Timer::Stop("CPU");

    Timer::Start("GPU");
    Array zg = gpu::opening_by_reconstruction(z, ir);
    Timer::Stop("GPU");

    bool ret = assert_almost_equal(zc, zg, tol);
    EXPECT_EQ(ret, true);
  }

  auto map = Timer::DumpDurations();
  for (auto [k, v] : map)
    RecordProperty(k, v);
}

TEST(GpuCpu, MorphologicalOperators)
{
  Timer::Clear();

  const float tol = 1e0f;
  const int   ir = 32;

  Array z = helper_generate_array();

  std::vector<MorphologyOperation> ops = {
      MorphologyOperation::MO_BORDER,
      MorphologyOperation::MO_CLOSING,
      MorphologyOperation::MO_DILATION,
      MorphologyOperation::MO_EROSION,
      MorphologyOperation::MO_OPENING,
      MorphologyOperation::MO_GRADIENT,
      MorphologyOperation::MO_TOP_HAT,
      MorphologyOperation::MO_BLACK_HAT,
      MorphologyOperation::MO_LAPLACIAN,
      MorphologyOperation::MO_CLOSING_BY_RECONSTRUCTION,
      MorphologyOperation::MO_OPENING_BY_RECONSTRUCTION,
  };

  for (const auto op : ops)
  {
    Timer::Start("CPU");
    Array zc = morphological_operators(z, ir, op);
    Timer::Stop("CPU");

    Timer::Start("GPU");
    Array zg = gpu::morphological_operators(z, ir, op);
    Timer::Stop("GPU");

    bool ret = assert_almost_equal(zc, zg, tol);
    EXPECT_EQ(ret, true);
  }

  auto map = Timer::DumpDurations();
  for (auto [k, v] : map)
    RecordProperty(k, v);
}

TEST(GpuCpu, HarmonicInterpolation)
{
  Timer::Clear();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  Array            mask = z;
  clamp_min(mask, 0.f);

  Timer::Start("CPU");
  Array zc = harmonic_interpolation(z, mask, 500, 1e-5f);
  Timer::Stop("CPU");

  Timer::Start("GPU");
  Array zg = gpu::harmonic_interpolation(z, mask, 500, 1e-5f);
  Timer::Stop("GPU");

  // Both CPU SOR and GPU Red-Black SOR solve the same Laplace problem
  bool ret = assert_almost_equal(zc, zg, 1e-2f);
  EXPECT_EQ(ret, true);
}

TEST(GpuCpu, WaterDepthFromMask)
{
  Timer::Clear();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z);

  Array mask = select_rivers(z, 1.f / shape.x, 10.f);
  remap(mask);

  const float mask_threshold = 0.1f;
  const int   iterations_max = 1000;
  const float tolerance = 1e-3f;

  Timer::Start("CPU");
  Array wc = water_depth_from_mask(z,
                                   mask,
                                   mask_threshold,
                                   iterations_max,
                                   tolerance);
  Timer::Stop("CPU");

  Timer::Start("GPU");
  Array wg = gpu::water_depth_from_mask(z,
                                        mask,
                                        mask_threshold,
                                        iterations_max,
                                        tolerance);
  Timer::Stop("GPU");

  bool ret = assert_almost_equal(wc, wg, 1e-2f);
  EXPECT_EQ(ret, true);
}

TEST(GpuCpu, Noise)
{
  const glm::ivec2 shape = {64, 64};
  const glm::vec2  kw = {4.f, 4.f};
  const uint       seed = 42;
  const float      tol = 1e-5f;

  Array nx = 0.25f * noise(NoiseType::PERLIN, shape, kw, seed + 10);
  Array ny = 0.25f * noise(NoiseType::PERLIN, shape, kw, seed + 20);

  std::vector<NoiseType> types = {
      NoiseType::PERLIN,
      NoiseType::PERLIN_BILLOW,
      NoiseType::PERLIN_HALF,
      NoiseType::SIMPLEX2,
      NoiseType::SIMPLEX2S,
      NoiseType::VALUE,
      NoiseType::VALUE_CUBIC,
      NoiseType::VALUE_LINEAR,
      NoiseType::WORLEY,
      NoiseType::WORLEY_DOUBLE,
  };

  for (const auto type : types)
  {
    Array zc = noise(type, shape, kw, seed, &nx, &ny);
    Array zg = gpu::noise(type, shape, kw, seed, &nx, &ny);

    EXPECT_TRUE(assert_almost_equal(zc, zg, tol))
        << "Mismatch for NoiseType " << static_cast<int>(type);
  }
}

TEST(GpuCpu, NoiseFbm)
{
  const glm::ivec2 shape = {64, 64};
  const glm::vec2  kw = {4.f, 4.f};
  const uint       seed = 42;
  const int        octaves = 6;
  const float      weight = 0.7f;
  const float      persistence = 0.5f;
  const float      lacunarity = 2.f;
  const float      tol = 1e-5f;

  Array ctrl = 0.5f * (1.f + noise(NoiseType::PERLIN, shape, kw, seed + 5));
  Array nx = 0.25f * noise(NoiseType::PERLIN, shape, kw, seed + 10);
  Array ny = 0.25f * noise(NoiseType::PERLIN, shape, kw, seed + 20);

  std::vector<NoiseType> types = {
      NoiseType::PERLIN,
      NoiseType::PERLIN_BILLOW,
      NoiseType::PERLIN_HALF,
      NoiseType::SIMPLEX2,
      NoiseType::SIMPLEX2S,
      NoiseType::VALUE,
      NoiseType::VALUE_CUBIC,
      NoiseType::VALUE_LINEAR,
      NoiseType::WORLEY,
      NoiseType::WORLEY_DOUBLE,
  };

  for (const auto type : types)
  {
    Array zc = noise_fbm(type,
                         shape,
                         kw,
                         seed,
                         octaves,
                         weight,
                         persistence,
                         lacunarity,
                         nullptr,
                         &nx,
                         &ny);
    Array zg = gpu::noise_fbm(type,
                              shape,
                              kw,
                              seed,
                              octaves,
                              weight,
                              persistence,
                              lacunarity,
                              nullptr,
                              &nx,
                              &ny);

    EXPECT_TRUE(assert_almost_equal(zc, zg, tol))
        << "Mismatch for FBM NoiseType " << static_cast<int>(type);
  }
}
