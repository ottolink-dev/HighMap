#include <algorithm>
#include <cstdint>

#include "highmap.hpp"
#include "highmap/internal/assert.hpp"
#include "highmap/internal/timer.hpp"

#include <bit>
#include <gtest/gtest.h>

using namespace hmap;

namespace
{
bool only_finite_values(const Array &array)
{
  // These tests also run in release builds with fast-math enabled.
  return std::all_of(array.vector.begin(),
                     array.vector.end(),
                     [](float value)
                     {
                       return (std::bit_cast<std::uint32_t>(value) &
                               0x7f800000u) != 0x7f800000u;
                     });
}
} // namespace

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

    EXPECT_TRUE(only_finite_values(zc)) << "CPU NoiseType " << int(type);
    EXPECT_TRUE(only_finite_values(zg)) << "GPU NoiseType " << int(type);
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

    EXPECT_TRUE(only_finite_values(zc)) << "CPU FBM NoiseType " << int(type);
    EXPECT_TRUE(only_finite_values(zg)) << "GPU FBM NoiseType " << int(type);
    EXPECT_TRUE(assert_almost_equal(zc, zg, tol))
        << "Mismatch for FBM NoiseType " << static_cast<int>(type);
  }
}

TEST(GpuCpu, CurvatureQuadric)
{
  const float tol = 1e-4f;
  Array       z = helper_generate_array();

  std::vector<CurvatureType> types = {
      CurvatureType::CT_MIN,
      CurvatureType::CT_MAX,
      CurvatureType::CT_MEAN,
      CurvatureType::CT_GAUSSIAN,
      CurvatureType::CT_PROFILE,
      CurvatureType::CT_CONTOUR,
      CurvatureType::CT_TANGENTIAL,
      CurvatureType::CT_ACCUMULATION,
      CurvatureType::CT_SHAPE_INDEX,
      CurvatureType::CT_UNSPHERICITY,
      CurvatureType::CT_RING,
      CurvatureType::CT_ROTOR,
  };

  for (int ir : {0, 2})
  {
    for (const auto type : types)
    {
      Array zc = curvature_quadric(z, ir, type);
      Array zg = gpu::curvature_quadric(z, ir, type, false);

      EXPECT_TRUE(assert_almost_equal(zc, zg, tol))
          << "Mismatch for CurvatureType " << static_cast<int>(type)
          << " with ir = " << ir;
    }
  }
}

TEST(GpuCpu, NoiseSlope)
{
  const glm::ivec2 shape = {64, 64};
  const glm::vec2  kw = {4.f, 4.f};
  const uint       seed = 42;
  const int        octaves = 6;
  const float      weight = 0.7f;
  const float      persistence = 0.5f;
  const float      lacunarity = 2.f;
  const float      gradient_scale = 0.5f;
  const float      slope_power = 2.f;
  const float      octave_blend = 0.5f;
  const float      tol = 1e-4f;

  // isotropic
  {
    Array zc = noise_slope(NoiseType::PERLIN,
                           shape,
                           kw,
                           seed,
                           octaves,
                           weight,
                           persistence,
                           lacunarity,
                           gradient_scale,
                           slope_power,
                           octave_blend);
    Array zg = gpu::noise_slope(NoiseType::PERLIN,
                                shape,
                                kw,
                                seed,
                                octaves,
                                weight,
                                persistence,
                                lacunarity,
                                gradient_scale,
                                slope_power,
                                octave_blend);

    EXPECT_TRUE(assert_almost_equal(zc, zg, tol));
  }

  // directional bias + stretching
  {
    const float bias_angle = 45.f;
    const float bias_strength = 0.8f;
    const float angle = 30.f;
    const float stretch = 2.f;
    const float stretch_decay = 0.2f;

    Array zc = noise_slope(NoiseType::PERLIN,
                           shape,
                           kw,
                           seed,
                           octaves,
                           weight,
                           persistence,
                           lacunarity,
                           gradient_scale,
                           slope_power,
                           octave_blend,
                           bias_angle,
                           bias_strength,
                           nullptr,
                           nullptr,
                           nullptr,
                           {0.f, 1.f, 0.f, 1.f},
                           angle,
                           stretch,
                           stretch_decay);
    Array zg = gpu::noise_slope(NoiseType::PERLIN,
                                shape,
                                kw,
                                seed,
                                octaves,
                                weight,
                                persistence,
                                lacunarity,
                                gradient_scale,
                                slope_power,
                                octave_blend,
                                bias_angle,
                                bias_strength,
                                nullptr,
                                nullptr,
                                nullptr,
                                {0.f, 1.f, 0.f, 1.f},
                                {0, 0},
                                angle,
                                stretch,
                                stretch_decay);

    EXPECT_TRUE(assert_almost_equal(zc, zg, tol));
  }
}
