/* Copyright (c) 2023 Otto Link. Distributed under the terms of the GNU General
 * Public License. The full license is in the file LICENSE, distributed with
 * this software. */
#include "highmap/array.hpp"
#include "highmap/primitives.hpp"

#include <gtest/gtest.h>

namespace hmap
{

TEST(NoiseStretch, CpuIsotropicPreserved)
{
  const glm::ivec2 shape = {64, 64};
  const glm::vec2  kw = {4.f, 4.f};
  const uint       seed = 42;

  Array z0 = noise(NoiseType::PERLIN, shape, kw, seed);
  Array z1 = noise(NoiseType::PERLIN,
                   shape,
                   kw,
                   seed,
                   nullptr,
                   nullptr,
                   {0.f, 1.f, 0.f, 1.f},
                   37.f,
                   1.f);
  EXPECT_EQ(z0.vector, z1.vector);

  Array f0 = noise_fbm(NoiseType::PERLIN, shape, kw, seed);
  Array f1 = noise_fbm(NoiseType::PERLIN,
                       shape,
                       kw,
                       seed,
                       8,
                       0.7f,
                       0.5f,
                       2.f,
                       nullptr,
                       nullptr,
                       nullptr,
                       {0.f, 1.f, 0.f, 1.f},
                       37.f,
                       1.f,
                       0.5f);
  EXPECT_EQ(f0.vector, f1.vector);
}

TEST(NoiseStretch, CpuDirectionalMatchesAxisScaling)
{
  const glm::ivec2 shape = {64, 64};
  const glm::vec2  kw = {4.f, 4.f};
  const uint       seed = 42;

  // stretch = 2 along angle = 0 deg should match halving kw.x
  Array z_ref_x = noise(NoiseType::PERLIN, shape, {2.f, 4.f}, seed);
  Array z_str_x = noise(NoiseType::PERLIN,
                        shape,
                        kw,
                        seed,
                        nullptr,
                        nullptr,
                        {0.f, 1.f, 0.f, 1.f},
                        0.f,
                        2.f);

  for (size_t i = 0; i < z_ref_x.vector.size(); ++i)
    EXPECT_NEAR(z_ref_x.vector[i], z_str_x.vector[i], 1e-5f);

  // stretch = 2 along angle = 90 deg should match halving kw.y
  Array z_ref_y = noise(NoiseType::PERLIN, shape, {4.f, 2.f}, seed);
  Array z_str_y = noise(NoiseType::PERLIN,
                        shape,
                        kw,
                        seed,
                        nullptr,
                        nullptr,
                        {0.f, 1.f, 0.f, 1.f},
                        90.f,
                        2.f);

  for (size_t i = 0; i < z_ref_y.vector.size(); ++i)
    EXPECT_NEAR(z_ref_y.vector[i], z_str_y.vector[i], 1e-5f);
}

TEST(NoiseStretch, CpuOctaveDecay)
{
  const glm::ivec2 shape = {64, 64};
  const glm::vec2  kw = {4.f, 4.f};
  const uint       seed = 42;

  Array f_no_decay = noise_fbm(NoiseType::PERLIN,
                               shape,
                               kw,
                               seed,
                               4,
                               0.7f,
                               0.5f,
                               2.f,
                               nullptr,
                               nullptr,
                               nullptr,
                               {0.f, 1.f, 0.f, 1.f},
                               45.f,
                               3.f,
                               0.f);
  Array f_with_decay = noise_fbm(NoiseType::PERLIN,
                                 shape,
                                 kw,
                                 seed,
                                 4,
                                 0.7f,
                                 0.5f,
                                 2.f,
                                 nullptr,
                                 nullptr,
                                 nullptr,
                                 {0.f, 1.f, 0.f, 1.f},
                                 45.f,
                                 3.f,
                                 0.5f);

  EXPECT_NE(f_no_decay.vector, f_with_decay.vector);

  // also verify other fractal noise variants compile and run with stretch
  Array iq = noise_iq(NoiseType::PERLIN,
                      shape,
                      kw,
                      seed,
                      4,
                      0.7f,
                      0.5f,
                      2.f,
                      0.05f,
                      nullptr,
                      nullptr,
                      nullptr,
                      {0.f, 1.f, 0.f, 1.f},
                      30.f,
                      2.f,
                      0.5f);
  EXPECT_EQ(iq.shape, shape);
}

TEST(NoiseStretch, GpuIsotropicAndDirectional)
{
  const glm::ivec2 shape = {64, 64};
  const glm::vec2  kw = {4.f, 4.f};
  const uint       seed = 42;

  Array g0 = gpu::noise(NoiseType::PERLIN, shape, kw, seed);
  Array g1 = gpu::noise(NoiseType::PERLIN,
                        shape,
                        kw,
                        seed,
                        nullptr,
                        nullptr,
                        {0.f, 1.f, 0.f, 1.f},
                        {0, 0},
                        37.f,
                        1.f);
  EXPECT_EQ(g0.vector, g1.vector);

  Array g_ref_x = gpu::noise(NoiseType::PERLIN, shape, {2.f, 4.f}, seed);
  Array g_str_x = gpu::noise(NoiseType::PERLIN,
                             shape,
                             kw,
                             seed,
                             nullptr,
                             nullptr,
                             {0.f, 1.f, 0.f, 1.f},
                             {0, 0},
                             0.f,
                             2.f);
  for (size_t i = 0; i < g_ref_x.vector.size(); ++i)
    EXPECT_NEAR(g_ref_x.vector[i], g_str_x.vector[i], 1e-5f);

  Array gf0 = gpu::noise_fbm(NoiseType::PERLIN, shape, kw, seed);
  Array gf1 = gpu::noise_fbm(NoiseType::PERLIN,
                             shape,
                             kw,
                             seed,
                             8,
                             0.7f,
                             0.5f,
                             2.f,
                             nullptr,
                             nullptr,
                             nullptr,
                             {0.f, 1.f, 0.f, 1.f},
                             {0, 0},
                             37.f,
                             1.f,
                             0.5f);
  EXPECT_EQ(gf0.vector, gf1.vector);

  Array gf_decay = gpu::noise_fbm(NoiseType::PERLIN,
                                  shape,
                                  kw,
                                  seed,
                                  4,
                                  0.7f,
                                  0.5f,
                                  2.f,
                                  nullptr,
                                  nullptr,
                                  nullptr,
                                  {0.f, 1.f, 0.f, 1.f},
                                  {0, 0},
                                  45.f,
                                  3.f,
                                  0.5f);
  EXPECT_NE(gf0.vector, gf_decay.vector);
}

TEST(NoiseStretch, AllNoiseTypesCpuAndGpu)
{
  const glm::ivec2             shape = {32, 32};
  const glm::vec2              kw = {4.f, 4.f};
  const uint                   seed = 42;
  const std::vector<NoiseType> types = {NoiseType::PERLIN,
                                        NoiseType::PERLIN_BILLOW,
                                        NoiseType::PERLIN_HALF,
                                        NoiseType::SIMPLEX2,
                                        NoiseType::SIMPLEX2S,
                                        NoiseType::VALUE,
                                        NoiseType::VALUE_CUBIC,
                                        NoiseType::VALUE_LINEAR,
                                        NoiseType::WORLEY,
                                        NoiseType::WORLEY_DOUBLE,
                                        NoiseType::WORLEY_VALUE};

  for (NoiseType nt : types)
  {
    Array c_n = noise(nt,
                      shape,
                      kw,
                      seed,
                      nullptr,
                      nullptr,
                      {0.f, 1.f, 0.f, 1.f},
                      30.f,
                      2.f);
    Array c_f = noise_fbm(nt,
                          shape,
                          kw,
                          seed,
                          4,
                          0.7f,
                          0.5f,
                          2.f,
                          nullptr,
                          nullptr,
                          nullptr,
                          {0.f, 1.f, 0.f, 1.f},
                          30.f,
                          2.f,
                          0.5f);
    Array g_n = gpu::noise(nt,
                           shape,
                           kw,
                           seed,
                           nullptr,
                           nullptr,
                           {0.f, 1.f, 0.f, 1.f},
                           {0, 0},
                           30.f,
                           2.f);
    Array g_f = gpu::noise_fbm(nt,
                               shape,
                               kw,
                               seed,
                               4,
                               0.7f,
                               0.5f,
                               2.f,
                               nullptr,
                               nullptr,
                               nullptr,
                               {0.f, 1.f, 0.f, 1.f},
                               {0, 0},
                               30.f,
                               2.f,
                               0.5f);

    EXPECT_EQ(c_n.shape, shape);
    EXPECT_EQ(c_f.shape, shape);
    EXPECT_EQ(g_n.shape, shape);
    EXPECT_EQ(g_f.shape, shape);
  }
}

} // namespace hmap
