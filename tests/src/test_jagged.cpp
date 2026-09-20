#include "highmap/array.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/filters.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/primitives.hpp"
#include "highmap/range.hpp"

#include <gtest/gtest.h>

using namespace hmap;

class JaggedTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    gpu::init_opencl();
  }
};

TEST_F(JaggedTest, ConstantFieldPreserved)
{
  Array input = Array(glm::ivec2(64, 64), 3.5f);

  Array out = gpu::jagged(input,
                          glm::vec2(4.f, 4.f),
                          1,
                          {0.5f, 0.5f},
                          1.f,
                          1.f);

  EXPECT_TRUE(assert_almost_equal(out, input, 1e-4f));
}

TEST_F(JaggedTest, GammaOneReturnsInput)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(input, 0.1f, 0.9f);

  // gamma = 1.0 means exponent = 1.0 (no effect)
  Array out_g1 = gpu::jagged(input,
                             glm::vec2(4.f, 4.f),
                             1,
                             {0.5f, 0.5f},
                             1.f,
                             1.f);

  EXPECT_EQ(out_g1.shape, shape);
  EXPECT_TRUE(assert_almost_equal(out_g1, input, 1e-4f));
}

TEST_F(JaggedTest, FactorZeroReturnsInput)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(input, 0.1f, 0.9f);

  // factor = 0.0 disables jagged effect
  Array out_f0 = gpu::jagged(input,
                             glm::vec2(4.f, 4.f),
                             1,
                             {0.5f, 0.5f},
                             0.5f,
                             1.f,
                             0.f);

  EXPECT_EQ(out_f0.shape, shape);
  EXPECT_TRUE(assert_almost_equal(out_f0, input, 1e-4f));
}

TEST_F(JaggedTest, ShapeGammaModulatesEffect)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(input, 0.1f, 0.9f);

  // Compare unmodulated shape_gamma = 0 vs modulated shape_gamma = 1.0 vs
  // shape_gamma = 2.0 with gamma = 2.0 (quadratic compression towards 0)
  Array out_sg0 = gpu::jagged(input, 6.f, 1, {0.5f, 0.5f}, 2.f, 0.f);
  Array out_sg1 = gpu::jagged(input, 6.f, 1, {0.5f, 0.5f}, 2.f, 1.f);
  Array out_sg2 = gpu::jagged(input, 6.f, 1, {0.5f, 0.5f}, 2.f, 2.f);

  EXPECT_FALSE(assert_almost_equal(out_sg0, out_sg1, 1e-3f));
  EXPECT_FALSE(assert_almost_equal(out_sg1, out_sg2, 1e-3f));
}

TEST_F(JaggedTest, MaskPreservesUnmaskedArea)
{
  const int  nx = 64;
  const int  ny = 64;
  glm::ivec2 shape = {nx, ny};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 10);
  Array      mask(shape, 0.f);

  // Mask top half only
  for (int j = ny / 2; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      mask(i, j) = 1.f;

  Array out = gpu::jagged(input,
                          glm::vec2(6.f, 6.f),
                          1,
                          {0.5f, 0.5f},
                          1.5f,
                          1.f,
                          1.f,
                          0.f,
                          &mask);

  // Bottom half (mask = 0) should remain unchanged
  for (int j = 0; j < ny / 2; ++j)
    for (int i = 0; i < nx; ++i)
      EXPECT_NEAR(out(i, j), input(i, j), 1e-5f);
}

TEST_F(JaggedTest, ScalarKwOverloadEquivalence)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 123);
  remap(input, 0.1f, 0.9f);

  Array out1 = gpu::jagged(input, 5.f, 7, {0.6f, 0.6f}, 1.5f, 1.5f);
  Array out2 = gpu::jagged(input,
                           glm::vec2(5.f, 5.f),
                           7,
                           {0.6f, 0.6f},
                           1.5f,
                           1.5f);

  EXPECT_TRUE(assert_almost_equal(out1, out2, 1e-5f));
}

TEST_F(JaggedTest, AngleRotatesPattern)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(input, 0.1f, 0.9f);

  Array out_a0 = gpu::jagged(input, 6.f, 1, {0.5f, 0.5f}, 1.5f, 1.f, 1.f, 0.f);
  Array out_a45 = gpu::jagged(input, 6.f, 1, {0.5f, 0.5f}, 1.5f, 1.f, 1.f, 45.f);

  EXPECT_FALSE(assert_almost_equal(out_a0, out_a45, 1e-3f));
}

TEST_F(JaggedTest, EmptyArrayReturnsEmpty)
{
  Array input;
  Array out = gpu::jagged(input, 4.f);
  EXPECT_TRUE(out.vector.empty());
}
