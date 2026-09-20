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

TEST_F(JaggedTest, AmpZeroReturnsInput)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);

  // amp = 0.0 means no bubble dome elevation added
  Array out_a0 = gpu::jagged(input,
                             glm::vec2(4.f, 4.f),
                             0.f,
                             1,
                             {0.5f, 0.5f},
                             1.f,
                             0.f);

  EXPECT_EQ(out_a0.shape, shape);
  EXPECT_TRUE(assert_almost_equal(out_a0, input, 1e-5f));
}

TEST_F(JaggedTest, AmpNonZeroModifies)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);

  Array out = gpu::jagged(input, glm::vec2(4.f, 4.f), 0.2f);

  EXPECT_FALSE(assert_almost_equal(out, input, 1e-3f));
}

TEST_F(JaggedTest, GammaModulatesDomeShape)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);

  Array out_g05 = gpu::jagged(input, 6.f, 0.2f, 1, {0.5f, 0.5f}, 0.5f);
  Array out_g10 = gpu::jagged(input, 6.f, 0.2f, 1, {0.5f, 0.5f}, 1.0f);
  Array out_g20 = gpu::jagged(input, 6.f, 0.2f, 1, {0.5f, 0.5f}, 2.0f);

  EXPECT_FALSE(assert_almost_equal(out_g05, out_g10, 1e-3f));
  EXPECT_FALSE(assert_almost_equal(out_g10, out_g20, 1e-3f));
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
                          0.2f,
                          1,
                          {0.5f, 0.5f},
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

  Array out1 = gpu::jagged(input, 5.f, 0.15f, 7, {0.6f, 0.6f}, 1.5f, 30.f);
  Array out2 = gpu::jagged(input,
                           glm::vec2(5.f, 5.f),
                           0.15f,
                           7,
                           {0.6f, 0.6f},
                           1.5f,
                           30.f);

  EXPECT_TRUE(assert_almost_equal(out1, out2, 1e-5f));
}

TEST_F(JaggedTest, AngleRotatesPattern)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);

  Array out_a0 = gpu::jagged(input, 6.f, 0.2f, 1, {0.5f, 0.5f}, 1.f, 0.f);
  Array out_a45 = gpu::jagged(input, 6.f, 0.2f, 1, {0.5f, 0.5f}, 1.f, 45.f);

  EXPECT_FALSE(assert_almost_equal(out_a0, out_a45, 1e-3f));
}

TEST_F(JaggedTest, EmptyArrayReturnsEmpty)
{
  Array input;
  Array out = gpu::jagged(input, 4.f);
  EXPECT_TRUE(out.vector.empty());
}
