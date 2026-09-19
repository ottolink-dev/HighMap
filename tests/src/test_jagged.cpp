#include "highmap/array.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/filters.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/primitives.hpp"

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

TEST_F(JaggedTest, FactorZeroReturnsVoronoiCells)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);

  // shape_factor = 0 disables edge falloff, returning pure Voronoi cell
  // constants when factor = 0
  Array out_f0 = gpu::jagged(input,
                             glm::vec2(4.f, 4.f),
                             1,
                             {0.5f, 0.5f},
                             0.f,
                             0.f);
  Array out_f1 = gpu::jagged(input,
                             glm::vec2(4.f, 4.f),
                             1,
                             {0.5f, 0.5f},
                             1.f,
                             0.f);

  EXPECT_EQ(out_f0.shape, shape);
  EXPECT_EQ(out_f1.shape, shape);
  EXPECT_FALSE(assert_almost_equal(out_f0, out_f1, 1e-3f));
}

TEST_F(JaggedTest, ShapeFactorModulatesEffect)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);

  // Compare unmodulated shape_factor = 0 vs modulated shape_factor = 1.0 vs
  // shape_factor = 2.0
  Array out_sf0 = gpu::jagged(input, 6.f, 1, {0.5f, 0.5f}, 1.f, 0.f);
  Array out_sf1 = gpu::jagged(input, 6.f, 1, {0.5f, 0.5f}, 1.f, 1.f);
  Array out_sf2 = gpu::jagged(input, 6.f, 1, {0.5f, 0.5f}, 1.f, 2.f);

  EXPECT_FALSE(assert_almost_equal(out_sf0, out_sf1, 1e-3f));
  EXPECT_FALSE(assert_almost_equal(out_sf1, out_sf2, 1e-3f));
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

  Array out =
      gpu::jagged(input, glm::vec2(6.f, 6.f), 1, {0.5f, 0.5f}, 1.f, 1.f, &mask);

  // Bottom half (mask = 0) should remain unchanged
  for (int j = 0; j < ny / 2; ++j)
    for (int i = 0; i < nx; ++i)
      EXPECT_NEAR(out(i, j), input(i, j), 1e-5f);
}

TEST_F(JaggedTest, ScalarKwOverloadEquivalence)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 123);

  Array out1 = gpu::jagged(input, 5.f, 7, {0.6f, 0.6f}, 1.f, 1.5f);
  Array out2 = gpu::jagged(input,
                           glm::vec2(5.f, 5.f),
                           7,
                           {0.6f, 0.6f},
                           1.f,
                           1.5f);

  EXPECT_TRUE(assert_almost_equal(out1, out2, 1e-5f));
}

TEST_F(JaggedTest, EmptyArrayReturnsEmpty)
{
  Array input;
  Array out = gpu::jagged(input, 4.f);
  EXPECT_TRUE(out.vector.empty());
}
