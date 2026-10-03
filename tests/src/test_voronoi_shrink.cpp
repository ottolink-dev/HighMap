#include "highmap/array.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/filters.hpp"
#include "highmap/opencl/gpu_opencl.hpp"
#include "highmap/primitives.hpp"

#include <gtest/gtest.h>

using namespace hmap;

class VoronoiShrinkTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    gpu::init_opencl();
  }
};

TEST_F(VoronoiShrinkTest, ShrinkFactorOnePreservesInput)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);

  // shrink_factor = 1.0 means no shrinking / no rifts
  Array out = gpu::voronoi_shrink(input,
                                  glm::vec2(4.f, 4.f),
                                  1.0f,
                                  -999.f,
                                  1,
                                  {0.5f, 0.5f},
                                  0.f);

  EXPECT_EQ(out.shape, shape);
  EXPECT_TRUE(assert_almost_equal(out, input, 1e-4f));
}

TEST_F(VoronoiShrinkTest, ShrinkingCreatesRiftsWithFillValue)
{
  glm::ivec2 shape = {128, 128};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  float      fill_val = -10.f;

  Array out = gpu::voronoi_shrink(input,
                                  glm::vec2(4.f, 4.f),
                                  0.5f,
                                  fill_val,
                                  1,
                                  {0.5f, 0.5f});

  // Verify that some pixels are filled with fill_val (rifts)
  int fill_count = 0;
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      if (std::abs(out(i, j) - fill_val) < 1e-4f)
      {
        fill_count++;
      }
    }

  EXPECT_GT(fill_count, 0);
  EXPECT_LT(fill_count, shape.x * shape.y);
}

TEST_F(VoronoiShrinkTest, ScalarKwOverloadEquivalence)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 123);

  Array out1 =
      gpu::voronoi_shrink(input, 5.f, 0.7f, 0.f, 7, {0.6f, 0.6f}, 30.f);
  Array out2 = gpu::voronoi_shrink(input,
                                   glm::vec2(5.f, 5.f),
                                   0.7f,
                                   0.f,
                                   7,
                                   {0.6f, 0.6f},
                                   30.f);

  EXPECT_TRUE(assert_almost_equal(out1, out2, 1e-5f));
}

TEST_F(VoronoiShrinkTest, MaskPreservesUnmaskedArea)
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

  Array out = gpu::voronoi_shrink(input,
                                  glm::vec2(6.f, 6.f),
                                  0.5f,
                                  -5.f,
                                  1,
                                  {0.5f, 0.5f},
                                  0.f,
                                  &mask);

  // Bottom half (mask = 0) should remain unchanged
  for (int j = 0; j < ny / 2; ++j)
    for (int i = 0; i < nx; ++i)
      EXPECT_NEAR(out(i, j), input(i, j), 1e-5f);
}

TEST_F(VoronoiShrinkTest, NoisePerturbationModifiesOutput)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  Array      noise_x = noise_fbm(NoiseType::PERLIN, shape, {4.f, 4.f}, 1);
  Array      noise_y = noise_fbm(NoiseType::PERLIN, shape, {4.f, 4.f}, 2);

  Array out_unperturbed = gpu::voronoi_shrink(input,
                                              glm::vec2(4.f, 4.f),
                                              0.7f,
                                              0.f,
                                              1,
                                              {0.5f, 0.5f});
  Array out_perturbed = gpu::voronoi_shrink(input,
                                            glm::vec2(4.f, 4.f),
                                            0.7f,
                                            0.f,
                                            1,
                                            {0.5f, 0.5f},
                                            0.f,
                                            nullptr,
                                            &noise_x,
                                            &noise_y);

  EXPECT_FALSE(assert_almost_equal(out_unperturbed, out_perturbed, 1e-3f));
}

TEST_F(VoronoiShrinkTest, EmptyArrayReturnsEmpty)
{
  Array input;
  Array out = gpu::voronoi_shrink(input, 4.f);
  EXPECT_TRUE(out.vector.empty());
}

// --- Point list / Cloud tests

TEST_F(VoronoiShrinkTest, PointsShrinkFactorOnePreservesInput)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);

  std::vector<glm::vec2> pts = {{0.25f, 0.25f},
                                {0.75f, 0.25f},
                                {0.25f, 0.75f},
                                {0.75f, 0.75f}};

  Array out = gpu::voronoi_shrink(input, pts, 1.0f, -999.f);

  EXPECT_EQ(out.shape, shape);
  EXPECT_TRUE(assert_almost_equal(out, input, 1e-4f));
}

TEST_F(VoronoiShrinkTest, PointsShrinkingCreatesRiftsWithFillValue)
{
  glm::ivec2 shape = {128, 128};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  float      fill_val = -10.f;

  std::vector<glm::vec2> pts = {{0.2f, 0.2f},
                                {0.8f, 0.2f},
                                {0.5f, 0.5f},
                                {0.2f, 0.8f},
                                {0.8f, 0.8f}};

  Array out = gpu::voronoi_shrink(input, pts, 0.5f, fill_val);

  int rifts_count = 0;
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      if (std::abs(out(i, j) - fill_val) < 1e-5f)
      {
        rifts_count++;
      }
    }

  EXPECT_GT(rifts_count, 0);
  EXPECT_LT(rifts_count, shape.x * shape.y);
}

TEST_F(VoronoiShrinkTest, PointsOverloadEquivalence)
{
  glm::ivec2 shape = {64, 64};
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);

  std::vector<glm::vec2> pts_vec = {{0.3f, 0.3f}, {0.7f, 0.7f}, {0.3f, 0.7f}};
  std::vector<Point>     pts_obj = {Point(pts_vec[0]),
                                    Point(pts_vec[1]),
                                    Point(pts_vec[2])};
  Cloud                  cloud(pts_obj);

  Array out_vec = gpu::voronoi_shrink(input, pts_vec, 0.6f, -1.f);
  Array out_obj = gpu::voronoi_shrink(input, pts_obj, 0.6f, -1.f);
  Array out_cld = gpu::voronoi_shrink(input, cloud, 0.6f, -1.f);

  EXPECT_TRUE(assert_almost_equal(out_vec, out_obj, 1e-5f));
  EXPECT_TRUE(assert_almost_equal(out_vec, out_cld, 1e-5f));
}

TEST_F(VoronoiShrinkTest, PointsMaskPreservesUnmaskedArea)
{
  glm::ivec2 shape = {64, 64};
  int        nx = shape.x;
  int        ny = shape.y;
  Array      input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);

  // Mask: 0 on bottom half, 1 on top half
  Array mask(shape);
  for (int j = 0; j < ny; ++j)
    for (int i = 0; i < nx; ++i)
      mask(i, j) = (j >= ny / 2) ? 1.f : 0.f;

  std::vector<glm::vec2> pts = {{0.3f, 0.3f}, {0.7f, 0.7f}};
  Array out = gpu::voronoi_shrink(input, pts, 0.5f, -5.f, &mask);

  // Bottom half (mask = 0) should remain unchanged
  for (int j = 0; j < ny / 2; ++j)
    for (int i = 0; i < nx; ++i)
      EXPECT_NEAR(out(i, j), input(i, j), 1e-5f);
}

TEST_F(VoronoiShrinkTest, PointsEmptyCloudReturnsUnchangedShape)
{
  Array input(glm::ivec2(32, 32));
  Cloud empty_cloud;
  Array out = gpu::voronoi_shrink(input, empty_cloud, 0.5f);
  EXPECT_EQ(out.shape, input.shape);
}
