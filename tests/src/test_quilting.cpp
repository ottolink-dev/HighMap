#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "highmap/operator.hpp"
#include "highmap/primitives.hpp"
#include "highmap/range.hpp"
#include "highmap/synthesis.hpp"

#include <gtest/gtest.h>

using namespace hmap;

namespace
{

Array make_sample_array(glm::ivec2 shape, std::uint32_t seed)
{
  Array z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(4.f, 4.f), seed);
  remap(z);
  return z;
}

} // namespace

TEST(Quilting, InvalidParametersReturnEmptyArray)
{
  Array input = make_sample_array({64, 64}, 1);

  EXPECT_TRUE(quilting({}, {16, 16}, {2, 2}, 0.5f, 1).vector.empty())
      << "empty input array list";
  EXPECT_TRUE(quilting({&input}, {0, 16}, {2, 2}, 0.5f, 1).vector.empty())
      << "invalid patch_base_shape";
  EXPECT_TRUE(quilting({&input}, {16, 16}, {0, 2}, 0.5f, 1).vector.empty())
      << "invalid tiling";

  Array empty_arr;
  EXPECT_TRUE(quilting({&empty_arr}, {16, 16}, {2, 2}, 0.5f, 1).vector.empty())
      << "empty array in input list";

  EXPECT_TRUE(quilting_blend({}, {16, 16}, 0.5f, 1).vector.empty());
  EXPECT_TRUE(
      quilting_expand(empty_arr, 2.f, {16, 16}, 0.5f, 1).vector.empty());
  EXPECT_TRUE(quilting_shuffle(empty_arr, {16, 16}, 0.5f, 1).vector.empty());
}

TEST(Quilting, BasicQuiltingAndShape)
{
  Array input = make_sample_array({64, 64}, 1);

  glm::ivec2 patch_base_shape = {16, 16};
  glm::ivec2 tiling = {3, 3};
  float      overlap = 0.5f;

  Array out = quilting({&input}, patch_base_shape, tiling, overlap, 42);

  glm::ivec2 patch_shape = {(int)(patch_base_shape.x * (1.f + overlap)),
                            (int)(patch_base_shape.y * (1.f + overlap))};
  glm::ivec2 expected_shape = {
      tiling.x * patch_base_shape.x + (patch_shape.x - patch_base_shape.x),
      tiling.y * patch_base_shape.y + (patch_shape.y - patch_base_shape.y)};

  EXPECT_EQ(out.shape.x, expected_shape.x);
  EXPECT_EQ(out.shape.y, expected_shape.y);

  for (float v : out.vector)
  {
    EXPECT_TRUE(std::isfinite(v));
  }
}

TEST(Quilting, SecondaryArrays)
{
  Array input = make_sample_array({64, 64}, 1);
  Array sec1 = make_sample_array({64, 64}, 2);
  Array sec2 = make_sample_array({64, 64}, 3);

  std::vector<Array *> sec_ptrs = {&sec1, &sec2};

  Array out = quilting({&input}, {16, 16}, {2, 2}, 0.5f, 10, sec_ptrs);

  EXPECT_FALSE(out.vector.empty());
  EXPECT_EQ(sec1.shape.x, out.shape.x);
  EXPECT_EQ(sec1.shape.y, out.shape.y);
  EXPECT_EQ(sec2.shape.x, out.shape.x);
  EXPECT_EQ(sec2.shape.y, out.shape.y);

  for (float v : sec1.vector)
    EXPECT_TRUE(std::isfinite(v));
  for (float v : sec2.vector)
    EXPECT_TRUE(std::isfinite(v));
}

TEST(Quilting, QuiltingBlend)
{
  Array a1 = make_sample_array({64, 64}, 1);
  Array a2 = make_sample_array({64, 64}, 2);

  Array out = quilting_blend({&a1, &a2}, {16, 16}, 0.5f, 100);

  EXPECT_EQ(out.shape.x, a2.shape.x);
  EXPECT_EQ(out.shape.y, a2.shape.y);

  for (float v : out.vector)
    EXPECT_TRUE(std::isfinite(v));
}

TEST(Quilting, QuiltingExpand)
{
  Array input = make_sample_array({64, 64}, 1);

  Array out_exp = quilting_expand(input, 2.f, {16, 16}, 0.5f, 7, {}, false);
  EXPECT_EQ(out_exp.shape.x, 128);
  EXPECT_EQ(out_exp.shape.y, 128);

  Array out_keep = quilting_expand(input, 2.f, {16, 16}, 0.5f, 7, {}, true);
  EXPECT_EQ(out_keep.shape.x, 64);
  EXPECT_EQ(out_keep.shape.y, 64);
}

TEST(Quilting, QuiltingShuffle)
{
  Array input = make_sample_array({64, 64}, 1);

  Array out = quilting_shuffle(input, {16, 16}, 0.5f, 50);

  EXPECT_EQ(out.shape.x, input.shape.x);
  EXPECT_EQ(out.shape.y, input.shape.y);

  for (float v : out.vector)
    EXPECT_TRUE(std::isfinite(v));
}

TEST(Quilting, FindVerticalCutPathOptimality)
{
  // Construct a known error matrix where the optimal path is unique
  // Cost matrix (5x4):
  // 10 10  0 10 10
  // 10  0 10 10 10
  //  0 10 10 10 10
  // 10  0 10 10 10
  Array error(glm::ivec2(5, 4));
  std::fill(error.vector.begin(), error.vector.end(), 10.f);
  error(2, 0) = 0.f;
  error(1, 1) = 0.f;
  error(0, 2) = 0.f;
  error(1, 3) = 0.f;

  std::vector<int> path;
  find_vertical_cut_path(error, path);

  ASSERT_EQ(path.size(), 4u);
  EXPECT_EQ(path[0], 2);
  EXPECT_EQ(path[1], 1);
  EXPECT_EQ(path[2], 0);
  EXPECT_EQ(path[3], 1);
}

TEST(Quilting, ZeroOverlapQuilting)
{
  Array input = make_sample_array({64, 64}, 1);

  Array out = quilting({&input}, {16, 16}, {2, 2}, 0.0f, 42);
  EXPECT_EQ(out.shape.x, 32);
  EXPECT_EQ(out.shape.y, 32);

  for (float v : out.vector)
    EXPECT_TRUE(std::isfinite(v));
}
