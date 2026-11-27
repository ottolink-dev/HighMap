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

TEST(Quilting, SecondaryChannelsShareGuideSeamsAndTransforms)
{
  for (glm::ivec2 shape : {glm::ivec2(48, 48), glm::ivec2(48, 32)})
    for (int flags = 0; flags < 8; ++flags)
      for (float overlap : {0.f, 0.5f, 0.9f})
      {
        SCOPED_TRACE(::testing::Message() << shape.x << "x" << shape.y << " flags="
                                          << flags << " overlap=" << overlap);
        Array      guide = make_sample_array(shape, 91);
        Array      source = make_sample_array(shape, 37);
        glm::ivec2 patch = shape.x == shape.y ? glm::ivec2(12, 12) : glm::ivec2(12, 8);
        // A separate one-channel call checks every secondary cell, including
        // horizontal and vertical seams, against a channel-count invariant.
        Array              expected = source;
        Array              expected_guide = quilting({&guide},
                                        patch,
                                                     {3, 4},
                                        overlap,
                                        12,
                                                     {&expected},
                                        flags & 1,
                                        flags & 2,
                                        flags & 4,
                                        0.5f);
        std::vector<Array> channels(20, source);
        channels[0] = guide;
        std::vector<Array *> pointers;
        for (auto &channel : channels)
          pointers.push_back(&channel);
        Array actual = quilting({&guide},
                                patch,
                                {3, 4},
                                overlap,
                                12,
                                pointers,
                                flags & 1,
                                flags & 2,
                                flags & 4,
                                0.5f);
        EXPECT_EQ(actual.vector, expected_guide.vector);
        EXPECT_EQ(channels[0].vector, actual.vector);
        for (size_t k = 1; k < channels.size(); ++k)
          EXPECT_EQ(channels[k].vector, expected.vector) << "channel " << k;
      }
}

TEST(Quilting, SecondaryAliasesAreNotOverwrittenUntilAllSourcesAreConsumed)
{
  for (int flags = 0; flags < 8; ++flags)
  {
    SCOPED_TRACE(flags);
    Array guide = make_sample_array({48, 32}, 91);
    Array secondary = make_sample_array({48, 32}, 37);
    Array expected_guide_channel = guide;
    Array expected_secondary = secondary;
    Array expected = quilting({&guide},
                              {12, 8},
                              {3, 4},
                              0.9f,
                              12,
                              {&expected_guide_channel, &expected_secondary},
                              flags & 1,
                              flags & 2,
                              flags & 4,
                              0.5f);
    // The guide itself is also a destination. Repeated secondary destinations
    // must still read their original input, even after earlier rows finish.
    Array actual = quilting({&guide},
                            {12, 8},
                            {3, 4},
                            0.9f,
                            12,
                            {&secondary, &guide, &secondary, &guide},
                            flags & 1,
                            flags & 2,
                            flags & 4,
                            0.5f);
    EXPECT_EQ(actual.vector, expected.vector);
    EXPECT_EQ(guide.vector, expected_guide_channel.vector);
    EXPECT_EQ(secondary.vector, expected_secondary.vector);
  }
}

TEST(Quilting, ExpandAndShufflePreserveSecondaryChannelAlignment)
{
  for (bool expand : {false, true})
    for (bool keep_shape : {false, true})
      for (float ratio : {1.5f, 2.f})
      {
        SCOPED_TRACE(::testing::Message() << "expand=" << expand << " keep=" << keep_shape
                                          << " ratio=" << ratio);
        Array                guide = make_sample_array({48, 32}, 91);
        std::vector<Array>   channels(4, guide);
        std::vector<Array *> pointers;
        for (auto &channel : channels)
          pointers.push_back(&channel);
        Array actual = expand ? quilting_expand(guide,
                                                ratio,
                                                {12, 8},
                                                0.9f,
                                                12,
                                                pointers,
                                                keep_shape,
                                                true,
                                                true,
                                                true,
                                                0.5f)
                              : quilting_shuffle(guide,
                                                 {12, 8},
                                                 0.9f,
                                                 12,
                                                 pointers,
                                                 true,
                                                 true,
                                                 true,
                                                 0.5f);
        ASSERT_FALSE(actual.vector.empty());
        for (const auto &channel : channels)
        {
          EXPECT_EQ(channel.shape, actual.shape);
          EXPECT_EQ(channel.vector, actual.vector);
          for (float value : channel.vector)
            EXPECT_TRUE(std::isfinite(value));
        }
      }
}

TEST(Quilting, MismatchedSecondaryShapesReturnEmptyWithoutMutation)
{
  Array guide = make_sample_array({48, 32}, 91);
  Array other_guide = make_sample_array({32, 32}, 92);
  Array secondary = make_sample_array({48, 32}, 37);
  Array small = make_sample_array({32, 32}, 38);
  const Array original_guide = guide;
  const Array original_secondary = secondary;
  const Array original_small = small;

  EXPECT_TRUE(quilting({&guide}, {12, 8}, {3, 4}, 0.9f, 12,
                       {&guide, &secondary, &small}).vector.empty());
  EXPECT_EQ(guide.vector, original_guide.vector);
  EXPECT_EQ(secondary.vector, original_secondary.vector);
  EXPECT_EQ(small.vector, original_small.vector);

  // Reject a mismatch with any possible guide, rather than depending on which
  // guide a particular random seed selects first.
  EXPECT_TRUE(quilting({&guide, &other_guide}, {12, 8}, {3, 4}, 0.9f, 12,
                       {&secondary}).vector.empty());
  EXPECT_EQ(secondary.vector, original_secondary.vector);
}

TEST(Quilting, RandomPatchSecondaryOwnershipPreservesSources)
{
  for (glm::ivec2 patch_shape : {glm::ivec2(12, 12), glm::ivec2(12, 8)})
    for (int flags = 0; flags < 8; ++flags)
    {
      Array guide = make_sample_array({48, 32}, 91);
      Array source = make_sample_array({48, 32}, 37);
      const Array original_source = source;
      std::mt19937 gen(12);
      std::vector<Array *> pointers(20, &source);
      std::vector<Array> patches(1, guide);
      get_random_patch(guide, patch_shape, gen, flags & 1, flags & 2,
                       flags & 4, &pointers, &patches);
      std::mt19937 reference_gen(12);
      Array expected = get_random_patch(source, patch_shape, reference_gen,
                                        flags & 1, flags & 2, flags & 4);
      ASSERT_EQ(patches.size(), pointers.size());
      for (const auto &patch : patches)
      {
        EXPECT_EQ(patch.shape, expected.shape);
        EXPECT_EQ(patch.vector, expected.vector);
      }
      EXPECT_EQ(source.vector, original_source.vector);
    }
}

TEST(Quilting, SubpixelSmoothingRadiusProducesFiniteSeams)
{
  for (float width : {0.f, 0.2f})
  {
    Array guide = make_sample_array({48, 48}, 91);
    std::vector<Array> channels(4, guide);
    std::vector<Array *> pointers;
    for (auto &channel : channels) pointers.push_back(&channel);
    // Expansion gives a six-cell patch with three overlap cells. Both widths
    // round the smoothing radius to zero, so the cut mask stays unsmoothed.
    Array output = quilting_expand(guide, 2.f, {12, 12}, 0.5f, 123,
                                    pointers, true, false, false, false, width);
    ASSERT_EQ(output.shape, guide.shape);
    for (float value : output.vector) EXPECT_TRUE(std::isfinite(value));
    for (const auto &channel : channels)
      EXPECT_EQ(channel.vector, output.vector);
  }
}
