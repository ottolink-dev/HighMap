#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "highmap/blending.hpp"
#include "highmap/coord_frame.hpp"
#include "highmap/filters.hpp"
#include "highmap/internal/assert.hpp"
#include "highmap/primitives.hpp"
#include "highmap/range.hpp"
#include "highmap/synthesis.hpp"
#include "highmap/virtual_array.hpp"

#include <gtest/gtest.h>

using namespace hmap;

namespace
{

Array make_exemplar(glm::ivec2 shape, std::uint32_t seed)
{
  Array z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(8.f, 8.f), seed);
  remap(z);
  return z;
}

// Same block-mean pooling the implementation uses to derive the
// low-resolution exemplar, so feeding the pooled exemplar back as the input
// must reproduce the exemplar.
Array block_mean_downsample(const Array &array, int factor)
{
  Array out(glm::ivec2(array.shape.x / factor, array.shape.y / factor));

  for (int j = 0; j < out.shape.y; j++)
    for (int i = 0; i < out.shape.x; i++)
    {
      float sum = 0.f;
      for (int q = 0; q < factor; q++)
        for (int p = 0; p < factor; p++)
          sum += array(i * factor + p, j * factor + q);
      out(i, j) = sum / (float)(factor * factor);
    }

  return out;
}

float max_abs_diff_interior(const Array &a, const Array &b, int margin)
{
  float err = 0.f;
  for (int j = margin; j < a.shape.y - margin; j++)
    for (int i = margin; i < a.shape.x - margin; i++)
      err = std::max(err, std::abs(a(i, j) - b(i, j)));
  return err;
}

} // namespace

TEST(TerrainSuperResolution, OutputShapeIsInputShapeTimesFactor)
{
  Array exemplar = make_exemplar({128, 128}, 1);
  Array input = make_exemplar({40, 24}, 2);

  Array out = terrain_super_resolution(input, exemplar, 4, 8, 4, 4);

  EXPECT_EQ(out.shape.x, 160);
  EXPECT_EQ(out.shape.y, 96);
}

TEST(TerrainSuperResolution, InvalidParametersReturnEmptyArray)
{
  Array exemplar = make_exemplar({64, 64}, 1);
  Array input = make_exemplar({16, 16}, 2);

  // valid baseline, so that the checks below fail for the intended reason
  ASSERT_FALSE(
      terrain_super_resolution(input, exemplar, 2, 8, 4, 4).vector.empty());

  EXPECT_TRUE(
      terrain_super_resolution(input, exemplar, 0, 8, 4, 4).vector.empty())
      << "factor < 1";
  EXPECT_TRUE(
      terrain_super_resolution(input, exemplar, 2, 1, 1, 1).vector.empty())
      << "patch_size < 2";
  EXPECT_TRUE(
      terrain_super_resolution(input, exemplar, 2, 8, 0, 4).vector.empty())
      << "analysis_stride < 1";
  EXPECT_TRUE(
      terrain_super_resolution(input, exemplar, 2, 8, 4, 9).vector.empty())
      << "synthesis_stride > patch_size";
  EXPECT_TRUE(
      terrain_super_resolution(input, exemplar, 2, 8, 4, 4, 0).vector.empty())
      << "sparsity < 1";
  // a 64x64 exemplar downsampled by 8 is 8x8, smaller than a 16-pixel patch
  EXPECT_TRUE(
      terrain_super_resolution(input, exemplar, 8, 16, 8, 8).vector.empty())
      << "exemplar smaller than patch";
  EXPECT_TRUE(
      terrain_super_resolution(Array(), exemplar, 2, 8, 4, 4).vector.empty())
      << "empty input";
  EXPECT_TRUE(
      terrain_super_resolution(input, Array(), 2, 8, 4, 4).vector.empty())
      << "empty exemplar";
}

TEST(TerrainSuperResolution, FlatExemplarReturnsEmptyArray)
{
  Array exemplar = Array(glm::ivec2(64, 64), 0.5f);
  Array input = make_exemplar({16, 16}, 2);

  EXPECT_TRUE(
      terrain_super_resolution(input, exemplar, 2, 8, 4, 4).vector.empty());
}

TEST(TerrainSuperResolution, ReconstructsExemplarFromItsOwnDownsampling)
{
  const int factor = 4;
  const int patch_size = 8;
  const int stride = 4;

  Array exemplar = make_exemplar({128, 128}, 3);
  Array input = block_mean_downsample(exemplar, factor);

  Array out = terrain_super_resolution(input,
                                       exemplar,
                                       factor,
                                       patch_size,
                                       stride,
                                       stride);

  ASSERT_EQ(out.shape.x, exemplar.shape.x);
  ASSERT_EQ(out.shape.y, exemplar.shape.y);

  // patches overlapping the padded border are not exemplar patches, so only
  // the interior is expected to match
  const int margin = patch_size * factor;
  EXPECT_LT(max_abs_diff_interior(out, exemplar, margin), 1e-3f);
}

TEST(TerrainSuperResolution, SparsityAboveOneStillReconstructsExactPatches)
{
  const int factor = 4;
  const int patch_size = 8;
  const int stride = 4;
  const int sparsity = 3;

  Array exemplar = make_exemplar({128, 128}, 3);
  Array input = block_mean_downsample(exemplar, factor);

  Array out = terrain_super_resolution(input,
                                       exemplar,
                                       factor,
                                       patch_size,
                                       stride,
                                       stride,
                                       sparsity);

  ASSERT_EQ(out.shape.x, exemplar.shape.x);
  ASSERT_EQ(out.shape.y, exemplar.shape.y);

  const int margin = patch_size * factor;
  EXPECT_LT(max_abs_diff_interior(out, exemplar, margin), 1e-3f);
}

TEST(TerrainSuperResolution, GenericRunIsFiniteAndDeterministic)
{
  Array exemplar = make_exemplar({256, 256}, 5);
  Array input = noise_fbm(NoiseType::PERLIN, {48, 32}, {2.f, 2.f}, 6, 3);
  remap(input);

  Array out1 = terrain_super_resolution(input, exemplar, 4, 16, 8, 8, 2);
  Array out2 = terrain_super_resolution(input, exemplar, 4, 16, 8, 8, 2);

  ASSERT_EQ(out1.vector.size(), (size_t)(192 * 128));
  ASSERT_EQ(out2.vector.size(), out1.vector.size());

  size_t n_not_finite = 0;
  size_t n_different = 0;
  for (size_t k = 0; k < out1.vector.size(); k++)
  {
    if (!std::isfinite(out1.vector[k])) n_not_finite++;
    if (out1.vector[k] != out2.vector[k]) n_different++;
  }

  EXPECT_EQ(n_not_finite, 0u);
  EXPECT_EQ(n_different, 0u);
}

TEST(TransferTerrainFeatures, PreservesMacroShapeAndAddsDetail)
{
  const glm::ivec2 shape(64, 64);
  CoordFrame       frame_target({0.f, 0.f}, {1.f, 1.f}, 0.f);
  CoordFrame       frame_source({5.f, 5.f}, {1.f, 1.f}, 0.f);

  Array target_base = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 10, 2);
  remap(target_base, 0.f, 1.f);

  Array source = noise_fbm(NoiseType::PERLIN, shape, {8.f, 8.f}, 20, 6);
  remap(source, 0.f, 1.f);

  Array target = target_base;
  transfer_terrain_features(target,
                            frame_target,
                            source,
                            frame_source,
                            8,
                            0.5f,
                            false,
                            2,
                            8,
                            4,
                            4,
                            2,
                            256,
                            1);

  // target must have been modified with high-frequency details
  EXPECT_FALSE(assert_almost_equal(target, target_base, 1e-4f));

  // large-scale low-pass structure must be preserved
  Array low_before = target_base;
  Array low_after = target;
  smooth_cpulse(low_before, 8);
  smooth_cpulse(low_after, 8);

  EXPECT_LT(max_abs_diff_interior(low_before, low_after, 8), 0.05f);
}

TEST(TransferTerrainFeatures, IgnoresOriginAndAccountsForCoordFrameScale)
{
  const glm::ivec2 shape(64, 64);
  CoordFrame       frame_target({0.f, 0.f}, {1.f, 1.f}, 0.f);
  CoordFrame       frame_src_a({0.f, 0.f}, {1.f, 1.f}, 0.f);
  CoordFrame       frame_src_b({10.f, -7.f}, {1.f, 1.f}, 45.f);
  CoordFrame       frame_src_scaled({0.f, 0.f}, {0.5f, 0.5f}, 0.f);

  Array base = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 1, 2);
  Array source = noise_fbm(NoiseType::PERLIN, shape, {8.f, 8.f}, 2, 5);

  Array out_a = base;
  transfer_terrain_features(out_a,
                            frame_target,
                            source,
                            frame_src_a,
                            8,
                            1.f,
                            false,
                            2,
                            8,
                            4,
                            4,
                            2);

  Array out_b = base;
  transfer_terrain_features(out_b,
                            frame_target,
                            source,
                            frame_src_b,
                            8,
                            1.f,
                            false,
                            2,
                            8,
                            4,
                            4,
                            2);

  // spatial position / origin does not affect exemplar transfer
  EXPECT_TRUE(assert_almost_equal(out_a, out_b, 1e-6f));

  Array out_scaled = base;
  transfer_terrain_features(out_scaled,
                            frame_target,
                            source,
                            frame_src_scaled,
                            8,
                            1.f,
                            false,
                            2,
                            8,
                            4,
                            4,
                            2);

  // different CoordFrame physical size changes the transferred feature scale
  EXPECT_FALSE(assert_almost_equal(out_a, out_scaled, 1e-4f));
}

TEST(TransferTerrainFeatures, MultiSourceWeightsAndBoundedDictionary)
{
  const glm::ivec2 shape(64, 64);
  CoordFrame       frame_target({0.f, 0.f}, {1.f, 1.f}, 0.f);
  CoordFrame       frame_s1({0.f, 0.f}, {1.f, 1.f}, 0.f);
  CoordFrame       frame_s2({0.f, 0.f}, {1.f, 1.f}, 0.f);

  Array base = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 3, 2);
  Array src1 = noise_fbm(NoiseType::PERLIN, shape, {6.f, 6.f}, 4, 5);
  Array src2 = noise_fbm(NoiseType::PERLIN, shape, {10.f, 10.f}, 5, 5);

  Array out_single = base;
  transfer_terrain_features(out_single,
                            frame_target,
                            {&src1},
                            {&frame_s1},
                            8,
                            1.f,
                            false,
                            {1.f},
                            2,
                            8,
                            4,
                            4,
                            2,
                            512,
                            42);

  Array out_zero_second = base;
  transfer_terrain_features(out_zero_second,
                            frame_target,
                            {&src1, &src2},
                            {&frame_s1, &frame_s2},
                            8,
                            1.f,
                            false,
                            {1.f, 0.f},
                            2,
                            8,
                            4,
                            4,
                            2,
                            512,
                            42);

  // source with weight 0 is skipped
  EXPECT_TRUE(assert_almost_equal(out_single, out_zero_second, 1e-6f));

  // small max_atoms reservoir sampling is deterministic for the same seed
  Array out_cap1 = base;
  Array out_cap2 = base;

  transfer_terrain_features(out_cap1,
                            frame_target,
                            {&src1, &src2},
                            {&frame_s1, &frame_s2},
                            8,
                            1.f,
                            false,
                            {1.f, 0.5f},
                            2,
                            8,
                            2,
                            4,
                            2,
                            10,
                            99);

  transfer_terrain_features(out_cap2,
                            frame_target,
                            {&src1, &src2},
                            {&frame_s1, &frame_s2},
                            8,
                            1.f,
                            false,
                            {1.f, 0.5f},
                            2,
                            8,
                            2,
                            4,
                            2,
                            10,
                            99);

  EXPECT_TRUE(assert_almost_equal(out_cap1, out_cap2, 1e-6f));
}

TEST(TransferTerrainFeatures, VirtualArrayConsistencyAndTiling)
{
  const glm::ivec2 shape(64, 64);
  const glm::ivec2 tile_shape(32, 32);
  CoordFrame       frame_target({0.f, 0.f}, {1.f, 1.f}, 0.f);
  CoordFrame       frame_source({0.f, 0.f}, {1.f, 1.f}, 0.f);

  Array base = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 7, 2);
  Array source = noise_fbm(NoiseType::PERLIN, shape, {8.f, 8.f}, 8, 5);

  Array out_arr = base;
  transfer_terrain_features(out_arr,
                            frame_target,
                            source,
                            frame_source,
                            4,
                            0.5f,
                            false,
                            2,
                            8,
                            4,
                            4,
                            2);

  // single-tile / VA_SINGLE_ARRAY matches Array exactly
  VirtualArray target_va_single(shape, shape, 0, StorageMode::VA_RAM);
  target_va_single.from_array(base, ComputeMode{});
  VirtualArray source_va(shape, tile_shape, 0, StorageMode::VA_RAM);
  source_va.from_array(source, ComputeMode{});

  transfer_terrain_features(target_va_single,
                            frame_target,
                            source_va,
                            frame_source,
                            4,
                            0.5f,
                            false,
                            2,
                            8,
                            4,
                            4,
                            2,
                            1024,
                            0,
                            ComputeMode{.mode = ForEachMode::VA_SINGLE_ARRAY});

  Array res_single = target_va_single.to_array(ComputeMode{});
  EXPECT_TRUE(assert_almost_equal(out_arr, res_single, 1e-5f));

  // multi-tile VA with halo synchronization in distributed mode
  VirtualArray target_va_tiled(shape, tile_shape, 4, StorageMode::VA_RAM);
  target_va_tiled.from_array(base, ComputeMode{});

  transfer_terrain_features(target_va_tiled,
                            frame_target,
                            source_va,
                            frame_source,
                            4,
                            0.5f,
                            false,
                            2,
                            8,
                            4,
                            4,
                            2,
                            1024,
                            0,
                            ComputeMode{.mode = ForEachMode::VA_DISTRIBUTED});

  Array res_tiled = target_va_tiled.to_array(ComputeMode{});
  EXPECT_EQ(res_tiled.shape, shape);
  EXPECT_FALSE(assert_almost_equal(res_tiled, base, 1e-4f));
}
