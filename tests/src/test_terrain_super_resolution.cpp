#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "highmap/primitives.hpp"
#include "highmap/range.hpp"
#include "highmap/synthesis.hpp"

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
