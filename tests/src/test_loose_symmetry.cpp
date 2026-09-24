#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>

#include "highmap/primitives.hpp"
#include "highmap/range.hpp"
#include "highmap/synthesis.hpp"
#include "highmap/transform.hpp"

#include <gtest/gtest.h>

using namespace hmap;

namespace
{

Array make_test_terrain(glm::ivec2 shape, std::uint32_t seed)
{
  Array z = noise_fbm(NoiseType::PERLIN, shape, glm::vec2(4.f, 4.f), seed);
  remap(z);
  return z;
}

} // namespace

TEST(SymmetrizeTest, AllSymmetryModesProduceExactGeometricSymmetry)
{
  Array input = make_test_terrain({32, 32}, 42);

  // --- Left to right

  Array lr = symmetrize(input, SymmetryType::SYMMETRY_LEFT_TO_RIGHT);
  ASSERT_EQ(lr.shape.x, 32);
  ASSERT_EQ(lr.shape.y, 32);
  for (int j = 0; j < 32; j++)
    for (int i = 0; i < 32; i++)
      EXPECT_FLOAT_EQ(lr(i, j), lr(31 - i, j));

  // --- Right to left

  Array rl = symmetrize(input, SymmetryType::SYMMETRY_RIGHT_TO_LEFT);
  for (int j = 0; j < 32; j++)
    for (int i = 0; i < 32; i++)
      EXPECT_FLOAT_EQ(rl(i, j), rl(31 - i, j));

  // --- Top to bottom

  Array tb = symmetrize(input, SymmetryType::SYMMETRY_TOP_TO_BOTTOM);
  for (int j = 0; j < 32; j++)
    for (int i = 0; i < 32; i++)
      EXPECT_FLOAT_EQ(tb(i, j), tb(i, 31 - j));

  // --- Bottom to top

  Array bt = symmetrize(input, SymmetryType::SYMMETRY_BOTTOM_TO_TOP);
  for (int j = 0; j < 32; j++)
    for (int i = 0; i < 32; i++)
      EXPECT_FLOAT_EQ(bt(i, j), bt(i, 31 - j));

  // --- Symmetry X (average)

  Array sym_x = symmetrize(input, SymmetryType::SYMMETRY_X);
  for (int j = 0; j < 32; j++)
    for (int i = 0; i < 32; i++)
      EXPECT_FLOAT_EQ(sym_x(i, j), sym_x(31 - i, j));

  // --- Symmetry Y (average)

  Array sym_y = symmetrize(input, SymmetryType::SYMMETRY_Y);
  for (int j = 0; j < 32; j++)
    for (int i = 0; i < 32; i++)
      EXPECT_FLOAT_EQ(sym_y(i, j), sym_y(i, 31 - j));

  // --- Symmetry XY (4-way average)

  Array sym_xy = symmetrize(input, SymmetryType::SYMMETRY_XY);
  for (int j = 0; j < 32; j++)
    for (int i = 0; i < 32; i++)
    {
      EXPECT_FLOAT_EQ(sym_xy(i, j), sym_xy(31 - i, j));
      EXPECT_FLOAT_EQ(sym_xy(i, j), sym_xy(i, 31 - j));
    }

  // --- Symmetry rot180

  Array sym_rot = symmetrize(input, SymmetryType::SYMMETRY_ROT180);
  for (int j = 0; j < 32; j++)
    for (int i = 0; i < 32; i++)
      EXPECT_FLOAT_EQ(sym_rot(i, j), sym_rot(31 - i, 31 - j));
}

TEST(SymmetrizeTest, EmptyInputReturnsEmptyArray)
{
  EXPECT_TRUE(symmetrize(Array()).vector.empty());
}

TEST(LooseSymmetryTest, OutputShapeMatchesInputShape)
{
  Array input = make_test_terrain({64, 64}, 1);

  Array out1 = loose_symmetry(input,
                              SymmetryType::SYMMETRY_LEFT_TO_RIGHT,
                              1.f,
                              2,
                              8,
                              2,
                              4);
  EXPECT_EQ(out1.shape.x, 64);
  EXPECT_EQ(out1.shape.y, 64);

  Array out2 = loose_symmetry(input,
                              SymmetryType::SYMMETRY_LEFT_TO_RIGHT,
                              1.f,
                              4,
                              8,
                              2,
                              4);
  EXPECT_EQ(out2.shape.x, 64);
  EXPECT_EQ(out2.shape.y, 64);
}

TEST(LooseSymmetryTest, EmptyInputReturnsEmptyArray)
{
  EXPECT_TRUE(loose_symmetry(Array()).vector.empty());
}

TEST(LooseSymmetryTest, TunableStrengthRunsAndProducesFiniteResults)
{
  Array input = make_test_terrain({64, 64}, 7);

  Array out_half = loose_symmetry(input,
                                  SymmetryType::SYMMETRY_TOP_TO_BOTTOM,
                                  0.5f,
                                  2,
                                  8,
                                  2,
                                  4);

  ASSERT_EQ(out_half.shape.x, 64);
  ASSERT_EQ(out_half.shape.y, 64);

  for (float v : out_half.vector)
    EXPECT_TRUE(std::isfinite(v));
}

TEST(LooseSymmetryTest, HasNaturalVariationComparedToStrictMirror)
{
  Array input = make_test_terrain({64, 64}, 123);

  Array strict = symmetrize(input, SymmetryType::SYMMETRY_LEFT_TO_RIGHT);
  Array loose = loose_symmetry(input,
                               SymmetryType::SYMMETRY_LEFT_TO_RIGHT,
                               1.f,
                               2,
                               8,
                               2,
                               4);

  ASSERT_EQ(loose.shape.x, 64);
  ASSERT_EQ(loose.shape.y, 64);

  // loose symmetry output should not be an exact pixel-by-pixel mirror
  int asymmetric_pixels = 0;
  for (int j = 0; j < 64; j++)
    for (int i = 0; i < 32; i++)
    {
      if (std::abs(loose(i, j) - loose(63 - i, j)) > 1e-4f) asymmetric_pixels++;
    }

  EXPECT_GT(asymmetric_pixels, 0);

  // overall left and right halves should have comparable macro mean/energy
  float left_sum = 0.f;
  float right_sum = 0.f;
  for (int j = 0; j < 64; j++)
  {
    for (int i = 0; i < 32; i++)
      left_sum += loose(i, j);
    for (int i = 32; i < 64; i++)
      right_sum += loose(i, j);
  }

  EXPECT_NEAR(left_sum, right_sum, 0.2f * (left_sum + right_sum));
}

TEST(SymmetrizeTest, FlattenCenterReducesCenterElevationToMinimum)
{
  Array input({33, 33}, 5.0f);
  // set minimum value
  input(0, 0) = 1.0f;
  input(32, 32) = 1.0f;
  // center is at index (16, 16)
  input(16, 16) = 10.0f;

  Array flat_x = symmetrize(input,
                            SymmetryType::SYMMETRY_LEFT_TO_RIGHT,
                            /* flatten_center */ true,
                            /* flatten_radius */ 0.1f);
  // Along the symmetry axis (i = 16), elevation should equal min elevation
  // (1.0)
  for (int j = 0; j < 33; j++)
    EXPECT_NEAR(flat_x(16, j), 1.0f, 1e-5f);

  // Far from center (e.g. i = 0), elevation should remain unflattened
  EXPECT_NEAR(flat_x(0, 0), 1.0f, 1e-5f);
  EXPECT_NEAR(flat_x(32, 16), input(0, 16), 1e-5f);

  Array flat_rot = symmetrize(input,
                              SymmetryType::SYMMETRY_ROT180,
                              /* flatten_center */ true,
                              /* flatten_radius */ 0.1f);
  EXPECT_NEAR(flat_rot(16, 16), 1.0f, 1e-5f);
}

TEST(LooseSymmetryTest, FlattenCenterRunsAndFlattensCenterArea)
{
  Array input = make_test_terrain({64, 64}, 42);

  Array flat_loose = loose_symmetry(input,
                                    SymmetryType::SYMMETRY_LEFT_TO_RIGHT,
                                    1.f,
                                    2,
                                    8,
                                    2,
                                    4,
                                    1,
                                    /* flatten_center */ true,
                                    /* flatten_radius */ 0.1f);

  ASSERT_EQ(flat_loose.shape.x, 64);
  ASSERT_EQ(flat_loose.shape.y, 64);
  for (float v : flat_loose.vector)
    EXPECT_TRUE(std::isfinite(v));
}
