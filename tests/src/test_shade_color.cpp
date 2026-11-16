#include "highmap/colorize.hpp"
#include "highmap/internal/assert.hpp"
#include "highmap/local_metrics.hpp"
#include "highmap/texture.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(ShadeColor, TextureBasicAndAlphaPreservation)
{
  glm::ivec2 shape = {4, 4};
  Array      r(shape, 0.5f);
  Array      g(shape, 0.4f);
  Array      b(shape, 0.3f);
  Array      a(shape, 0.9f);

  Texture tex(r, g, b, a);
  EXPECT_EQ(tex.num_channels(), 4);

  Array factor(shape, 0.5f);
  factor(0, 0) = 0.0f;
  factor(1, 1) = 1.0f;

  shade_color(tex, factor);

  // verify pixel at (0, 0)
  glm::vec3 expected_00 = shade_color(glm::vec3(0.5f, 0.4f, 0.3f), 0.0f);
  EXPECT_NEAR(tex(0, 0, 0), expected_00.r, 1e-6f);
  EXPECT_NEAR(tex(0, 0, 1), expected_00.g, 1e-6f);
  EXPECT_NEAR(tex(0, 0, 2), expected_00.b, 1e-6f);

  // verify pixel at (1, 1)
  glm::vec3 expected_11 = shade_color(glm::vec3(0.5f, 0.4f, 0.3f), 1.0f);
  EXPECT_NEAR(tex(1, 1, 0), expected_11.r, 1e-6f);
  EXPECT_NEAR(tex(1, 1, 1), expected_11.g, 1e-6f);
  EXPECT_NEAR(tex(1, 1, 2), expected_11.b, 1e-6f);

  // verify pixel at (0, 1) with factor 0.5
  glm::vec3 expected_01 = shade_color(glm::vec3(0.5f, 0.4f, 0.3f), 0.5f);
  EXPECT_NEAR(tex(0, 1, 0), expected_01.r, 1e-6f);
  EXPECT_NEAR(tex(0, 1, 1), expected_01.g, 1e-6f);
  EXPECT_NEAR(tex(0, 1, 2), expected_01.b, 1e-6f);

  // verify alpha channel was preserved
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      EXPECT_EQ(tex(i, j, 3), 0.9f);
}

TEST(ShadeColor, EmptyOrMismatchedInputs)
{
  Texture empty_tex;
  Array   factor({4, 4}, 0.5f);

  // Should gracefully return without throwing or crashing
  shade_color(empty_tex, factor);
  EXPECT_TRUE(empty_tex.is_empty());

  Texture valid_tex({4, 4}, 3, 0.5f);
  Array   empty_factor;
  shade_color(valid_tex, empty_factor);
  EXPECT_EQ(valid_tex(0, 0, 0), 0.5f);

  Array mismatched_factor({2, 2}, 0.5f);
  shade_color(valid_tex, mismatched_factor);
  EXPECT_EQ(valid_tex(0, 0, 0), 0.5f);

  Texture one_channel_tex({4, 4}, 1, 0.5f);
  shade_color(one_channel_tex, factor);
  EXPECT_EQ(one_channel_tex(0, 0, 0), 0.5f);
}

TEST(ShadeTopographicResidual, FlatTerrainProducesNeutralShade)
{
  glm::ivec2 shape = {16, 16};
  Array      z(shape, 100.0f); // completely flat terrain

  Texture tex(shape, 3, 0.6f);

  shade_topographic_residual(tex, z, 3, 1.0f);

  // For flat terrain, residual is 0 everywhere, so factor is 0.5 + 0.5 * 0 =
  // 0.5
  glm::vec3 expected = shade_color(glm::vec3(0.6f, 0.6f, 0.6f), 0.5f);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_NEAR(tex(i, j, 0), expected.r, 1e-5f);
      EXPECT_NEAR(tex(i, j, 1), expected.g, 1e-5f);
      EXPECT_NEAR(tex(i, j, 2), expected.b, 1e-5f);
    }
}

TEST(ShadeTopographicResidual, PeakAndValley)
{
  glm::ivec2 shape = {17, 17};
  Array      z(shape, 0.0f);

  // Create a peak in the center
  z(8, 8) = 10.0f;

  Texture tex(shape, 3, 0.5f);
  shade_topographic_residual(tex, z, 3, 1.0f);

  // Peak has positive residual, so factor > 0.5
  float     center_r = tex(8, 8, 0);
  glm::vec3 neutral = shade_color(glm::vec3(0.5f, 0.5f, 0.5f), 0.5f);

  // With t > 0.5, color is lighter
  EXPECT_GT(center_r, neutral.r);
}

TEST(ShadeTopographicResidual, InvalidInputs)
{
  Texture empty_tex;
  Array   z({8, 8}, 1.0f);

  shade_topographic_residual(empty_tex, z, 2);
  EXPECT_TRUE(empty_tex.is_empty());

  Texture valid_tex({8, 8}, 3, 0.5f);
  Array   empty_z;
  shade_topographic_residual(valid_tex, empty_z, 2);
  EXPECT_EQ(valid_tex(0, 0, 0), 0.5f);

  Array mismatched_z({4, 4}, 1.0f);
  shade_topographic_residual(valid_tex, mismatched_z, 2);
  EXPECT_EQ(valid_tex(0, 0, 0), 0.5f);
}
