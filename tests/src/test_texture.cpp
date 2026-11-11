#include "highmap/blending.hpp"
#include "highmap/colorize.hpp"
#include "highmap/colormaps.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/primitives.hpp"
#include "highmap/texture.hpp"
#include "highmap/transform.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(TextureIO, PngRoundTripOrientationAndData)
{
  glm::ivec2 shape = {16, 32};

  // create a 3-channel texture with a distinct coordinate/color pattern to
  // check orientation
  Texture tex(shape, 3);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      // Red channel depends on X, Green on Y, Blue on diagonal
      float r = float(i) / float(shape.x - 1);
      float g = float(j) / float(shape.y - 1);
      float b = float(i + j) / float(shape.x + shape.y - 2);
      tex.set_pixel(i, j, glm::vec3(r, g, b));
    }

  const std::string fname = "test_texture_roundtrip.png";

  // to_png flips vertically (cv::flip(mat, mat, 0)) to match typical graphics
  // conventions
  tex.to_png(fname, CV_8U);

  // since to_png flips y, we construct the loaded texture with flip_j = true to
  // invert it back
  Texture loaded(fname, true);

  EXPECT_EQ(loaded.shape.x, tex.shape.x);
  EXPECT_EQ(loaded.shape.y, tex.shape.y);

  // file load constructor creates 4 channels (RGBA)
  EXPECT_EQ(loaded.num_channels(), 4);

  // verify that the data and orientation match within CV_8U quantization error
  // (1.0 / 255.0)
  float tolerance = 1.5f / 255.f;

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      glm::vec3 orig = tex.get_pixel3(i, j);
      glm::vec4 read_val = loaded.get_pixel4(i, j);

      EXPECT_NEAR(orig.x, read_val.x, tolerance);
      EXPECT_NEAR(orig.y, read_val.y, tolerance);
      EXPECT_NEAR(orig.z, read_val.z, tolerance);
      EXPECT_NEAR(read_val.w,
                  1.0f,
                  tolerance); // Alpha should be filled with 1.f
    }
}

TEST(TextureColorize, CustomColormapAndOperations)
{
  glm::ivec2 shape = {8, 8};
  Array      level = white(shape, 0.f, 1.f, 0);

  // Red, Green, Blue at positions 0.0, 0.5, 1.0
  std::vector<float>     positions = {0.f, 0.5f, 1.f};
  std::vector<glm::vec3> colors = {glm::vec3(1.f, 0.f, 0.f),
                                   glm::vec3(0.f, 1.f, 0.f),
                                   glm::vec3(0.f, 0.f, 1.f)};

  Texture col = colorize(level, 0.f, 1.f, positions, colors, false, nullptr);
  EXPECT_EQ(col.num_channels(), 3);
  EXPECT_EQ(col.shape, shape);

  // Check luminance
  Array lum = luminance(col);
  EXPECT_EQ(lum.shape, shape);
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      glm::vec3 rgb = col.get_pixel3(i, j);
      float     expected_lum = 0.299f * rgb.x + 0.587f * rgb.y + 0.114f * rgb.z;
      EXPECT_NEAR(lum(i, j), expected_lum, 1e-5f);
    }
  }

  // Check mixing
  Texture t1(shape, 4, 0.f);
  Texture t2(shape, 4, 1.f);
  // Set alpha channels
  t1[3] = Array(shape, 1.f);
  t2[3] = Array(shape, 1.f);

  Texture t_mixed = mix(t1, t2, MixMethod::MM_LINEAR);
  EXPECT_EQ(t_mixed.num_channels(), 4);
  // Mixing equal alphas (1.f, 1.f) with mixing factor t = 1.0 / (1.0 + 1.0 *
  // 0.0) = 1.0 Result should equal t2
  for (int c = 0; c < 3; ++c)
  {
    for (int idx = 0; idx < shape.x * shape.y; ++idx)
    {
      EXPECT_NEAR(t_mixed[c].vector[idx], 1.f, 1e-5f);
    }
  }

  // Also test Mixbox mixing compiles and runs
  Texture t_mixbox = mix(t1, t2, MixMethod::MM_MIXBOX);
  EXPECT_EQ(t_mixbox.num_channels(), 4);

  // Check mix_normal_map basic execution
  Texture n1(shape, 4, 0.5f); // flat normals
  Texture n2(shape, 4, 0.5f);
  n1[2] = Array(shape, 1.f);
  n2[2] = Array(shape, 1.f);

  Texture blended = mix_normal_map(n1,
                                   n2,
                                   1.0f,
                                   NormalMapBlendingMethod::NMAP_LINEAR);
  EXPECT_EQ(blended.num_channels(), 4);
}

TEST(TextureColorize, MixWithMask)
{
  glm::ivec2 shape = {4, 4};
  Texture    t1(shape, 3, 0.0f); // all black
  Texture    t2(shape, 3, 1.0f); // all white

  // 1. Mask = 0.0 -> result should be t1 (0.0)
  Array   mask_0(shape, 0.0f);
  Texture mixed_0 = mix(t1, t2, mask_0, MixMethod::MM_LINEAR);
  EXPECT_EQ(mixed_0.num_channels(), 3);
  for (int c = 0; c < 3; ++c)
    for (int idx = 0; idx < shape.x * shape.y; ++idx)
      EXPECT_NEAR(mixed_0[c].vector[idx], 0.0f, 1e-5f);

  // 2. Mask = 1.0 -> result should be t2 (1.0)
  Array   mask_1(shape, 1.0f);
  Texture mixed_1 = mix(t1, t2, mask_1, MixMethod::MM_LINEAR);
  EXPECT_EQ(mixed_1.num_channels(), 3);
  for (int c = 0; c < 3; ++c)
    for (int idx = 0; idx < shape.x * shape.y; ++idx)
      EXPECT_NEAR(mixed_1[c].vector[idx], 1.0f, 1e-5f);

  // 3. Mask = 0.5 with MM_LINEAR -> result should be 0.5
  Array   mask_half(shape, 0.5f);
  Texture mixed_half = mix(t1, t2, mask_half, MixMethod::MM_LINEAR);
  EXPECT_EQ(mixed_half.num_channels(), 3);
  for (int c = 0; c < 3; ++c)
    for (int idx = 0; idx < shape.x * shape.y; ++idx)
      EXPECT_NEAR(mixed_half[c].vector[idx], 0.5f, 1e-5f);

  // 4. MM_SQRT_AVG with 0.5
  Texture mixed_sqrt = mix(t1, t2, mask_half, MixMethod::MM_SQRT_AVG);
  EXPECT_EQ(mixed_sqrt.num_channels(), 3);
  float expected_sqrt = std::sqrt(0.5f * 0.f * 0.f + 0.5f * 1.f * 1.f);
  for (int c = 0; c < 3; ++c)
    for (int idx = 0; idx < shape.x * shape.y; ++idx)
      EXPECT_NEAR(mixed_sqrt[c].vector[idx], expected_sqrt, 1e-5f);

  // 5. MM_MIXBOX
  Texture mixed_mixbox = mix(t1, t2, mask_half, MixMethod::MM_MIXBOX);
  EXPECT_EQ(mixed_mixbox.num_channels(), 3);

  // 6. Test with gain on mask
  // With mask = 0.25 and gain = 2.0:
  // gain(0.25, 2.0) = 0.5 * (2 * 0.25)^2 = 0.5 * 0.25 = 0.125
  Array   mask_quarter(shape, 0.25f);
  Texture mixed_gain = mix(t1, t2, mask_quarter, MixMethod::MM_LINEAR, 2.0f);
  for (int c = 0; c < 3; ++c)
    for (int idx = 0; idx < shape.x * shape.y; ++idx)
      EXPECT_NEAR(mixed_gain[c].vector[idx], 0.125f, 1e-5f);

  // 7. 4-channel textures
  Texture t1_rgba(shape, 4, 0.2f);
  Texture t2_rgba(shape, 4, 0.8f);
  t1_rgba[3] = Array(shape, 0.4f);
  t2_rgba[3] = Array(shape, 0.9f);
  Texture mixed_rgba = mix(t1_rgba, t2_rgba, mask_half, MixMethod::MM_LINEAR);
  EXPECT_EQ(mixed_rgba.num_channels(), 4);
  for (int c = 0; c < 3; ++c)
    for (int idx = 0; idx < shape.x * shape.y; ++idx)
      EXPECT_NEAR(mixed_rgba[c].vector[idx], 0.5f, 1e-5f);
  for (int idx = 0; idx < shape.x * shape.y; ++idx)
    EXPECT_NEAR(mixed_rgba[3].vector[idx], 0.65f, 1e-5f);

  // 8. Invalid inputs
  Texture empty_tex;
  EXPECT_EQ(mix(empty_tex, t2, mask_half).num_channels(), 0);
  EXPECT_EQ(mix(t1, empty_tex, mask_half).num_channels(), 0);
  Array empty_mask;
  EXPECT_EQ(mix(t1, t2, empty_mask).num_channels(), 0);
  Array mismatched_mask({8, 8}, 0.5f);
  EXPECT_EQ(mix(t1, t2, mismatched_mask).num_channels(), 0);
}

TEST(TextureColorize, BivariateReverseAndNoise)
{
  glm::ivec2 shape = {4, 4};

  Array a1(shape, 0.f);
  Array a2(shape, 0.f);

  Array noise1(shape, 0.f);
  Array noise2(shape, 0.f);

  // Move both normalized values to 1.0 through noise.
  noise1 = 1.f;
  noise2 = 1.f;

  const std::vector<float> positions = {0.f, 1.f};

  const std::vector<glm::vec3> colors1 = {glm::vec3(1.f, 0.f, 0.f),
                                          glm::vec3(0.f, 0.f, 1.f)};

  const std::vector<glm::vec3> colors2 = {glm::vec3(0.f, 1.f, 0.f),
                                          glm::vec3(1.f, 1.f, 1.f)};

  Texture result = colorize_bivariate(a1,
                                      a2,
                                      {0.f, 1.f},
                                      {0.f, 1.f},
                                      positions,
                                      positions,
                                      colors1,
                                      colors2,
                                      MixMethod::MM_LINEAR,
                                      true,
                                      true,
                                      &noise1,
                                      &noise2);

  // With noise, both values are 1.0. With reverse=true, both map to the
  // first color of their respective colormaps.
  // color1 = red, color2 = green -> linear mix = yellow.
  glm::vec3 color = result.get_pixel3(0, 0);

  EXPECT_NEAR(color.r, 0.5f, 1e-5f);
  EXPECT_NEAR(color.g, 0.5f, 1e-5f);
  EXPECT_NEAR(color.b, 0.f, 1e-5f);
}

TEST(TextureColorize, TerrainTextureBasic)
{
  glm::ivec2 shape = {16, 16};
  Array      z(shape, 10.f);
  z(0, 0) = 5.f;
  z(15, 15) = 15.f;

  TerrainTextureParams params;

  // Default computation: should automatically compute zmin and zmax
  Texture tex_default = terrain_texture(z, nullptr, nullptr, params);
  EXPECT_EQ(tex_default.shape, shape);
  EXPECT_EQ(tex_default.num_channels(), 4);

  // Invalid custom noise shape should fail validation and return empty Texture
  Array   wrong_noise(glm::ivec2(8, 8), 0.5f);
  Texture tex_invalid = terrain_texture(z,
                                        nullptr,
                                        nullptr,
                                        params,
                                        glm::vec4(0.f, 0.f, 1.f, 1.f),
                                        &wrong_noise);
  EXPECT_TRUE(tex_invalid.is_empty());

  // Custom noise bypassing procedural noise
  Array   custom_noise(shape, 0.5f);
  Texture tex_custom_noise = terrain_texture(z,
                                             nullptr,
                                             nullptr,
                                             params,
                                             glm::vec4(0.f, 0.f, 1.f, 1.f),
                                             &custom_noise);
  EXPECT_EQ(tex_custom_noise.shape, shape);
  EXPECT_EQ(tex_custom_noise.num_channels(), 4);

  // Deactivated color noise (color_noise <= 0)
  params.color_noise = 0.f;
  Texture tex_no_color_noise = terrain_texture(z, nullptr, nullptr, params);
  EXPECT_EQ(tex_no_color_noise.shape, shape);
  EXPECT_EQ(tex_no_color_noise.num_channels(), 4);
}

TEST(TextureColorize, TerrainTextureSnowDepth)
{
  glm::ivec2 shape = {32, 32};
  Array      z(shape, 0.f);

  TerrainTextureParams params;
  params.snow_amount = 1.0f;
  params.snow_line = 0.99f; // high snow line so elevation-based snow is absent
  params.color_noise = 0.f;

  // 1. Without snow depth and elevation below snow line -> no snow
  Texture tex_no_snow = terrain_texture(z, nullptr, nullptr, params);
  EXPECT_EQ(tex_no_snow.shape, shape);

  // 2. With snow depth array provided (unnormalized depth values, e.g. 5.0f)
  Array snow_depth(shape, 0.f);
  for (int j = 0; j < 16; ++j)
    for (int i = 0; i < 32; ++i)
      snow_depth(i, j) = 5.0f;

  Texture tex_snow = terrain_texture(z, nullptr, &snow_depth, params);
  EXPECT_EQ(tex_snow.shape, shape);

  // Snowy half should differ from non-snowy half and approach params.snow color
  glm::vec3 sample_snowy(tex_snow[0](16, 8),
                         tex_snow[1](16, 8),
                         tex_snow[2](16, 8));
  glm::vec3 sample_no_snow(tex_snow[0](16, 24),
                           tex_snow[1](16, 24),
                           tex_snow[2](16, 24));
  EXPECT_GT(sample_snowy.r, sample_no_snow.r);
  EXPECT_NEAR(sample_snowy.r, params.snow.r, 0.1f);

  // 3. Snow amount acts as overall multiplier (snow_amount = 0 eliminates snow)
  params.snow_amount = 0.f;
  Texture tex_zero_snow = terrain_texture(z, nullptr, &snow_depth, params);
  for (int c = 0; c < 3; ++c)
  {
    EXPECT_TRUE(assert_almost_equal(tex_zero_snow[c], tex_no_snow[c], 1e-4f));
  }
}

TEST(TextureTransform, BasicTransforms)
{
  glm::ivec2 shape = {4, 8};
  Texture    tex(shape, 3);

  // Initialize with distinct values
  for (int c = 0; c < 3; ++c)
  {
    for (int j = 0; j < shape.y; ++j)
    {
      for (int i = 0; i < shape.x; ++i)
      {
        tex[c](i, j) = static_cast<float>(c * 100 + j * 10 + i);
      }
    }
  }

  // 1. Test flip_lr
  {
    Texture temp = tex;
    flip_lr(temp);
    EXPECT_EQ(temp.shape, shape);
    for (int c = 0; c < 3; ++c)
    {
      for (int j = 0; j < shape.y; ++j)
      {
        for (int i = 0; i < shape.x; ++i)
        {
          EXPECT_FLOAT_EQ(temp[c](i, j), tex[c](shape.x - 1 - i, j));
        }
      }
    }
  }

  // 2. Test flip_ud
  {
    Texture temp = tex;
    flip_ud(temp);
    EXPECT_EQ(temp.shape, shape);
    for (int c = 0; c < 3; ++c)
    {
      for (int j = 0; j < shape.y; ++j)
      {
        for (int i = 0; i < shape.x; ++i)
        {
          EXPECT_FLOAT_EQ(temp[c](i, j), tex[c](i, shape.y - 1 - j));
        }
      }
    }
  }

  // 3. Test rot180
  {
    Texture temp = tex;
    rot180(temp);
    EXPECT_EQ(temp.shape, shape);
    for (int c = 0; c < 3; ++c)
    {
      for (int j = 0; j < shape.y; ++j)
      {
        for (int i = 0; i < shape.x; ++i)
        {
          EXPECT_FLOAT_EQ(temp[c](i, j),
                          tex[c](shape.x - 1 - i, shape.y - 1 - j));
        }
      }
    }
  }

  // 4. Test transpose
  {
    Texture    temp = transpose(tex);
    glm::ivec2 expected_shape = {shape.y, shape.x};
    EXPECT_EQ(temp.shape, expected_shape);
    for (int c = 0; c < 3; ++c)
    {
      for (int j = 0; j < shape.y; ++j)
      {
        for (int i = 0; i < shape.x; ++i)
        {
          EXPECT_FLOAT_EQ(temp[c](j, i), tex[c](i, j));
        }
      }
    }
  }

  // 5. Test rot90
  {
    Texture temp = tex;
    rot90(temp);
    glm::ivec2 expected_shape = {shape.y, shape.x};
    EXPECT_EQ(temp.shape, expected_shape);
    for (int c = 0; c < 3; ++c)
    {
      for (int j = 0; j < shape.y; ++j)
      {
        for (int i = 0; i < shape.x; ++i)
        {
          EXPECT_FLOAT_EQ(temp[c](j, shape.x - 1 - i), tex[c](i, j));
        }
      }
    }
  }

  // 6. Test rot270
  {
    Texture temp = tex;
    rot270(temp);
    glm::ivec2 expected_shape = {shape.y, shape.x};
    EXPECT_EQ(temp.shape, expected_shape);
    for (int c = 0; c < 3; ++c)
    {
      for (int j = 0; j < shape.y; ++j)
      {
        for (int i = 0; i < shape.x; ++i)
        {
          EXPECT_FLOAT_EQ(temp[c](shape.y - 1 - j, i), tex[c](i, j));
        }
      }
    }
  }
}

TEST(TextureBlending, BlendPoissonBf)
{
  glm::ivec2 shape = {32, 32};
  glm::vec2  kw = {2.f, 2.f};

  Array   z1_r = noise_fbm(NoiseType::PERLIN, shape, kw, 1);
  Array   z1_g = noise_fbm(NoiseType::PERLIN, shape, kw, 2);
  Array   z1_b = noise_fbm(NoiseType::PERLIN, shape, kw, 3);
  Texture tex1(z1_r, z1_g, z1_b);

  Array   z2_r = noise_fbm(NoiseType::WORLEY, shape, 2.f * kw, 4);
  Array   z2_g = noise_fbm(NoiseType::WORLEY, shape, 2.f * kw, 5);
  Array   z2_b = noise_fbm(NoiseType::WORLEY, shape, 2.f * kw, 6);
  Texture tex2(z2_r, z2_g, z2_b);

  const int iterations = 50;

  // Test global Poisson blending
  Texture blended = gpu::blend_poisson_bf(tex1, tex2, iterations);
  EXPECT_EQ(blended.shape, shape);
  EXPECT_EQ(blended.num_channels(), 3);

  for (int c = 0; c < 3; ++c)
  {
    Array expected_ch = gpu::blend_poisson_bf(tex1[c], tex2[c], iterations);
    EXPECT_TRUE(assert_almost_equal(blended[c], expected_ch, 1e-5f));
  }

  // Test masked Poisson blending
  Array   mask(shape, 0.5f);
  Texture blended_mask = gpu::blend_poisson_bf(tex1, tex2, iterations, &mask);
  EXPECT_EQ(blended_mask.shape, shape);
  EXPECT_EQ(blended_mask.num_channels(), 3);

  for (int c = 0; c < 3; ++c)
  {
    Array expected_ch = gpu::blend_poisson_bf(tex1[c],
                                              tex2[c],
                                              iterations,
                                              &mask);
    EXPECT_TRUE(assert_almost_equal(blended_mask[c], expected_ch, 1e-5f));
  }

  // Invalid / mismatched inputs
  Texture empty_tex;
  EXPECT_EQ(gpu::blend_poisson_bf(empty_tex, tex2, iterations).num_channels(),
            0);

  Texture mismatched_channels(shape, 2);
  EXPECT_EQ(gpu::blend_poisson_bf(tex1, mismatched_channels, iterations)
                .num_channels(),
            0);

  Texture mismatched_shape({16, 16}, 3);
  EXPECT_EQ(
      gpu::blend_poisson_bf(tex1, mismatched_shape, iterations).num_channels(),
      0);

  Array mismatched_mask({16, 16}, 1.f);
  EXPECT_EQ(gpu::blend_poisson_bf(tex1, tex2, iterations, &mismatched_mask)
                .num_channels(),
            0);
}
