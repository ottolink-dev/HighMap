#include "highmap/array.hpp"
#include "highmap/internal/assert.hpp"
#include "highmap/primitives.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(ArrayIO, ConstructorExrNoRemapFlag)
{
  Array a = white({32, 16}, 0.f, 1.f, 0);

  const std::string fname = "test_array.exr";
  a.to_exr(fname);

  Array loaded(fname, false);

  EXPECT_EQ(loaded.shape.x, a.shape.x);
  EXPECT_EQ(loaded.shape.y, a.shape.y);

  EXPECT_TRUE(assert_almost_equal(a, loaded, 1e-4f));
}

TEST(ArrayIO, ConstructorNonExrRemapEnabled)
{
  Array a = white({32, 16}, 0.f, 1.f, 0);

  const std::string fname = "test_array.png";
  a.to_png_grayscale(fname, CV_32F);

  Array loaded(fname, true);

  EXPECT_EQ(loaded.shape.x, a.shape.x);
  EXPECT_EQ(loaded.shape.y, a.shape.y);

  // values may be normalized, so only structure check is strict
  EXPECT_EQ(loaded.vector.size(), a.vector.size());
}

TEST(ArrayIO, ConstructorTiffRoundTrip)
{
  Array a = white({32, 16}, 0.f, 1.f, 0);

  const std::string fname = "test_array.tiff";
  a.to_tiff(fname);

  Array loaded(fname, true);

  EXPECT_EQ(loaded.shape.x, a.shape.x);
  EXPECT_EQ(loaded.shape.y, a.shape.y);

  // TIFF is remapped -> only structural sanity
  EXPECT_EQ(loaded.vector.size(), a.vector.size());
}

TEST(ArrayIO, ExrVsPngDynamicRange)
{
  Array a(glm::ivec2(8, 8));
  for (int i = 0; i < 64; ++i)
    a.vector[i] = float(i) / 63.f;

  const std::string exr = "dyn.exr";
  const std::string png = "dyn.png";

  a.to_exr(exr);
  a.to_png_grayscale(png, CV_32F);

  Array exr_loaded(exr, false);
  Array png_loaded(png, true);

  float exr_range = exr_loaded.max() - exr_loaded.min();
  float png_range = png_loaded.max() - png_loaded.min();

  EXPECT_GT(exr_range, 0.5f);
  EXPECT_LE(png_range, 1.0f);
}

TEST(ArrayIO, ConstantArrayConstructor)
{
  Array a(glm::ivec2(8, 8));
  std::fill(a.vector.begin(), a.vector.end(), 3.14f);

  const std::string fname = "const.exr";
  a.to_exr(fname);

  Array loaded(fname, false);

  EXPECT_TRUE(assert_almost_equal(a, loaded, 1e-3f));
}
