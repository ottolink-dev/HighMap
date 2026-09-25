#include "highmap/boundary.hpp"
#include "highmap/primitives.hpp"

#include <gtest/gtest.h>

TEST(BoundaryTest, DomainBoundaryEnumValues)
{
  EXPECT_EQ(hmap::DomainBoundary::BOUNDARY_LEFT, 0);
  EXPECT_EQ(hmap::DomainBoundary::BOUNDARY_RIGHT, 1);
  EXPECT_EQ(hmap::DomainBoundary::BOUNDARY_TOP, 2);
  EXPECT_EQ(hmap::DomainBoundary::BOUNDARY_BOTTOM, 3);
}

TEST(BoundaryTest, PickBoundaryCellFlatArray)
{
  hmap::Array z(glm::ivec2(16, 16), 5.0f);

  glm::ivec2 cell_b = hmap::pick_boundary_cell(
      z,
      hmap::DomainBoundary::BOUNDARY_BOTTOM,
      42);
  EXPECT_EQ(cell_b.y, 0);
  EXPECT_GE(cell_b.x, 0);
  EXPECT_LT(cell_b.x, 16);

  glm::ivec2 cell_t = hmap::pick_boundary_cell(
      z,
      hmap::DomainBoundary::BOUNDARY_TOP,
      42);
  EXPECT_EQ(cell_t.y, 15);
  EXPECT_GE(cell_t.x, 0);
  EXPECT_LT(cell_t.x, 16);

  glm::ivec2 cell_l = hmap::pick_boundary_cell(
      z,
      hmap::DomainBoundary::BOUNDARY_LEFT,
      42);
  EXPECT_EQ(cell_l.x, 0);
  EXPECT_GE(cell_l.y, 0);
  EXPECT_LT(cell_l.y, 16);

  glm::ivec2 cell_r = hmap::pick_boundary_cell(
      z,
      hmap::DomainBoundary::BOUNDARY_RIGHT,
      42);
  EXPECT_EQ(cell_r.x, 15);
  EXPECT_GE(cell_r.y, 0);
  EXPECT_LT(cell_r.y, 16);
}

TEST(BoundaryTest, GenerateBufferedArrayPivotCentered)
{
  hmap::Array arr(glm::ivec2(5, 5));
  for (int j = 0; j < 5; ++j)
    for (int i = 0; i < 5; ++i)
      arr(i, j) = float(i + j * 10);

  glm::ivec4  buffers(2, 2, 2, 2);
  hmap::Array buffered = hmap::generate_buffered_array(arr, buffers, false);

  EXPECT_EQ(buffered.shape.x, 9);
  EXPECT_EQ(buffered.shape.y, 9);

  // Interior is in [2, 6] (West pivot is 2, East pivot is 6)
  // For West buffer:
  // cell 1 (distance 1 to pivot 2) mirrors cell 3:
  EXPECT_FLOAT_EQ(buffered(1, 3), buffered(3, 3));
  // cell 0 (distance 2 to pivot 2) mirrors cell 4:
  EXPECT_FLOAT_EQ(buffered(0, 3), buffered(4, 3));

  // For East buffer:
  // cell 7 (distance 1 to pivot 6) mirrors cell 5:
  EXPECT_FLOAT_EQ(buffered(7, 3), buffered(5, 3));
  // cell 8 (distance 2 to pivot 6) mirrors cell 4:
  EXPECT_FLOAT_EQ(buffered(8, 3), buffered(4, 3));
}

TEST(BoundaryTest, FalloffWithAndWithoutNoise)
{
  hmap::Array a1(glm::ivec2(16, 16), 1.f);
  hmap::falloff(a1, 1.f);

  hmap::Array noise(glm::ivec2(16, 16), 0.f);
  hmap::Array a2(glm::ivec2(16, 16), 1.f);
  hmap::falloff(a2, 1.f, hmap::DistanceFunction::EUCLIDIAN, &noise);

  for (int j = 0; j < 16; ++j)
    for (int i = 0; i < 16; ++i)
      EXPECT_NEAR(a1(i, j), a2(i, j), 1e-5f);
}

TEST(BoundaryTest, MakePeriodicStitchingContinuity)
{
  hmap::Array z = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                  glm::ivec2(64, 64),
                                  glm::vec2(4.f, 4.f),
                                  1);
  hmap::Array zp = hmap::make_periodic_stitching(z, 0.5f);

  EXPECT_EQ(zp.shape.x, 64);
  EXPECT_EQ(zp.shape.y, 64);

  // Check continuity across boundary
  for (int j = 0; j < zp.shape.y; ++j)
  {
    EXPECT_NEAR(zp(0, j), zp(zp.shape.x - 1, j), 1e-2f);
  }
  for (int i = 0; i < zp.shape.x; ++i)
  {
    EXPECT_NEAR(zp(i, 0), zp(i, zp.shape.y - 1), 1e-2f);
  }
}

TEST(BoundaryTest, SetBordersBboxFullDomain)
{
  hmap::Array a(glm::ivec2(65, 65), 1.0f);
  glm::vec4   border_values = {0.0f, 2.0f, -1.0f, 3.0f};
  glm::vec4   buffer_sizes = {0.2f, 0.2f, 0.2f, 0.2f};
  glm::vec4   bbox = {0.0f, 1.0f, 0.0f, 1.0f};

  hmap::set_borders(a, border_values, buffer_sizes, bbox);

  // West boundary (i = 0, x = 0)
  EXPECT_NEAR(a(0, 32), 0.0f, 1e-5f);
  // East boundary (i = 64, x = 64/65)
  EXPECT_NEAR(a(64, 32), 2.0f, 0.05f);
  // South boundary (j = 0, y = 0)
  EXPECT_NEAR(a(32, 0), -1.0f, 1e-5f);
  // North boundary (j = 64, y = 64/65)
  EXPECT_NEAR(a(32, 64), 3.0f, 0.05f);
  // Center (x = 0.5, y = 0.5) should be untouched
  EXPECT_NEAR(a(32, 32), 1.0f, 1e-5f);
}

TEST(BoundaryTest, SetBordersBboxSubdomain)
{
  // Top-right tile: x in [0.5, 1.0], y in [0.5, 1.0]
  hmap::Array a_tr(glm::ivec2(33, 33), 10.0f);
  glm::vec4   border_values = {0.0f, 2.0f, -1.0f, 3.0f};
  glm::vec4   buffer_sizes = {0.2f, 0.2f, 0.2f, 0.2f};
  glm::vec4   bbox_tr = {0.5f, 1.0f, 0.5f, 1.0f};

  hmap::set_borders(a_tr, border_values, buffer_sizes, bbox_tr);

  // Left side of tile (x = 0.5) is far from west boundary (x = 0), so untouched
  EXPECT_NEAR(a_tr(0, 16), 10.0f, 1e-5f);
  // Bottom side of tile (y = 0.5) is far from south boundary (y = 0), so
  // untouched
  EXPECT_NEAR(a_tr(16, 0), 10.0f, 1e-5f);
  // Right side of tile (x = 65/66) touches east boundary
  EXPECT_NEAR(a_tr(32, 16), 2.0f, 0.2f);
  // Top side of tile (y = 65/66) touches north boundary
  EXPECT_NEAR(a_tr(16, 32), 3.0f, 0.2f);

  // Interior tile: x in [0.3, 0.7], y in [0.3, 0.7] (all edges > 0.2 away from
  // domain boundaries)
  hmap::Array a_mid(glm::ivec2(33, 33), 5.0f);
  glm::vec4   bbox_mid = {0.3f, 0.7f, 0.3f, 0.7f};
  hmap::set_borders(a_mid, border_values, buffer_sizes, bbox_mid);

  for (int j = 0; j < 33; ++j)
    for (int i = 0; i < 33; ++i)
      EXPECT_NEAR(a_mid(i, j), 5.0f, 1e-5f);

  // Test uniform float overload
  hmap::Array a_uniform(glm::ivec2(65, 65), 1.0f);
  glm::vec4   bbox = {0.0f, 1.0f, 0.0f, 1.0f};
  hmap::set_borders(a_uniform, 0.0f, 0.2f, bbox);
  EXPECT_NEAR(a_uniform(0, 32), 0.0f, 1e-5f);
  EXPECT_NEAR(a_uniform(64, 32), 0.0f, 0.05f);
  EXPECT_NEAR(a_uniform(32, 0), 0.0f, 1e-5f);
  EXPECT_NEAR(a_uniform(32, 64), 0.0f, 0.05f);
  EXPECT_NEAR(a_uniform(32, 32), 1.0f, 1e-5f);
}
