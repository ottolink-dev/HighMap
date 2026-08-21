#include "highmap/authoring.hpp"
#include "highmap/dbg/assert.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(ElevationFromDistanceFields, ConstraintsPreservedExactly)
{
  glm::ivec2 shape = {64, 64};

  // Mountain ridge at center horizontal line
  Array mountains(shape, 0.f);
  for (int i = 0; i < shape.x; ++i)
  {
    mountains(i, 32) = 1.f;
  }

  // Coastline at row 16 and 48
  Array coastline(shape, 0.f);
  for (int i = 0; i < shape.x; ++i)
  {
    coastline(i, 16) = 1.f;
    coastline(i, 48) = 1.f;
  }

  // Boundary at top and bottom edges (row 0 and 63)
  Array boundary(shape, 0.f);
  for (int i = 0; i < shape.x; ++i)
  {
    boundary(i, 0) = 1.f;
    boundary(i, shape.y - 1) = 1.f;
  }

  // With smoothing_radius = 0.f (exact singular check)
  Array z = elevation_from_distance_fields(
      mountains,
      &boundary,
      &coastline,
      1.0f,
      0.0f, // smoothing_radius = 0 for exact checks
      1.f,
      -1.f,
      0.f);

  // Verify exact constraint values
  for (int i = 0; i < shape.x; ++i)
  {
    EXPECT_FLOAT_EQ(z(i, 32), 1.0f);
    EXPECT_FLOAT_EQ(z(i, 16), 0.0f);
    EXPECT_FLOAT_EQ(z(i, 48), 0.0f);
    EXPECT_FLOAT_EQ(z(i, 0), -1.0f);
    EXPECT_FLOAT_EQ(z(i, shape.y - 1), -1.0f);
  }

  // Verify intermediate values: between mountains (+1) and coast (0),
  // elevations should be > 0 and < 1
  for (int j = 17; j < 32; ++j)
  {
    EXPECT_GT(z(32, j), 0.0f);
    EXPECT_LT(z(32, j), 1.0f);
  }

  // Verify intermediate values: between coast (0) and boundary (-1), elevations
  // should be < 0 and > -1
  for (int j = 1; j < 16; ++j)
  {
    EXPECT_LT(z(32, j), 0.0f);
    EXPECT_GT(z(32, j), -1.0f);
  }
}

TEST(ElevationFromDistanceFields, SoftenedHarmonicWeightingSmoothTransitions)
{
  glm::ivec2 shape = {64, 64};

  Array mountains(shape, 0.f);
  mountains(32, 32) = 1.f;

  Array coastline(shape, 0.f);
  for (int i = 0; i < shape.x; ++i)
  {
    coastline(i, 16) = 1.f;
  }

  // Generate with default smoothing_radius = 2.0f
  Array z = elevation_from_distance_fields(mountains,
                                           nullptr,
                                           &coastline,
                                           1.0f,
                                           2.0f);

  // Ensure bounded within [-1, 1]
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_GE(z(i, j), -1.0f);
      EXPECT_LE(z(i, j), 1.0f);
    }
}

TEST(ElevationFromDistanceFields, TwoConstraintsHarmonicMean)
{
  glm::ivec2 shape = {65, 65};

  // Mountain at (32, 32)
  Array mountains(shape, 0.f);
  mountains(32, 32) = 1.f;

  // Default boundary around perimeter
  Array z = elevation_from_distance_fields(mountains,
                                           nullptr,
                                           nullptr,
                                           1.0f,
                                           0.0f,
                                           1.f,
                                           -1.f,
                                           0.f);

  // At center mountain point
  EXPECT_FLOAT_EQ(z(32, 32), 1.0f);

  // At corners (boundary)
  EXPECT_FLOAT_EQ(z(0, 0), -1.0f);
  EXPECT_FLOAT_EQ(z(64, 64), -1.0f);

  // Monotonically decreases away from the mountain towards the boundary
  for (int i = 32; i < 64; ++i)
  {
    EXPECT_GE(z(i, 32), z(i + 1, 32));
  }
}

TEST(ElevationFromDistanceFields, ExponentEffect)
{
  glm::ivec2 shape = {65, 65};

  Array mountains(shape, 0.f);
  mountains(32, 32) = 1.f;

  Array z_p05 = elevation_from_distance_fields(mountains,
                                               nullptr,
                                               nullptr,
                                               0.5f,
                                               0.0f);
  Array z_p20 = elevation_from_distance_fields(mountains,
                                               nullptr,
                                               nullptr,
                                               2.0f,
                                               0.0f);

  EXPECT_FLOAT_EQ(z_p05(32, 32), 1.0f);
  EXPECT_FLOAT_EQ(z_p20(32, 32), 1.0f);
  EXPECT_FLOAT_EQ(z_p05(0, 0), -1.0f);
  EXPECT_FLOAT_EQ(z_p20(0, 0), -1.0f);
}
