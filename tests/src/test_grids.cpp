#include "highmap.hpp"

#include <gtest/gtest.h>

TEST(Grids, ConvertLengthToPixel)
{
  EXPECT_EQ(hmap::convert_length_to_pixel(0.5f, 100, true, false, 1.f), 50);
  EXPECT_EQ(hmap::convert_length_to_pixel(0.0f, 100, true, false, 1.f), 1);
  EXPECT_EQ(hmap::convert_length_to_pixel(0.0f, 100, false, false, 1.f), 0);
  EXPECT_EQ(hmap::convert_length_to_pixel(1.5f, 100, false, true, 1.f), 99);
}

TEST(Grids, GridXyVector)
{
  std::vector<float> x, y;
  hmap::grid_xy_vector(x, y, {4, 3}, {0.f, 1.f, 0.f, 1.f}, true);

  ASSERT_EQ(x.size(), 4u);
  ASSERT_EQ(y.size(), 3u);

  EXPECT_NEAR(x[0], 0.f, 1e-5f);
  EXPECT_NEAR(x[3], 1.f, 1e-5f);
  EXPECT_NEAR(y[0], 0.f, 1e-5f);
  EXPECT_NEAR(y[2], 1.f, 1e-5f);
}

TEST(Grids, RescaleGridFromUnitSquareToBbox)
{
  std::vector<float> x = {0.f, 0.5f, 1.f};
  std::vector<float> y = {0.f, 1.f};

  hmap::rescale_grid_from_unit_square_to_bbox(x, y, {10.f, 20.f, 5.f, 15.f});

  EXPECT_NEAR(x[0], 10.f, 1e-5f);
  EXPECT_NEAR(x[1], 15.f, 1e-5f);
  EXPECT_NEAR(x[2], 20.f, 1e-5f);

  EXPECT_NEAR(y[0], 5.f, 1e-5f);
  EXPECT_NEAR(y[1], 15.f, 1e-5f);

  // single element edge cases (no div by zero)
  std::vector<float> single_x = {0.f};
  std::vector<float> single_y = {0.f};
  hmap::rescale_grid_from_unit_square_to_bbox(single_x,
                                              single_y,
                                              {10.f, 20.f, 5.f, 15.f});
  EXPECT_NEAR(single_x[0], 10.f, 1e-5f);
  EXPECT_NEAR(single_y[0], 5.f, 1e-5f);
}

TEST(Grids, RescalePointsToUnitSquare)
{
  std::vector<float> x = {10.f, 15.f, 20.f};
  std::vector<float> y = {5.f, 15.f};

  hmap::rescale_points_to_unit_square(x, y, {10.f, 20.f, 5.f, 15.f});

  EXPECT_NEAR(x[0], 0.f, 1e-5f);
  EXPECT_NEAR(x[1], 0.5f, 1e-5f);
  EXPECT_NEAR(x[2], 1.f, 1e-5f);

  EXPECT_NEAR(y[0], 0.f, 1e-5f);
  EXPECT_NEAR(y[1], 1.f, 1e-5f);

  // degenerate bounding box (no div by zero)
  std::vector<float> dx = {10.f};
  std::vector<float> dy = {5.f};
  hmap::rescale_points_to_unit_square(dx, dy, {10.f, 10.f, 5.f, 5.f});
  EXPECT_FALSE(std::isnan(dx[0]));
  EXPECT_FALSE(std::isnan(dy[0]));
}

TEST(Grids, NeighborhoodStructures)
{
  using namespace hmap::neighborhood;

  // VON_NEUMANN_4
  EXPECT_EQ(VON_NEUMANN_4.size(), 4u);
  for (const auto &nbr : VON_NEUMANN_4)
  {
    EXPECT_NEAR(nbr.distance * nbr.inv_distance, 1.0f, 1e-5f);
    EXPECT_EQ(std::abs(nbr.offset.x) + std::abs(nbr.offset.y), 1);
  }

  // MOORE_8
  EXPECT_EQ(MOORE_8.size(), 8u);
  for (size_t i = 0; i < 4; ++i)
  {
    EXPECT_NEAR(MOORE_8[i].distance, 1.0f, 1e-5f);
    EXPECT_NEAR(MOORE_8[i].inv_distance, 1.0f, 1e-5f);
  }
  for (size_t i = 4; i < 8; ++i)
  {
    EXPECT_NEAR(MOORE_8[i].distance, std::numbers::sqrt2_v<float>, 1e-5f);
    EXPECT_NEAR(MOORE_8[i].inv_distance,
                1.0f / std::numbers::sqrt2_v<float>,
                1e-5f);
  }

  // D8_FLOW
  EXPECT_EQ(D8_FLOW.size(), 8u);
  for (const auto &nbr : D8_FLOW)
  {
    EXPECT_NEAR(nbr.distance * nbr.inv_distance, 1.0f, 1e-5f);
  }

  // DINF_8
  EXPECT_EQ(DINF_8.size(), 8u);
  for (size_t i = 0; i < 4; ++i)
  {
    EXPECT_NEAR(DINF_8[i].weight, 0.500f, 1e-5f);
  }
  for (size_t i = 4; i < 8; ++i)
  {
    EXPECT_NEAR(DINF_8[i].weight, 0.354f, 1e-5f);
  }

  // GEOMORPHON_8
  EXPECT_EQ(GEOMORPHON_8.size(), 8u);
  for (const auto &nbr : GEOMORPHON_8)
  {
    EXPECT_NEAR(nbr.distance * nbr.inv_distance, 1.0f, 1e-5f);
  }
}
