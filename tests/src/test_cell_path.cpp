#include "highmap.hpp"

#include <gtest/gtest.h>

TEST(CellPath, ConstructionAndContainerInterface)
{
  hmap::CellPath empty_path;
  EXPECT_TRUE(empty_path.empty());
  EXPECT_EQ(empty_path.size(), 0u);

  std::vector<glm::ivec2> indices = {{0, 0}, {1, 2}, {3, 4}};
  hmap::CellPath          path(indices);

  EXPECT_FALSE(path.empty());
  EXPECT_EQ(path.size(), 3u);
  EXPECT_EQ(path.front(), glm::ivec2(0, 0));
  EXPECT_EQ(path.back(), glm::ivec2(3, 4));
  EXPECT_EQ(path[1], glm::ivec2(1, 2));

  path.push_back({5, 6});
  EXPECT_EQ(path.size(), 4u);
  EXPECT_EQ(path.back(), glm::ivec2(5, 6));

  path.clear();
  EXPECT_TRUE(path.empty());
}

TEST(CellPath, BresenhamLine)
{
  hmap::CellPath line;
  hmap::add_line_bresenham(line, {0, 0}, {3, 3});

  ASSERT_EQ(line.size(), 4u);
  EXPECT_EQ(line[0], glm::ivec2(0, 0));
  EXPECT_EQ(line[1], glm::ivec2(1, 1));
  EXPECT_EQ(line[2], glm::ivec2(2, 2));
  EXPECT_EQ(line[3], glm::ivec2(3, 3));
}

TEST(CellPath, AdjacencyCheckAndEnforcement)
{
  hmap::CellPath connected(
      std::vector<glm::ivec2>{{0, 0}, {1, 1}, {2, 1}, {2, 2}});
  EXPECT_TRUE(hmap::is_path_adjacent(connected));

  hmap::CellPath with_gaps(std::vector<glm::ivec2>{{0, 0}, {0, 3}, {4, 3}});
  EXPECT_FALSE(hmap::is_path_adjacent(with_gaps));

  hmap::enforce_path_adjacency(with_gaps);
  EXPECT_TRUE(hmap::is_path_adjacent(with_gaps));
  EXPECT_EQ(with_gaps.front(), glm::ivec2(0, 0));
  EXPECT_EQ(with_gaps.back(), glm::ivec2(4, 3));

  // ensure no consecutive duplicate points were generated
  for (size_t i = 1; i < with_gaps.size(); ++i)
  {
    EXPECT_NE(with_gaps[i], with_gaps[i - 1]);
  }
}

TEST(CellPath, AddNoiseDeformation)
{
  hmap::CellPath line(std::vector<glm::ivec2>{{0, 10}, {50, 10}});
  hmap::enforce_path_adjacency(line);

  glm::ivec2 shape = {64, 64};
  hmap::add_noise(line, 42, 0.1f, 5.f, shape);

  EXPECT_TRUE(hmap::is_path_adjacent(line));
  for (const auto &p : line)
  {
    EXPECT_GE(p.x, 0);
    EXPECT_LT(p.x, shape.x);
    EXPECT_GE(p.y, 0);
    EXPECT_LT(p.y, shape.y);
  }
}
