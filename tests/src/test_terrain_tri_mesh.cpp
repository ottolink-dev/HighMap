#include <cmath>
#include <vector>

#include "highmap/terrain_tri_mesh.hpp"
#include "highmap/virtual_array/virtual_array.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(TerrainTriMeshTest, EmptyMeshSanity)
{
  TerrainTriMesh mesh;
  EXPECT_EQ(mesh.size(), 0u);
  EXPECT_TRUE(mesh.get_points().empty());
  EXPECT_TRUE(mesh.get_triangles().empty());
  EXPECT_TRUE(mesh.get_convex_hull().empty());

  int last_tri = 0;
  EXPECT_EQ(mesh.find_triangle({0.5f, 0.5f}), -1);
  EXPECT_FLOAT_EQ(mesh.interpolate_z_linear({0.5f, 0.5f}, last_tri, -1.f),
                  -1.f);
  EXPECT_FLOAT_EQ(mesh.interpolate_z_nearest({0.5f, 0.5f}), 0.f);

  auto bbox = mesh.get_bbox();
  EXPECT_FALSE(bbox.contains({0.5f, 0.5f}));

  mesh.remap_z(0.f, 1.f);
  mesh.relax_xy(0.5f, 1);
  mesh.relax_xyz(0.5f, 1);
}

TEST(TerrainTriMeshTest, BasicTriangulationAndBarycentric)
{
  std::vector<glm::vec3> pts = {{0.f, 0.f, 0.f},
                                {1.f, 0.f, 0.5f},
                                {1.f, 1.f, 1.0f},
                                {0.f, 1.f, 0.5f}};

  TerrainTriMesh mesh(pts);
  EXPECT_EQ(mesh.size(), 4u);
  EXPECT_EQ(mesh.get_triangles().size(), 2u);
  EXPECT_EQ(mesh.get_convex_hull().size(), 4u);

  int   last_tri = 0;
  float z_center = mesh.interpolate_z_linear({0.5f, 0.5f}, last_tri);
  EXPECT_NEAR(z_center, 0.5f, 1e-4f);

  float z_near = mesh.interpolate_z_nearest({0.1f, 0.1f});
  EXPECT_FLOAT_EQ(z_near, 0.f);

  float z_approx = mesh.interpolate_z_nearest_approx({0.1f, 0.1f}, last_tri);
  EXPECT_FLOAT_EQ(z_approx, 0.f);
}

TEST(TerrainTriMeshTest, RelaxAndSmoothing)
{
  std::vector<glm::vec3> pts = {{0.f, 0.f, 0.f},
                                {1.f, 0.f, 0.f},
                                {1.f, 1.f, 0.f},
                                {0.f, 1.f, 0.f},
                                {0.5f, 0.5f, 1.0f}};

  TerrainTriMesh mesh(pts);
  EXPECT_EQ(mesh.size(), 5u);

  mesh.relax_xyz(0.5f, 2, true);
  EXPECT_LT(mesh.get_points()[4].z, 1.0f);
  EXPECT_GT(mesh.get_points()[4].z, 0.0f);

  EXPECT_FLOAT_EQ(mesh.get_points()[0].z, 0.f);
  EXPECT_FLOAT_EQ(mesh.get_points()[1].z, 0.f);
  EXPECT_FLOAT_EQ(mesh.get_points()[2].z, 0.f);
  EXPECT_FLOAT_EQ(mesh.get_points()[3].z, 0.f);
}

TEST(TerrainTriMeshTest, SlopeLimiter)
{
  std::vector<glm::vec3> pts = {{0.f, 0.f, 0.f},
                                {1.f, 0.f, 0.f},
                                {1.f, 1.f, 0.f},
                                {0.f, 1.f, 0.f},
                                {0.5f, 0.5f, 10.0f}};

  TerrainTriMesh mesh(pts);
  mesh.slope_limiter(1.0f, 20, 0.5f);

  EXPECT_LT(mesh.get_points()[4].z, 10.0f);
}

TEST(TerrainTriMeshTest, GenerateFromHeightmapRandomTiledVsNonTiled)
{
  const glm::ivec2 shape(128, 128);
  const int        control_points_count = 150;
  const uint32_t   seed = 42;

  Array z(shape);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
      z(i, j) = std::sin(0.05f * i) * std::cos(0.07f * j) + 0.01f * i;

  TerrainTriMesh mesh_ref = generate_terrain_tri_mesh_from_heightmap_random(
      z,
      control_points_count,
      seed);

  VirtualArray va_z(shape, glm::ivec2(64, 64), 4);
  va_z.from_array(z, ComputeMode{});

  for (ForEachMode mode :
       {ForEachMode::VA_SEQUENTIAL, ForEachMode::VA_DISTRIBUTED})
  {
    ComputeMode cm;
    cm.mode = mode;

    TerrainTriMesh mesh_va;

    hmap::for_each_tile(
        {&va_z},
        {},
        [&](std::vector<const hmap::Array *> p_arrays_in,
            std::vector<hmap::Array *>,
            const hmap::UvTile &region,
            ICommunicator      &comm)
        {
          auto [pa_z] = unpack<1>(p_arrays_in);

          TerrainTriMesh tile_mesh =
              hmap::generate_terrain_tri_mesh_from_heightmap_random(
                  *pa_z,
                  control_points_count,
                  seed,
                  region,
                  &comm);

          if (comm.rank() == 0) mesh_va = std::move(tile_mesh);
        },
        cm);

    ASSERT_EQ(mesh_va.size(), mesh_ref.size());
    ASSERT_EQ(mesh_va.get_triangles().size(), mesh_ref.get_triangles().size());

    const auto &pts_ref = mesh_ref.get_points();
    const auto &pts_va = mesh_va.get_points();
    for (size_t i = 0; i < pts_ref.size(); ++i)
    {
      EXPECT_NEAR(pts_va[i].x, pts_ref[i].x, 1e-5f);
      EXPECT_NEAR(pts_va[i].y, pts_ref[i].y, 1e-5f);
      EXPECT_NEAR(pts_va[i].z, pts_ref[i].z, 5e-2f);
    }
  }
}
