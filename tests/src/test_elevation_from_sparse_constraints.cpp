#include "highmap/authoring.hpp"
#include "highmap/internal/assert.hpp"
#include "highmap/opencl/gpu_opencl.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(ElevationFromSparseConstraints, ExactConstraintsPreservation)
{
  gpu::init_opencl();

  glm::ivec2 shape = {64, 64};

  // Mountain ridge along center row 32
  Array mountains(shape, 0.0f);
  for (int i = 0; i < shape.x; ++i)
    mountains(i, 32) = 1.0f;

  // Coastlines at rows 16 and 48
  Array coastline(shape, 0.0f);
  for (int i = 0; i < shape.x; ++i)
  {
    coastline(i, 16) = 1.0f;
    coastline(i, 48) = 1.0f;
  }

  // Boundary edges at row 0 and 63
  Array boundary(shape, 0.0f);
  for (int i = 0; i < shape.x; ++i)
  {
    boundary(i, 0) = 1.0f;
    boundary(i, shape.y - 1) = 1.0f;
  }

  Array z = elevation_from_sparse_constraints(mountains,
                                              &coastline,
                                              &boundary,
                                              1.0f,
                                              0.0f,
                                              -1.0f,
                                              500,
                                              1e-5f);

  // Exact Dirichlet constraint preservation
  for (int i = 0; i < shape.x; ++i)
  {
    EXPECT_NEAR(z(i, 32), 1.0f, 1e-4f);
    EXPECT_NEAR(z(i, 16), 0.0f, 1e-4f);
    EXPECT_NEAR(z(i, 48), 0.0f, 1e-4f);
    EXPECT_NEAR(z(i, 0), -1.0f, 1e-4f);
    EXPECT_NEAR(z(i, shape.y - 1), -1.0f, 1e-4f);
  }

  // Intermediate values are strictly monotonically bounded
  // (Between mountain +1.0 and coastline 0.0)
  for (int j = 17; j < 32; ++j)
  {
    EXPECT_GT(z(32, j), 0.0f);
    EXPECT_LT(z(32, j), 1.0f);
  }

  // (Between coastline 0.0 and boundary -1.0)
  for (int j = 1; j < 16; ++j)
  {
    EXPECT_LT(z(32, j), 0.0f);
    EXPECT_GT(z(32, j), -1.0f);
  }
}

TEST(ElevationFromSparseConstraints, DefaultPerimeterBoundary)
{
  gpu::init_opencl();

  glm::ivec2 shape = {65, 65};

  // Mountain peak at center
  Array mountains(shape, 0.0f);
  mountains(32, 32) = 1.0f;

  Array z = elevation_from_sparse_constraints(mountains,
                                              nullptr,
                                              nullptr,
                                              1.0f,
                                              0.0f,
                                              -1.0f,
                                              500,
                                              1e-5f);

  EXPECT_NEAR(z(32, 32), 1.0f, 1e-3f);
  EXPECT_NEAR(z(0, 0), -1.0f, 1e-4f);
  EXPECT_NEAR(z(64, 64), -1.0f, 1e-4f);

  // Elevation decreases monotonically away from center peak toward boundary
  for (int i = 32; i < 64; ++i)
  {
    EXPECT_GE(z(i, 32), z(i + 1, 32) - 1e-5f);
  }
}

TEST(HarmonicInterpolation, VaryingDiffusionAnisotropy)
{
  gpu::init_opencl();

  glm::ivec2 shape = {65, 65};

  Array u_init(shape, 0.0f);
  Array mask(shape, 0.0f);
  u_init(32, 32) = 1.0f;
  mask(32, 32) = 1.0f;

  for (int i = 0; i < shape.x; ++i)
  {
    u_init(i, 0) = -1.0f;
    mask(i, 0) = 1.0f;
    u_init(i, shape.y - 1) = -1.0f;
    mask(i, shape.y - 1) = 1.0f;
  }
  for (int j = 0; j < shape.y; ++j)
  {
    u_init(0, j) = -1.0f;
    mask(0, j) = 1.0f;
    u_init(shape.x - 1, j) = -1.0f;
    mask(shape.x - 1, j) = 1.0f;
  }

  // Dx >> Dy (strong horizontal diffusion)
  Array dx_high(shape, 10.0f);
  Array dy_low(shape, 1.0f);

  Array z_aniso = gpu::harmonic_interpolation(u_init,
                                              mask,
                                              dx_high,
                                              dy_low,
                                              1000,
                                              1e-5f);

  // Peak preserved
  EXPECT_NEAR(z_aniso(32, 32), 1.0f, 1e-3f);
  EXPECT_NEAR(z_aniso(0, 0), -1.0f, 1e-4f);

  // Horizontal decay is slower than vertical decay due to Dx >> Dy
  // z(32 + 10, 32) > z(32, 32 + 10)
  EXPECT_GT(z_aniso(42, 32), z_aniso(32, 42));
}

TEST(HarmonicInterpolation, CpuGpuConsistencyWithDxDy)
{
  gpu::init_opencl();

  glm::ivec2 shape = {33, 33};
  Array      u_init(shape, 0.0f);
  Array      mask(shape, 0.0f);

  u_init(16, 16) = 1.0f;
  mask(16, 16) = 1.0f;

  for (int i = 0; i < shape.x; ++i)
  {
    u_init(i, 0) = -1.0f;
    mask(i, 0) = 1.0f;
    u_init(i, shape.y - 1) = -1.0f;
    mask(i, shape.y - 1) = 1.0f;
  }
  for (int j = 0; j < shape.y; ++j)
  {
    u_init(0, j) = -1.0f;
    mask(0, j) = 1.0f;
    u_init(shape.x - 1, j) = -1.0f;
    mask(shape.x - 1, j) = 1.0f;
  }

  // Spatially varying Dx, Dy fields
  Array dx(shape, 1.0f);
  Array dy(shape, 1.0f);
  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      dx(i, j) = 1.0f + 2.0f * static_cast<float>(i) / 33.0f;
      dy(i, j) = 0.5f + 1.5f * static_cast<float>(j) / 33.0f;
    }

  Array z_cpu = harmonic_interpolation(u_init, mask, dx, dy, 1000, 1e-6f);
  Array z_gpu = gpu::harmonic_interpolation(u_init, mask, dx, dy, 1000, 1e-6f);

  for (int j = 0; j < shape.y; ++j)
    for (int i = 0; i < shape.x; ++i)
    {
      EXPECT_NEAR(z_cpu(i, j), z_gpu(i, j), 1e-3f);
    }
}
