#include "highmap.hpp"
#include "highmap/internal/assert.hpp"

#include <gtest/gtest.h>

using namespace hmap;

namespace
{

Array helper_generate_array()
{
  const glm::ivec2 shape = {512, 512};
  const glm::vec2  kw = {4.f, 4.f};
  const int        seed = 0;

  Array z = noise_fbm(NoiseType::PERLIN, shape, kw, seed);
  return z;
}

} // namespace

TEST(Advection, MultiField)
{
  Array z = helper_generate_array();

  // Create two distinct arrays of the same shape
  Array f1(z.shape);
  Array f2(z.shape);
  for (int j = 0; j < z.shape.x; ++j)
  {
    for (int i = 0; i < z.shape.y; ++i)
    {
      f1(i, j) = (i % 2 == 0) ? 1.f : 0.f;
      f2(i, j) = (j % 2 == 0) ? 0.8f : 0.2f;
    }
  }

  int           nparticles = 1;
  std::uint32_t seed = 42;
  bool          reverse = false;
  bool          post_filter = true;
  float         post_filter_sigma = 0.125f;
  float         advection_length = 0.1f;
  float         value_persistence = 0.99f;
  float         inertia = 0.2f;

  // 1. Advect individually
  Array out1_single_a = gpu::advection_particle(z,
                                                f1,
                                                nparticles,
                                                seed,
                                                reverse,
                                                post_filter,
                                                post_filter_sigma,
                                                advection_length,
                                                value_persistence,
                                                inertia);
  Array out1_single_b = gpu::advection_particle(z,
                                                f1,
                                                nparticles,
                                                seed,
                                                reverse,
                                                post_filter,
                                                post_filter_sigma,
                                                advection_length,
                                                value_persistence,
                                                inertia);

  float max_diff_det = 0.f;
  for (int idx = 0; idx < z.shape.x * z.shape.y; ++idx)
  {
    max_diff_det = std::max(
        max_diff_det,
        std::abs(out1_single_a.vector[idx] - out1_single_b.vector[idx]));
  }
  std::cout << "Max deterministic diff (1 particle): " << max_diff_det
            << std::endl;

  bool ret_det = assert_almost_equal(out1_single_a, out1_single_b, 1e-5f);
  EXPECT_EQ(ret_det, true);

  // 2. Test multi-field vs single-field with 1 particle
  std::vector<Array> fields = {f1, f2};
  std::vector<Array> out_multi = gpu::advection_particle(z,
                                                         fields,
                                                         nparticles,
                                                         seed,
                                                         reverse,
                                                         post_filter,
                                                         post_filter_sigma,
                                                         advection_length,
                                                         value_persistence,
                                                         inertia);
  Array              out2_single_a = gpu::advection_particle(z,
                                                f2,
                                                nparticles,
                                                seed,
                                                reverse,
                                                post_filter,
                                                post_filter_sigma,
                                                advection_length,
                                                value_persistence,
                                                inertia);

  bool ret_multi1 = assert_almost_equal(out1_single_a, out_multi[0], 1e-5f);
  EXPECT_EQ(ret_multi1, true);

  bool ret_multi2 = assert_almost_equal(out2_single_a, out_multi[1], 1e-5f);
  EXPECT_EQ(ret_multi2, true);

  // 3. Test with many particles to ensure stability (no crashes/hangs)
  int                nparticles_many = 50000;
  std::vector<Array> out_multi_many = gpu::advection_particle(z,
                                                              fields,
                                                              nparticles_many,
                                                              seed,
                                                              reverse,
                                                              post_filter,
                                                              post_filter_sigma,
                                                              advection_length,
                                                              value_persistence,
                                                              inertia);
  EXPECT_EQ(out_multi_many.size(), 2);
}

TEST(Advection, Multiscale)
{
  Array z = helper_generate_array();

  Array f1(z.shape);
  Array f2(z.shape);
  for (int j = 0; j < z.shape.x; ++j)
  {
    for (int i = 0; i < z.shape.y; ++i)
    {
      f1(i, j) = (i % 4 == 0) ? 1.f : 0.f;
      f2(i, j) = (j % 4 == 0) ? 0.9f : 0.1f;
    }
  }

  std::uint32_t    seed = 42;
  std::vector<int> steps_per_level = {2, 1};
  float            particles_ratio = 0.1f;

  // 1. Single field multiscale
  Array out_single = gpu::advection_particle_multiscale(z,
                                                        f1,
                                                        seed,
                                                        steps_per_level,
                                                        particles_ratio);

  EXPECT_EQ(out_single.shape, z.shape);
  EXPECT_FALSE(out_single.vector.empty());

  // 2. Multi-field multiscale
  std::vector<Array> fields = {f1, f2};
  std::vector<Array> out_multi = gpu::advection_particle_multiscale(
      z,
      fields,
      seed,
      steps_per_level,
      particles_ratio);

  EXPECT_EQ(out_multi.size(), 2);
  EXPECT_EQ(out_multi[0].shape, z.shape);
  EXPECT_EQ(out_multi[1].shape, z.shape);

  // Check sanity: values should be non-empty and within bounds [0, 1]
  for (int i = 0; i < z.shape.x * z.shape.y; ++i)
  {
    EXPECT_GE(out_single.vector[i], 0.f);
    EXPECT_LE(out_single.vector[i], 1.f);
    EXPECT_GE(out_multi[0].vector[i], 0.f);
    EXPECT_LE(out_multi[0].vector[i], 1.f);
    EXPECT_GE(out_multi[1].vector[i], 0.f);
    EXPECT_LE(out_multi[1].vector[i], 1.f);
  }
}
