#include "highmap.hpp"
#include "highmap/dbg/assert.hpp"

#include <gtest/gtest.h>

using namespace hmap;

TEST(SnowSimulation, ZeroIterationsReturnsEmpty)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z);

  Array fall_map(shape, 1.f);
  Array melting_map(shape, 0.f);
  Array talus(shape, 2.5f / shape.x);

  Array s =
      gpu::snow_simulation(z, 0.1f, fall_map, melting_map, talus, 0, 0.5f);

  EXPECT_EQ(s.shape, shape);
  EXPECT_FLOAT_EQ(s.max(), 0.f);
  EXPECT_FLOAT_EQ(s.min(), 0.f);
}

TEST(SnowSimulation, IterationParityAndDeterminism)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z);

  Array fall_map(shape, 1.f);
  Array melting_map(shape, 0.f);
  Array talus(shape, 2.5f / shape.x);

  Array s_even =
      gpu::snow_simulation(z, 0.1f, fall_map, melting_map, talus, 50, 0.5f);
  Array s_odd =
      gpu::snow_simulation(z, 0.1f, fall_map, melting_map, talus, 51, 0.5f);
  Array s_even_again =
      gpu::snow_simulation(z, 0.1f, fall_map, melting_map, talus, 50, 0.5f);

  // deterministic within float precision
  EXPECT_NEAR(s_even.sum(), s_even_again.sum(), 1e-2f);

  // evolution check
  const float sum_even = s_even.sum();
  const float sum_odd = s_odd.sum();
  EXPECT_GT(sum_even, 0.f);
  EXPECT_NEAR(sum_odd, sum_even, 0.05f * sum_even);
  EXPECT_GE(s_odd.min(), 0.f);
}

TEST(SnowSimulation, SnowAccumulationWithoutMelting)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 7);
  remap(z);

  Array fall_map(shape, 1.f);
  Array melting_map(shape, 0.f);
  Array talus(shape, 2.5f / shape.x);

  Array s =
      gpu::snow_simulation(z, 0.1f, fall_map, melting_map, talus, 100, 0.5f);

  EXPECT_GT(s.sum(), 0.f);
  EXPECT_GE(s.min(), 0.f);
}

TEST(SnowSimulation, MeltingReducesSnowVolume)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {64, 64};
  Array            z = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z);

  Array fall_map(shape, 1.f);
  Array no_melt(shape, 0.f);
  Array melt(shape, 1.f);
  Array talus(shape, 2.5f / shape.x);

  Array s_no_melt =
      gpu::snow_simulation(z, 0.1f, fall_map, no_melt, talus, 100, 0.5f);
  Array s_melt =
      gpu::snow_simulation(z, 0.1f, fall_map, melt, talus, 100, 0.5f);

  EXPECT_GT(s_no_melt.sum(), s_melt.sum());
  EXPECT_GE(s_melt.min(), 0.f);
}

TEST(SnowSimulation, TiledMinStorageLruMin)
{
  gpu::init_opencl();

  const glm::ivec2 shape = {128, 128};
  const glm::ivec2 tile_shape = {32, 32}; // 4x4 = 16 tiles
  const int        halo = 8;

  Array z_arr = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z_arr);

  for (auto mode : {ForEachMode::VA_SEQUENTIAL, ForEachMode::VA_DISTRIBUTED})
  {
    ComputeMode cm{.mode = mode, .trim_storage = true};

    VirtualArray z(shape, tile_shape, halo, StorageMode::VA_DISK_LRU_MIN);
    VirtualArray fall_map(shape,
                          tile_shape,
                          halo,
                          StorageMode::VA_DISK_LRU_MIN);
    VirtualArray melting_map(shape,
                             tile_shape,
                             halo,
                             StorageMode::VA_DISK_LRU_MIN);
    VirtualArray talus_map(shape,
                           tile_shape,
                           halo,
                           StorageMode::VA_DISK_LRU_MIN);
    VirtualArray z_out(shape, tile_shape, halo, StorageMode::VA_DISK_LRU_MIN);
    VirtualArray snow_depth_out(shape,
                                tile_shape,
                                halo,
                                StorageMode::VA_DISK_LRU_MIN);

    z.from_array(z_arr, cm);
    fall_map.fill(1.f, cm);
    melting_map.fill(0.f, cm);
    talus_map.fill(1.5f / float(shape.x), cm);

    for_each_tile(
        {&z, &fall_map, &melting_map, &talus_map},
        {&z_out, &snow_depth_out},
        [&](const std::vector<const Array *> &p_arrays_in,
            std::vector<Array *>             &p_arrays_out,
            const TileRegion                 &region,
            ICommunicator                    &comm)
        {
          auto [pa_z, pa_fall, pa_melt, pa_talus] = unpack<4>(p_arrays_in);
          auto [pa_z_out, pa_snow] = unpack<2>(p_arrays_out);

          *pa_snow = gpu::snow_simulation(*pa_z,
                                          0.1f,
                                          *pa_fall,
                                          *pa_melt,
                                          *pa_talus,
                                          40,
                                          0.5f,
                                          1.f,
                                          0.5f,
                                          0.2f,
                                          0.8f,
                                          0.8f,
                                          1.5f,
                                          0.1f,
                                          /* post_filter */ true,
                                          0.8f,
                                          false,
                                          0.f,
                                          true,
                                          10,
                                          region,
                                          comm);

          *pa_z_out = *pa_z + *pa_snow;
        },
        cm);

    Array s_res = snow_depth_out.to_array(cm);
    EXPECT_GT(s_res.sum(), 0.f);
    EXPECT_GE(s_res.min(), -1e-5f);
    EXPECT_LE(snow_depth_out.live_tile_count(), 2u);
  }
}
