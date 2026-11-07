#include "highmap.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/virtual_array/tile_storage.hpp"

#include <gtest/gtest.h>

using namespace hmap;

namespace
{

const std::vector<ForEachMode> k_compute_modes = {
    ForEachMode::VA_SEQUENTIAL,
    ForEachMode::VA_DISTRIBUTED,
    ForEachMode::VA_SINGLE_ARRAY};

const std::vector<StorageMode> k_storage_modes = {
    StorageMode::VA_RAM,
    StorageMode::VA_DISK_LRU,
    StorageMode::VA_DISK_LRU_MIN};

const std::vector<bool> k_trim_options = {false, true};

} // namespace

class ComputeModeMatrixTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    gpu::init_opencl();
  }
};

// --- Remap matrix test

TEST_F(ComputeModeMatrixTest, RemapMatrix)
{
  const glm::ivec2 shape{128, 128};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 0;

  Array input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 1234);

  // Reference monolithic result
  Array ref = input;
  remap(ref, 0.f, 1.f);

  for (auto mode : k_compute_modes)
  {
    for (auto storage_mode : k_storage_modes)
    {
      for (bool trim : k_trim_options)
      {
        ComputeMode cm{.mode = mode, .trim_storage = trim};

        VirtualArray va(shape, tile_shape, halo, storage_mode);
        va.from_array(input, cm);

        for_each_tile(
            va,
            [](Array &tile, const TileRegion &, ICommunicator &comm)
            { remap(tile, 0.f, 1.f, comm); },
            cm);

        Array res = va.to_array(cm);
        Array diff = abs(res - ref);

        EXPECT_LT(diff.max(), 1e-4f)
            << "Remap failed for mode: " << to_string(mode)
            << ", storage: " << to_string(storage_mode)
            << ", trim: " << trim;
      }
    }
  }
}

// --- SmoothCPulse matrix test

TEST_F(ComputeModeMatrixTest, SmoothCPulseMatrix)
{
  const glm::ivec2 shape{128, 128};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 8;
  const int        ir = 4;

  Array input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 5678);
  remap(input, 0.f, 1.f);

  // Reference single-array result with comm
  Array ref = input;
  smooth_cpulse(ref, ir, nullptr, TileRegion::dummy(), SingleArrayComm::instance());

  for (auto mode : k_compute_modes)
  {
    for (auto storage_mode : k_storage_modes)
    {
      for (bool trim : k_trim_options)
      {
        ComputeMode cm{.mode = mode, .trim_storage = trim};

        VirtualArray va(shape, tile_shape, halo, storage_mode);
        va.from_array(input, cm);

        for_each_tile(
            va,
            [ir](Array &tile, const TileRegion &region, ICommunicator &comm)
            { smooth_cpulse(tile, ir, nullptr, region, comm); },
            cm);

        Array res = va.to_array(cm);
        Array diff = abs(res - ref);

        EXPECT_LT(diff.max(), 5e-4f)
            << "SmoothCPulse failed for mode: " << to_string(mode)
            << ", storage: " << to_string(storage_mode)
            << ", trim: " << trim;
      }
    }
  }
}

// --- DistanceTransform matrix test

TEST_F(ComputeModeMatrixTest, DistanceTransformMatrix)
{
  const glm::ivec2 shape{128, 128};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 8;

  Array mask = white_sparse(shape, 0.f, 1.f, 0.03f, 9999);

  // Reference monolithic result
  Array ref = distance_transform(mask, false, TileRegion::dummy(), SingleArrayComm::instance());

  for (auto mode : k_compute_modes)
  {
    for (auto storage_mode : k_storage_modes)
    {
      for (bool trim : k_trim_options)
      {
        ComputeMode cm{.mode = mode, .trim_storage = trim};

        VirtualArray va(shape, tile_shape, halo, storage_mode);
        va.from_array(mask, cm);

        for_each_tile(
            va,
            [](Array &tile, const TileRegion &region, ICommunicator &comm)
            { tile = distance_transform(tile, false, region, comm); },
            cm);

        Array res = va.to_array(cm);
        Array diff = abs(res - ref);

        EXPECT_LT(diff.max(), 1e-4f)
            << "DistanceTransform failed for mode: " << to_string(mode)
            << ", storage: " << to_string(storage_mode)
            << ", trim: " << trim;
      }
    }
  }
}

// --- SnowSimulation matrix test

TEST_F(ComputeModeMatrixTest, SnowSimulationMatrix)
{
  const glm::ivec2 shape{128, 128};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 8;

  Array z_arr = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 42);
  remap(z_arr, 0.f, 1.f);

  // Reference tiled result (VA_DISTRIBUTED + VA_RAM)
  Array ref_snow_tiled;
  Array ref_snow_single;

  for (auto mode : k_compute_modes)
  {
    for (auto storage_mode : k_storage_modes)
    {
      for (bool trim : k_trim_options)
      {
        ComputeMode cm{.mode = mode, .trim_storage = trim};

        VirtualArray z(shape, tile_shape, halo, storage_mode);
        VirtualArray fall_map(shape, tile_shape, halo, storage_mode);
        VirtualArray melting_map(shape, tile_shape, halo, storage_mode);
        VirtualArray talus_map(shape, tile_shape, halo, storage_mode);
        VirtualArray snow_depth_out(shape, tile_shape, halo, storage_mode);

        z.from_array(z_arr, cm);
        fall_map.fill(1.f, cm);
        melting_map.fill(0.f, cm);
        talus_map.fill(1.5f / float(shape.x), cm);

        for_each_tile(
            {&z, &fall_map, &melting_map, &talus_map},
            {&snow_depth_out},
            [](std::vector<const Array *> p_arrays_in,
               std::vector<Array *>       p_arrays_out,
               const TileRegion          &region,
               ICommunicator             &comm)
            {
              auto [pa_z, pa_fall, pa_melt, pa_talus] = unpack<4>(p_arrays_in);
              auto [pa_snow] = unpack<1>(p_arrays_out);

              *pa_snow = gpu::snow_simulation(*pa_z,
                                              0.1f,
                                              *pa_fall,
                                              *pa_melt,
                                              *pa_talus,
                                              30,
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
            },
            cm);

        Array res = snow_depth_out.to_array(cm);
        EXPECT_GT(res.sum(), 0.f);

        if (mode == ForEachMode::VA_SINGLE_ARRAY)
        {
          if (ref_snow_single.shape.x == 0)
          {
            ref_snow_single = res;
          }
          Array diff = abs(res - ref_snow_single);
          EXPECT_LT(diff.max(), 1e-4f)
              << "SnowSimulation single-array mismatch for storage: "
              << to_string(storage_mode) << ", trim: " << trim;
        }
        else
        {
          if (ref_snow_tiled.shape.x == 0)
          {
            ref_snow_tiled = res;
          }
          Array diff = abs(res - ref_snow_tiled);
          EXPECT_LT(diff.max(), 1e-4f)
              << "SnowSimulation tiled mismatch for mode: " << to_string(mode)
              << ", storage: " << to_string(storage_mode)
              << ", trim: " << trim;
        }
      }
    }
  }

  // Cross-verification: tiled simulation with halo blending vs un-tiled single array
  if (ref_snow_tiled.shape.x > 0 && ref_snow_single.shape.x > 0)
  {
    Array diff_tiled_single = abs(ref_snow_tiled - ref_snow_single);
    EXPECT_LT(diff_tiled_single.max(), 0.03f);
  }
}
