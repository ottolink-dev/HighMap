#include "highmap.hpp"
#include "highmap/dbg/assert.hpp"
#include "highmap/virtual_array/tile_storage.hpp"

#include <gtest/gtest.h>

using namespace hmap;

namespace
{

const std::vector<ForEachMode> k_compute_modes = {ForEachMode::VA_SEQUENTIAL,
                                                  ForEachMode::VA_DISTRIBUTED,
                                                  ForEachMode::VA_SINGLE_ARRAY};

const std::vector<StorageMode> k_storage_modes = {StorageMode::VA_RAM,
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

// --- Array Methods with Communicator matrix test (min, max, mean, sum, std,
// ptp, range, unique_values, normalize)

TEST_F(ComputeModeMatrixTest, ArrayMethodsWithCommMatrix)
{
  const glm::ivec2 shape{128, 128};
  const glm::ivec2 tile_shape{32, 32};
  const int        halo = 0;

  Array input = noise_fbm(NoiseType::PERLIN, shape, {2.f, 2.f}, 777);
  // Introduce a few discrete values to test unique_values
  for (int j = 0; j < shape.y; ++j)
  {
    for (int i = 0; i < shape.x; ++i)
    {
      input(i, j) = std::round(input(i, j) * 10.f);
    }
  }

  // Monolithic reference metrics
  const float              ref_min = input.min();
  const float              ref_max = input.max();
  const float              ref_mean = input.mean();
  const float              ref_sum = input.sum();
  const float              ref_std = input.std();
  const float              ref_ptp = input.ptp();
  const glm::vec2          ref_range = input.range();
  const std::vector<float> ref_unique = input.unique_values();

  Array ref_norm = input;
  ref_norm.normalize();

  for (auto mode : k_compute_modes)
  {
    for (auto storage_mode : k_storage_modes)
    {
      for (bool trim : k_trim_options)
      {
        ComputeMode cm{.mode = mode, .trim_storage = trim};

        VirtualArray va(shape, tile_shape, halo, storage_mode);
        va.from_array(input, cm);

        // 1. Reductions & metrics using ICommunicator
        std::atomic<int>   tile_count{0};
        float              calc_min = 0.f;
        float              calc_max = 0.f;
        float              calc_mean = 0.f;
        float              calc_sum = 0.f;
        float              calc_std = 0.f;
        float              calc_ptp = 0.f;
        glm::vec2          calc_range{0.f, 0.f};
        std::vector<float> calc_unique;

        for_each_tile(
            va,
            [&](const Array &tile, const TileRegion &, ICommunicator &comm)
            {
              float              t_min = tile.min(comm);
              float              t_max = tile.max(comm);
              float              t_mean = tile.mean(comm);
              float              t_sum = tile.sum(comm);
              float              t_std = tile.std(comm);
              float              t_ptp = tile.ptp(comm);
              glm::vec2          t_range = tile.range(comm);
              std::vector<float> t_unique = tile.unique_values(comm);

              if (comm.rank() == 0)
              {
                calc_min = t_min;
                calc_max = t_max;
                calc_mean = t_mean;
                calc_sum = t_sum;
                calc_std = t_std;
                calc_ptp = t_ptp;
                calc_range = t_range;
                calc_unique = t_unique;
              }
              ++tile_count;
            },
            cm);

        EXPECT_NEAR(calc_min, ref_min, 1e-4f)
            << "min failed for mode=" << to_string(mode)
            << " storage=" << to_string(storage_mode) << " trim=" << trim;
        EXPECT_NEAR(calc_max, ref_max, 1e-4f)
            << "max failed for mode=" << to_string(mode)
            << " storage=" << to_string(storage_mode) << " trim=" << trim;
        EXPECT_NEAR(calc_mean, ref_mean, 1e-4f)
            << "mean failed for mode=" << to_string(mode)
            << " storage=" << to_string(storage_mode) << " trim=" << trim;
        EXPECT_NEAR(calc_sum, ref_sum, 1e-2f)
            << "sum failed for mode=" << to_string(mode)
            << " storage=" << to_string(storage_mode) << " trim=" << trim;
        EXPECT_NEAR(calc_std, ref_std, 1e-4f)
            << "std failed for mode=" << to_string(mode)
            << " storage=" << to_string(storage_mode) << " trim=" << trim;
        EXPECT_NEAR(calc_ptp, ref_ptp, 1e-4f)
            << "ptp failed for mode=" << to_string(mode)
            << " storage=" << to_string(storage_mode) << " trim=" << trim;
        EXPECT_NEAR(calc_range.x, ref_range.x, 1e-4f)
            << "range.x failed for mode=" << to_string(mode)
            << " storage=" << to_string(storage_mode) << " trim=" << trim;
        EXPECT_NEAR(calc_range.y, ref_range.y, 1e-4f)
            << "range.y failed for mode=" << to_string(mode)
            << " storage=" << to_string(storage_mode) << " trim=" << trim;
        EXPECT_EQ(calc_unique, ref_unique)
            << "unique_values failed for mode=" << to_string(mode)
            << " storage=" << to_string(storage_mode) << " trim=" << trim;

        // 2. In-place normalize using ICommunicator
        VirtualArray va_norm(shape, tile_shape, halo, storage_mode);
        va_norm.from_array(input, cm);

        for_each_tile(
            va_norm,
            [](Array &tile, const TileRegion &, ICommunicator &comm)
            { tile.normalize(comm); },
            cm);

        Array res_norm = va_norm.to_array(cm);
        Array diff_norm = abs(res_norm - ref_norm);

        EXPECT_LT(diff_norm.max(), 1e-4f)
            << "normalize failed for mode=" << to_string(mode)
            << " storage=" << to_string(storage_mode) << " trim=" << trim;
      }
    }
  }
}

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
            << ", storage: " << to_string(storage_mode) << ", trim: " << trim;
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
  smooth_cpulse(ref,
                ir,
                nullptr,
                TileRegion::dummy(),
                SingleArrayComm::instance());

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
            << ", storage: " << to_string(storage_mode) << ", trim: " << trim;
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
  Array ref = distance_transform(mask,
                                 false,
                                 TileRegion::dummy(),
                                 SingleArrayComm::instance());

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
            << ", storage: " << to_string(storage_mode) << ", trim: " << trim;
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
              << ", storage: " << to_string(storage_mode) << ", trim: " << trim;
        }
      }
    }
  }

  // Cross-verification: tiled simulation with halo blending vs un-tiled single
  // array
  if (ref_snow_tiled.shape.x > 0 && ref_snow_single.shape.x > 0)
  {
    Array diff_tiled_single = abs(ref_snow_tiled - ref_snow_single);
    EXPECT_LT(diff_tiled_single.max(), 0.03f);
  }
}
TEST_F(ComputeModeMatrixTest, SnowSimulation4KMinimalMemory)
{
  std::cout << "Starting SnowRepro (shape 4096, tile 1024, VA_DISK_LRU_MIN)..." << std::endl;
  glm::ivec2 shape = {4096, 4096};
  glm::ivec2 tile_shape = {1024, 1024}; // 4x4 = 16 tiles
  int        halo = 32;

  VirtualArray va_z(shape, {-1.f, 1.f, -1.f, 1.f}, tile_shape, halo, StorageMode::VA_DISK_LRU_MIN);
  VirtualArray va_snow(shape, {-1.f, 1.f, -1.f, 1.f}, tile_shape, halo, StorageMode::VA_DISK_LRU_MIN);
  VirtualArray va_depth(shape, {-1.f, 1.f, -1.f, 1.f}, tile_shape, halo, StorageMode::VA_DISK_LRU_MIN);
  VirtualArray va_talus(shape, {-1.f, 1.f, -1.f, 1.f}, tile_shape, halo, StorageMode::VA_DISK_LRU_MIN);
  VirtualArray va_melting(shape, {-1.f, 1.f, -1.f, 1.f}, tile_shape, halo, StorageMode::VA_DISK_LRU_MIN);

  ComputeMode cm;
  cm.mode = ForEachMode::VA_SEQUENTIAL;
  cm.trim_storage = true;

  std::cout << "Init va_z with noise_fbm..." << std::endl;
  for_each_tile(
      {},
      {&va_z},
      [](auto, auto out, const TileRegion &region, auto &)
      {
        *out[0] = gpu::noise_fbm(NoiseType::SIMPLEX2,
                                 out[0]->shape,
                                 glm::vec2(4.f, 4.f),
                                 1337,
                                 4,
                                 1.f,
                                 0.5f,
                                 2.0f,
                                 nullptr,
                                 nullptr,
                                 nullptr,
                                 region.bbox);
      },
      cm);

  std::cout << "Init other VAs..." << std::endl;
  for_each_tile(
      {},
      {&va_depth, &va_talus, &va_melting},
      [](auto, auto out, const auto &, auto &)
      {
        *out[0] = Array(out[0]->shape, 1.0f);
        *out[1] = Array(out[1]->shape, 0.5f);
        *out[2] = Array(out[2]->shape, 0.0f);
      },
      cm);

  std::cout << "Running snow_simulation..." << std::endl;
  for_each_tile(
      {&va_z, &va_depth, &va_melting, &va_talus},
      {&va_snow},
      [&](std::vector<const Array *> p_arrays_in,
          std::vector<Array *>       p_arrays_out,
          const TileRegion          &region,
          ICommunicator             &comm)
      {
        std::cout << "In tile (" << region.key.tx << ", " << region.key.ty << ")" << std::endl;
        auto [pa_z, pa_depth_map, pa_melting_map, pa_talus_map] =
            unpack<4>(p_arrays_in);
        auto [pa_snow_depth] = unpack<1>(p_arrays_out);

        const int sync_interval = 10;
        *pa_snow_depth =
            gpu::snow_simulation(*pa_z,
                                 10.0f,
                                 *pa_depth_map,
                                 *pa_melting_map,
                                 *pa_talus_map,
                                 30, // 30 iterations
                                 0.5f,
                                 1.0f,
                                 0.1f,
                                 0.1f,
                                 0.1f,
                                 0.1f,
                                 0.1f,
                                 0.1f,
                                 true,
                                 0.5f,
                                 true,
                                 0.0f,
                                 true,
                                 sync_interval,
                                 region,
                                 comm);
        std::cout << "Tile finished (" << region.key.tx << ", " << region.key.ty << ")" << std::endl;
      },
      cm);

  std::cout << "Done!" << std::endl;
}
