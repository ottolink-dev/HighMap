#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "highmap.hpp"

#define EXPECT(cond, msg)                                                      \
  do                                                                           \
  {                                                                            \
    if (!(cond))                                                               \
      std::cerr << "[FAIL] " << msg << "\n";                                   \
    else                                                                       \
      std::cout << "[ OK ] " << msg << "\n";                                   \
  } while (0)

int main()
{
  hmap::gpu::init_opencl();

  // ===========================================================================
  // Configuration
  // ===========================================================================

  const glm::vec4 bbox{0.f, 1.f, 0.f, 1.f};

  const glm::ivec2 shape{1024, 1024};
  const glm::ivec2 tile_shape{220, 312};
  const size_t     halo = 32;

  // const auto storage_mode = hmap::StorageMode::VA_RAM;
  const auto storage_mode = hmap::StorageMode::VA_DISK_LRU;
  // const auto storage_mode = hmap::StorageMode::VA_DISK_LRU_MIN;
  // const auto storage_mode = hmap::StorageMode::VA_DISK_SEQUENTIAL;

  const auto foreach_mode = hmap::ForEachMode::VA_DISTRIBUTED;
  // const auto foreach_mode = hmap::ForEachMode::VA_SEQUENTIAL;
  // const auto foreach_mode = hmap::ForEachMode::VA_SINGLE_ARRAY;
  // const auto foreach_mode = hmap::ForEachMode::VA_SINGLE_ARRAY_DOWNSCALED;

  const hmap::ComputeMode cm{.mode = foreach_mode,
                             .trim_storage = true,
                             .stride = 1};

  // ===========================================================================
  // Virtual arrays
  // ===========================================================================

  hmap::VirtualArray varray(shape, bbox, tile_shape, halo, storage_mode);
  hmap::VirtualArray varray_binary(shape, bbox, tile_shape, halo, storage_mode);

  // ===========================================================================
  // Tile operations
  // ===========================================================================

  const auto generate_noise =
      [](hmap::Array &tile, const hmap::TileRegion &region)
  {
    tile = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                           region.shape,
                           {2.f, 2.f},
                           0,
                           8,
                           0.7f,
                           0.5f,
                           2.f,
                           nullptr,
                           nullptr,
                           nullptr,
                           region.bbox);
  };

  const auto apply_gamma = [](hmap::Array &tile, const hmap::TileRegion &)
  { hmap::gamma_correction(tile, 2.f); };

  const auto smooth = [](hmap::Array &tile, const hmap::TileRegion &)
  { hmap::smooth_cpulse(tile, 64); };

  hmap::for_each_tile(varray, generate_noise, cm);
  hmap::for_each_tile(varray, apply_gamma, cm);
  // hmap::for_each_tile(varray, smooth, cm);

  // specific
  {
    const auto erosion = [](hmap::Array &tile, const hmap::TileRegion &)
    { hmap::gpu::hydraulic_stream_log(tile, 0.2f, 1e0f); };

    hmap::ComputeMode cm_dwn = {
        .mode = hmap::ForEachMode::VA_SINGLE_ARRAY_DOWNSCALED,
        .trim_storage = false,
        .k_cutoff = 0.3f};

    hmap::for_each_tile(varray, erosion, cm_dwn);
  }

  const hmap::ComputeMode cm_local = {.mode = hmap::ForEachMode::VA_SEQUENTIAL,
                                      .trim_storage = false};

  varray.to_array(cm_local).to_png("out_noise.png", hmap::Cmap::JET);

  // ===========================================================================
  // Multi-array operation
  // ===========================================================================

  const auto binarize_copy =
      [](const std::vector<const hmap::Array *> &in,
         std::vector<hmap::Array *>             &out,
         const hmap::TileRegion                 &)
  {
    auto &src = *in[0];
    auto &dst = *out[0];

    dst = src;
    hmap::make_binary(dst, 0.5f);
  };

  varray.trim_storage();

  hmap::for_each_tile({&varray}, {&varray_binary}, binarize_copy, cm);

  varray_binary.to_array(cm_local).to_png("out_binary.png", hmap::Cmap::JET);

  // ===========================================================================
  // Global operations
  // ===========================================================================

  varray.sync_overlap_buffers();
  varray.remap(0.f, 1.f, cm);

  // ===========================================================================
  // Tests
  // ===========================================================================

  // --- Determinism

  {
    auto a0 = varray.to_array(cm);
    auto a1 = varray.to_array(cm);

    float max_diff = 0.f;
    for (int j = 0; j < a0.shape.y; ++j)
      for (int i = 0; i < a0.shape.x; ++i)
        max_diff = std::max(max_diff, std::abs(a0(i, j) - a1(i, j)));

    EXPECT(max_diff == 0.f, "deterministic to_array()");
  }

  // --- Range invariant

  {
    float min_v = varray.min(cm);
    float max_v = varray.max(cm);

    EXPECT(min_v >= 0.f, "min >= 0 after remap");
    EXPECT(max_v <= 1.f, "max <= 1 after remap");

    glm::vec2 range = varray.range(cm);

    EXPECT(min_v == range.x, "min_v == range.x");
    EXPECT(max_v == range.y, "max_v == range.y");
  }

  {
    float low_p = 0.05f;
    float high_p = 0.95f;

    glm::vec2 range_p = varray.range_percentile(low_p, high_p, cm);
    glm::vec2 range_p_ref = varray.to_array_dbg().range_percentile(low_p,
                                                                   high_p);

    EXPECT(range_p.x == range_p_ref.x, "range_p.x == range_p_ref.x");
    EXPECT(range_p.y == range_p_ref.y, "range_p.y == range_p_ref.y");
  }

  // --- Binary validity

  {
    auto values = varray_binary.unique_values(cm);

    bool ok = true;
    for (float v : values)
      ok &= (v == 0.f || v == 1.f);

    EXPECT(ok, "binary array contains only {0,1}");
  }

  // --- Reduction consistency

  {
    auto arr = varray.to_array(cm);

    double sum = 0.0;
    for (int j = 0; j < arr.shape.y; ++j)
      for (int i = 0; i < arr.shape.x; ++i)
        sum += arr(i, j);

    double mean_ref = sum / (arr.shape.x * arr.shape.y);
    double mean_va = varray.mean(cm);

    EXPECT(std::abs(mean_ref - mean_va) < 1e-6f, "mean == sum / N");
  }

  // --- from_array round-trip

  {
    hmap::Array src = hmap::noise_fbm(hmap::NoiseType::PERLIN,
                                      shape,
                                      {2.f, 2.f},
                                      0);

    hmap::VirtualArray v(shape, bbox, tile_shape, halo, storage_mode);
    v.from_array(src, cm);

    auto dst = v.to_array(cm);

    float max_diff = 0.f;
    for (int j = 0; j < shape.y; ++j)
      for (int i = 0; i < shape.x; ++i)
        max_diff = std::max(max_diff, std::abs(src(i, j) - dst(i, j)));

    EXPECT(max_diff < 1e-6f, "from_array() round-trip");
  }

  // --- Interpolation sanity

  {
    float v0 = varray.get_nearest(bbox.x, bbox.z);
    float v1 = varray.get_bilinear(bbox.y, bbox.w);

    EXPECT(std::isfinite(v0), "nearest interpolation finite");
    EXPECT(std::isfinite(v1), "bilinear interpolation finite");
  }

  // ===========================================================================
  // Small tile stress test (edge case)
  // ===========================================================================

  {
    const glm::ivec2 small_tile{8, 8};
    const size_t     big_halo = 7;

    hmap::VirtualArray v_small(shape, bbox, small_tile, big_halo, storage_mode);

    hmap::for_each_tile(v_small, generate_noise, cm);
    v_small.sync_overlap_buffers();
    v_small.remap(0.f, 1.f, cm);

    float min_v = v_small.min(cm);
    float max_v = v_small.max(cm);

    EXPECT(min_v >= 0.f, "small tile min >= 0");
    EXPECT(max_v <= 1.f, "small tile max <= 1");

    auto arr = v_small.to_array(cm);
    EXPECT(arr.shape == shape, "small tile shape preserved");
  }

  return 0;
}
