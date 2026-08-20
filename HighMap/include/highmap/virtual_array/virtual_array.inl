/* Copyright (c) 2025 Otto Link. Distributed under the terms of the GNU General
   Public License. The full license is in the file LICENSE, distributed with
   this software. */
#pragma once
#if defined(__APPLE__)
#include <sys/sysctl.h>
#endif

struct TileAccess
{
  std::vector<const VirtualArray *> inputs;
  std::vector<VirtualArray *>       outputs;

  const VirtualArray *ref_va() const
  {
    return this->outputs.front()
               ? this->outputs.front()
               : (!this->inputs.empty() ? this->inputs.front() : nullptr);
  }
};

template <typename Func>
void for_each_tile(const TileAccess &access, Func &&func, const ComputeMode &cm)
{
  if (access.outputs.empty())
  {
    LOG_ERROR("for_each_tile: no output VirtualArray");
    return;
  }

  const VirtualArray *ref_va = access.ref_va();
  if (!ref_va) return;

  // --- validate compatibility

  auto check = [&](const VirtualArray *va)
  {
    if (!va) return true;
    return va->shape == ref_va->shape && va->tile_shape == ref_va->tile_shape;
  };

  for (auto *va : access.inputs)
    if (!check(va))
    {
      LOG_ERROR("Incompatible input VirtualArray");
      return;
    }

  for (auto *va : access.outputs)
    if (!check(va))
    {
      LOG_ERROR("Incompatible output VirtualArray");
      return;
    }

  // --- region dispatcher for tiled computations

  auto region_dispatcher = [&](const TileRegion &region)
  {
    std::vector<const Array *> in_tiles;
    std::vector<Array *>       out_tiles;

    in_tiles.reserve(access.inputs.size());
    out_tiles.reserve(access.outputs.size());

    for (auto *va : access.inputs)
      in_tiles.push_back(va ? &va->storage->get_tile(region) : nullptr);

    for (auto *va : access.outputs)
      out_tiles.push_back(va ? &va->storage->get_tile(region) : nullptr);

    func(in_tiles, out_tiles, region);

    for (auto *va : access.inputs)
      if (va) va->storage->release_tile(region);

    for (auto *va : access.outputs)
      if (va) va->storage->release_tile(region);
  };

  // reuse your existing modes
  switch (cm.mode)
  {
  case ForEachMode::VA_SEQUENTIAL:
    sequential_tile_loop(*ref_va, region_dispatcher);
    break;
  case ForEachMode::VA_DISTRIBUTED:
    distributed_tile_loop(*ref_va, region_dispatcher);
    break;
  case ForEachMode::VA_SINGLE_ARRAY:
    single_array_compute(access, func, cm.stride);
    break;
  case ForEachMode::VA_SINGLE_ARRAY_DOWNSCALED:
    single_array_downscaled_compute(access, func, cm.k_cutoff);
    break;
  }

  if (cm.trim_storage)
  {
    for (auto *va : access.inputs)
      if (va) va->trim_storage();
    for (auto *va : access.outputs)
      if (va) va->trim_storage();
  }
}

// fake const...
template <typename Func>
void for_each_tile(const VirtualArray &va, Func &&func, const ComputeMode &cm)
{
  auto &mutable_va = const_cast<VirtualArray &>(va);

  TileAccess acc;
  acc.outputs = {&mutable_va};

  for_each_tile(
      acc,
      [&](const std::vector<const Array *> &,
          std::vector<Array *> &out,
          const TileRegion     &region) { func(*out[0], region); },
      cm);
}

// no-const
template <typename Func>
void for_each_tile(VirtualArray &va, Func &&func, const ComputeMode &cm)
{
  TileAccess acc;
  acc.outputs = {&va};

  for_each_tile(
      acc,
      [&](const std::vector<const Array *> &,
          std::vector<Array *> &out,
          const TileRegion     &region) { func(*out[0], region); },
      cm);
}

template <typename Func>
void for_each_tile(const std::vector<VirtualArray *> &outputs,
                   Func                             &&func,
                   const ComputeMode                 &cm)
{
  TileAccess acc;
  acc.outputs = outputs;

  for_each_tile(
      acc,
      [&](const std::vector<const Array *> &,
          std::vector<Array *> &out,
          const TileRegion     &region) { func(out, region); },
      cm);
}

template <typename Func>
void for_each_tile(const std::vector<const VirtualArray *> &inputs,
                   const std::vector<VirtualArray *>       &outputs,
                   Func                                   &&func,
                   const ComputeMode                       &cm)
{
  TileAccess acc;
  acc.inputs = inputs;
  acc.outputs = outputs;

  for_each_tile(
      acc,
      [&](const std::vector<const Array *> &in,
          std::vector<Array *>             &out,
          const TileRegion                 &region) { func(in, out, region); },
      cm);
}

// --- COMPUTE

template <typename RegionDispatcher>
void sequential_tile_loop(const VirtualArray &ref_va,
                          RegionDispatcher  &&dispatcher)
{
  const int nx = ceil_div(ref_va.shape.x, ref_va.tile_shape.x);
  const int ny = ceil_div(ref_va.shape.y, ref_va.tile_shape.y);

  for (int ty = 0; ty < ny; ++ty)
    for (int tx = 0; tx < nx; ++tx)
    {
      TileRegion region = ref_va.tile_region_from_tile_coords(tx, ty);
      dispatcher(region);
    }
}

/// Return a sensible default worker count for tile-parallel computations.
///
/// On Apple Silicon (heterogeneous P/E cores) the pool is sized after the
/// performance cores: long per-tile tasks scheduled on efficiency cores
/// become stragglers that delay the whole batch. Other platforms fall back
/// to the full hardware concurrency.
inline int recommended_thread_count()
{
#if defined(__APPLE__)
  int    nperf = 0;
  size_t len = sizeof(nperf);
  if (sysctlbyname("hw.perflevel0.physicalcpu", &nperf, &len, nullptr, 0) ==
          0 &&
      nperf > 0)
    return nperf;
#endif
  return int(std::thread::hardware_concurrency());
}

template <typename RegionDispatcher>
void distributed_tile_loop(const VirtualArray &ref_va,
                           RegionDispatcher  &&dispatcher,
                           int                 nthreads = 0)
{
  const int nx = ceil_div(ref_va.shape.x, ref_va.tile_shape.x);
  const int ny = ceil_div(ref_va.shape.y, ref_va.tile_shape.y);
  const int ntasks = nx * ny;

  if (nthreads <= 0) nthreads = recommended_thread_count();

  nthreads = std::min(nthreads, ntasks);

  // dynamic scheduling: workers pull the next tile through a shared atomic
  // counter. This balances uneven per-tile costs (data-dependent erosion)
  // and keeps heterogeneous cores (e.g. Apple Silicon P/E clusters) busy
  // until the last task, unlike the previous static strided partitioning
  std::atomic<int> next_task{0};

  auto worker = [&]()
  {
    while (true)
    {
      const int k = next_task.fetch_add(1, std::memory_order_relaxed);
      if (k >= ntasks) break;

      const int ty = k / nx;
      const int tx = k % nx;

      TileRegion region = ref_va.tile_region_from_tile_coords(tx, ty);
      dispatcher(region);
    }
  };

  std::vector<std::future<void>> futures;
  futures.reserve(nthreads);

  for (int t = 0; t < nthreads; ++t)
    futures.emplace_back(std::async(std::launch::async, worker));

  for (auto &f : futures)
    f.get();
}

template <typename Func>
void single_array_compute(const TileAccess &access, Func &&func, int stride = 1)
{
  const VirtualArray *ref_va = access.ref_va();
  if (!ref_va) return;

  // create temporary arrays
  std::vector<Array> arrays_in;
  std::vector<Array> arrays_out;
  arrays_in.reserve(access.inputs.size());
  arrays_out.reserve(access.outputs.size());

  // force a CPU mode for data gathering operations
  const ComputeMode cm_local = {.mode = ForEachMode::VA_DISTRIBUTED,
                                .trim_storage = false};

  glm::ivec2 shape_wrk = ref_va->shape;

  if (stride > 1)
  {
    int nx = std::max(2, int(ref_va->shape.x / float(stride)));
    int ny = std::max(2, int(ref_va->shape.y / float(stride)));
    shape_wrk = glm::ivec2(nx, ny);
  }

  // project
  for (auto *va : access.inputs)
    arrays_in.push_back(va ? va->to_array(cm_local).resample_to_shape(shape_wrk)
                           : Array());

  for (auto *va : access.outputs)
    // TODO remove to_array(), left for retor-fit but misusage
    arrays_out.push_back(
        va ? va->to_array(cm_local).resample_to_shape(shape_wrk) : Array());

  // create pointer vector for user func
  std::vector<const Array *> in_tiles;
  std::vector<Array *>       out_tiles;
  in_tiles.reserve(access.inputs.size());
  out_tiles.reserve(access.outputs.size());

  for (size_t i = 0; i < access.inputs.size(); ++i)
    in_tiles.push_back(access.inputs[i] ? &arrays_in[i] : nullptr);

  for (size_t i = 0; i < access.outputs.size(); ++i)
    out_tiles.push_back(access.outputs[i] ? &arrays_out[i] : nullptr);

  // call user operation
  const TileRegion va_region = TileRegion(TileKey(),
                                          ref_va->bbox,
                                          shape_wrk,
                                          {0, 0, 0, 0});

  func(in_tiles, out_tiles, va_region);

  // copy back to VirtualArrays (only outputs do not need to update const
  // inputs)
  for (size_t i = 0; i < access.outputs.size(); ++i)
    if (access.outputs[i])
    {
      if (shape_wrk != ref_va->shape)
        access.outputs[i]->from_array_bicubic(arrays_out[i], cm_local);
      else
        access.outputs[i]->from_array(arrays_out[i], cm_local);
    }
}

template <typename Func>
void single_array_downscaled_compute(const TileAccess &access,
                                     Func            &&func,
                                     float             kc)
{
  const VirtualArray *ref_va = access.ref_va();
  if (!ref_va) return;

  // create temporary arrays
  std::vector<Array> arrays_in;
  std::vector<Array> arrays_out;
  arrays_in.reserve(access.inputs.size());
  arrays_out.reserve(access.outputs.size());

  // force a CPU mode for data gathering operations
  const ComputeMode cm_local = {.mode = ForEachMode::VA_DISTRIBUTED,
                                .trim_storage = false};

  glm::ivec2 shape_wrk = ref_va->shape;

  if (kc != 1.f)
  {
    int nx = std::max(2, int(ref_va->shape.x * kc));
    int ny = std::max(2, int(ref_va->shape.y * kc));
    shape_wrk = glm::ivec2(nx, ny);
  }

  // project
  for (auto *va : access.inputs)
    arrays_in.push_back(va ? va->to_array(cm_local).resample_to_shape(shape_wrk)
                           : Array());

  for (auto *va : access.outputs)
    // TODO remove to_array(), left for retor-fit but misusage
    arrays_out.push_back(
        va ? va->to_array(cm_local).resample_to_shape(shape_wrk) : Array());

  // create pointer vector for user func
  std::vector<const Array *> in_tiles;
  std::vector<Array *>       out_tiles;
  in_tiles.reserve(access.inputs.size());
  out_tiles.reserve(access.outputs.size());

  for (size_t i = 0; i < access.inputs.size(); ++i)
    in_tiles.push_back(access.inputs[i] ? &arrays_in[i] : nullptr);

  for (size_t i = 0; i < access.outputs.size(); ++i)
    out_tiles.push_back(access.outputs[i] ? &arrays_out[i] : nullptr);

  // call user operation
  const TileRegion va_region = TileRegion(TileKey(),
                                          ref_va->bbox,
                                          shape_wrk,
                                          {0, 0, 0, 0});

  func(in_tiles, out_tiles, va_region);

  // copy back to VirtualArrays (only outputs do not need to update const
  // inputs)
  for (size_t i = 0; i < access.outputs.size(); ++i)
    if (access.outputs[i])
    {
      // resample to initial shape and add higher frequency
      Array out_full = access.outputs[i]->to_array(cm_local);
      Array out_low_pass = out_full.resample_to_shape(shape_wrk)
                               .resample_to_shape_bicubic(ref_va->shape);

      // resample coarse result and use it to replace the initial
      // low-frequency content
      arrays_out[i] = arrays_out[i].resample_to_shape_bicubic(ref_va->shape);

      // TODO add post-filter
      // laplace(arrays_out[i]);

      access.outputs[i]->from_array(arrays_out[i] + (out_full - out_low_pass),
                                    cm_local);
    }
}

template <typename Func>
void for_each_cell(VirtualArray &va, Func &&func, const ComputeMode &cm)
{
  for_each_tile(
      va,
      [&](Array &tile, const TileRegion &region)
      {
        for (int j = 0; j < region.shape.y; ++j)
          for (int i = 0; i < region.shape.x; ++i)
            func(tile(i, j), i, j, region);
      },
      cm);
}
