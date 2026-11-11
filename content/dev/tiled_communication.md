# Tiled Communication and Overlap Synchronization

This document explains the architecture of the tiled communication system in HighMap, how collective directives operate across execution modes, and how halo overlap buffers are synchronized during iterative simulations.

## 1. Overview and Purpose

Terrain processing algorithms operating on large spatial domains must frequently scale beyond available GPU VRAM or system RAM. HighMap achieves out-of-core scalability through the `VirtualArray` container, which subdivides heightmaps into a regular grid of overlapping tiles.

While local kernel operations (e.g. stencil filters with radius smaller than halo) require only local neighborhood information, global and iterative algorithms require communication across tile boundaries:

- Global reductions (computing global minimum, maximum, volume, or convergence residuals).
- Global unions (collecting feature points, active wavefronts, or sink lists).
- Execution synchronization barriers.
- Iterative boundary halo exchanges (diffusion, wave propagation, flooding, snow creep, viscous flow, Poisson solvers).

HighMap provides an abstract communicator interface, `ICommunicator`, decoupling simulation code from the underlying concurrency model.

## 2. Communicator Interface (`ICommunicator`)

The `ICommunicator` interface provides an MPI-like abstraction tailored to 2D raster tile topologies:

```cpp
class ICommunicator
{
public:
  virtual ~ICommunicator() = default;

  virtual float              all_reduce(float local_val, ReduceOp op) = 0;
  virtual std::vector<float> all_union(const std::vector<float> &local_vals) = 0;
  virtual void               barrier() = 0;
  virtual int                rank() const = 0;
  virtual int                size() const = 0;
  virtual void               sync_overlap_buffers(
                                Array            &local_tile,
                                const TileRegion &region,
                                SyncOperation     op = SyncOperation::SmoothBlend) = 0;
};
```

### Supported Operations

- **`all_reduce(local_val, op)`**: Aggregates a scalar across all tiles using `ReduceOp::Min`, `ReduceOp::Max`, or `ReduceOp::Sum`. Returns the identical aggregated value to all calling tiles.
- **`all_union(local_vals)`**: Gathers and combines lists of scalar items (such as active seeds or coordinates) across all tiles.
- **`barrier()`**: Blocks execution until all tiles in the current collective phase reach the synchronization point.
- **`sync_overlap_buffers(local_tile, region, op)`**: Synchronizes seam buffers along tile overlaps using a selected `SyncOperation`.

### Overlap Synchronization Modes (`SyncOperation`)

- **`SmoothBlend`**: Blends overlapping cells smoothly using a quintic polynomial ($C^2$ smoothstep). Used by diffusion, flow simulations, snow settling, and Poisson boundary relaxations to prevent sharp seam artifacts.
- **`Average` / `Mean`**: Computes the arithmetic mean between corresponding overlap cells.
- **`Min`**: Takes the element-wise minimum across overlapping regions. Used by distance transforms and geodesic propagation.
- **`Max`**: Takes the element-wise maximum across overlapping regions. Used by priority-flood depression filling and water elevation surfaces.
- **`CopyFirst` / `CopySecond`**: Directional copy from one tile to another along boundary interfaces.

## 3. Execution Engines

The communicator adapts dynamically to the compute configuration specified in `ComputeMode`:

### Single-Array Mode (`SingleArrayComm`)

When an algorithm is invoked directly on a monolithic `Array` without tiling, `SingleArrayComm::instance()` serves as the default communicator:
- `size()` returns `1`, and `rank()` returns `0`.
- `all_reduce` returns the input scalar directly with zero overhead.
- `all_union` returns the input vector unchanged.
- `barrier()` and `sync_overlap_buffers()` are immediate no-ops.

Algorithms written with communicator parameters incur zero memory or synchronization penalty when run monolithically.

### Multithreaded Distributed Mode (`ThreadedTileComm`)

In `ForEachMode::VA_DISTRIBUTED` when sufficient memory exists to hold all tiles simultaneously, tiles execute concurrently across a thread pool:

1. **Registration**: Every thread registers its tile pointer and `TileRegion` in a shared `ThreadedCollectiveState`.
2. **Barrier**: Threads enter a condition-variable barrier waiting for all ranks to register.
3. **Seam Exchange**: Rank 0 iterates through the horizontal ($x$) and vertical ($y$) tile seams, reading and writing the overlapping memory regions between neighboring `Array` buffers according to the requested `SyncOperation`.
4. **Resumption**: A second barrier releases all threads once seam synchronization is complete.

### Out-of-Core Sequential Mode (`SequentialTileComm`)

When memory constraints dictate strict out-of-core execution (e.g. `StorageMode::VA_DISK_LRU_MIN`, where at most 2 live tiles are kept in memory at any instant), threads cannot execute in parallel.

HighMap implements a cooperative sequential scheduler (`SequentialThreadScheduler`):

1. **Tile Fibers**: Each tile runs inside a cooperative thread. Exactly one tile executes at any time.
2. **Staging Cache**: When a tile calls `sync_overlap_buffers`, its updated contents are flushed to a temporary disk-backed `VirtualArray` (`sync_va`).
3. **Turn Passing**: The active tile suspends, freeing its RAM tile allocation, and yields execution to the next rank in round-robin order.
4. **Global Seam Sweep**: Once all tiles reach the synchronization point, the scheduler executes `sync_va->sync_overlap_buffers(op)` using the low-memory 2-tile sliding window cache.
5. **Tile Resumption**: Ranks wake up in sequence and reload the updated, synchronized tile data from `sync_va`.

This mechanism guarantees strict adherence to low memory ceilings (e.g. maximum 2 live tiles) while preserving multi-pass synchronization consistency across the entire grid.

## 4. How Overlap Buffering Works

Each tile in a `VirtualArray` consists of a core nominal area augmented by a border of width $H$ (the halo):

```text
+-------------------------------------------------------+
|                 Top Halo (North)                      |
|   +-----------------------------------------------+   |
| L |                                               | R |
| e |                                               | i |
| f |                  Core Tile                    | g |
| t |             (Nominal Domain)                  | h |
|   |                                               | t |
| H |                                               |   |
| a |                                               | H |
| l |                                               | a |
| o |                                               | l |
|   +-----------------------------------------------+   |
|                Bottom Halo (South)                    |
+-------------------------------------------------------+
```

When two adjacent tiles $T_0$ and $T_1$ are placed side-by-side along the $x$-axis:
- The right halo of $T_0$ covers the left core edge of $T_1$.
- The left halo of $T_1$ covers the right core edge of $T_0$.
- The total overlap width between $T_0$ and $T_1$ is $2H$.

During `SmoothBlend` seam stitching, cells across the overlap $[0, H)$ are blended smoothly:

$$r = \text{smoothstep5}\left(\frac{p}{H - 1}\right)$$

$$\text{value}(p) = (1 - r) \cdot T_0(p_{\text{buf}}) + r \cdot T_1(p)$$

The resulting blended values are written to both tiles, ensuring continuous derivatives across the seam.

## 5. Iteration Batching in Domain Algorithms

Iterative physical simulations (such as hydraulic erosion, shallow water flow, snow creep, or Poisson pressure solvers) run for tens to hundreds of iterations. Synchronizing seams after every individual step would incur prohibitive synchronization overhead.

HighMap algorithms expose a `sync_interval` parameter enabling **iteration batching**:

```text
Total Iterations: [================================================]
Batch 0:          [==========] -> sync_overlap_buffers
Batch 1:                     [==========] -> sync_overlap_buffers
Batch 2:                                [==========] -> sync_overlap_buffers
Final:                                              [==========] -> sync_overlap_buffers
```

### Batching Loop Pattern

Algorithms follow a standard batching structure:

```cpp
const int batch_size = (comm.size() > 1 && sync_interval > 0)
                           ? sync_interval
                           : iterations;

for (int it_start = 0; it_start < iterations; it_start += batch_size)
{
  const int it_end = std::min(iterations, it_start + batch_size);

  {
    // 1. Scope GPU kernel execution for [it_start, it_end)
    auto run = clwrapper::Run("kernel_name");
    ...
    for (int it = it_start; it < it_end; ++it)
    {
      run.execute_async(...);
    }
    run.finish();

    // 2. Read back device buffers to host Array
    run.read_imagef("out");
  }

  // 3. Synchronize halo boundaries between batches
  if (it_end < iterations)
  {
    comm.sync_overlap_buffers(local_tile, region, SyncOperation::SmoothBlend);
  }
}

// 4. Final boundary synchronization
comm.sync_overlap_buffers(local_tile, region, SyncOperation::SmoothBlend);
```

By scoping OpenCL `Run` objects to each batch, device resources and command queues are flushed and synchronized with host memory before collective communicator operations execute.

## 6. Integration and Usage Example

To write a filter or simulation supporting both monolithic and tiled execution:

```cpp
void my_simulation(Array            &z,
                   int               iterations,
                   int               sync_interval = 0,
                   const TileRegion &region = TileRegion::dummy(),
                   ICommunicator    &comm = SingleArrayComm::instance())
{
  if (!validate_non_empty(z)) return;

  const int batch_size = (comm.size() > 1 && sync_interval > 0)
                             ? sync_interval
                             : iterations;

  for (int it_start = 0; it_start < iterations; it_start += batch_size)
  {
    const int it_end = std::min(iterations, it_start + batch_size);

    for (int it = it_start; it < it_end; ++it)
    {
      // perform local simulation step on z
    }

    if (it_end < iterations)
    {
      comm.sync_overlap_buffers(z, region, SyncOperation::SmoothBlend);
    }
  }

  comm.sync_overlap_buffers(z, region, SyncOperation::SmoothBlend);
}
```

Invoking through `for_each_tile` automatically injects the active communicator and region metadata:

```cpp
ComputeMode cm{.mode = ForEachMode::VA_DISTRIBUTED};

for_each_tile(
    va,
    [](Array &tile, const TileRegion &region, ICommunicator &comm)
    {
      my_simulation(tile, 100, 10, region, comm);
    },
    cm);
```

## 7. Buffer Shape Invariance and Nested Multiscale Pipelines

A fundamental invariant of the tiled communication model is that any `Array` buffer passed to `sync_overlap_buffers(local_tile, region, op)` must match the full tile geometry described by `region.shape` (nominal tile size plus halo margins).

### Incompatibility with Nested Multiscale Pyramids

Algorithms that construct internal resolution pyramids (such as multiscale filters or cascades that downsample local tiles to coarse ladder shapes $L_i$) cannot propagate the outer `TileRegion` or `ICommunicator` into the coarse levels:

- The downsampled buffer $z_{\text{coarse}}$ possesses dimensions $L_i$, whereas `region.shape` and `region.halo` correspond to the full tile resolution.
- Seam indexing logic assumes the halo margins and stride match `region.shape`. If invoked on a downsampled buffer, halo offsets exceed the array boundaries, leading to heap corruption.
- Intermediate downsampled levels represent purely local approximations used to calculate local delta adjustments (e.g. coarse valley carving), and are not decomposed into valid global sub-domains.

In such multiscale pipelines:
- Intermediate/coarse levels must run **locally** using `TileRegion::dummy()` and `SingleArrayComm::instance()`.
- The active `TileRegion` and `ICommunicator` must only be passed to the final resolution level where the buffer dimensions match `region.shape`.

### Runtime Shape Validation

Both `ThreadedTileComm` and `SequentialTileComm` validate that `local_tile.shape == region.shape` at the beginning of `sync_overlap_buffers`. If a mismatch is detected, an error is logged and synchronization is skipped to safeguard memory integrity:

```text
sync_overlap_buffers: buffer shape (x, y) does not match region shape (rx, ry).
This may happen if sync is invoked on an internal or nested downsampled buffer (e.g. multiscale pyramid).
Skipping synchronization to prevent memory corruption.
```
