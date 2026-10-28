# BVH2D Spatial Indexing for Path SDF Computation

This document records the performance investigation and implementation findings for accelerating 2D Path Signed Distance Field (SDF) evaluation using a Bounding Volume Hierarchy (`BVH2D`).

## 1. Context and Motivation

Computing the distance field of a 2D polyline or polygon across a discrete heightmap grid is a recurring operation in terrain authoring, hydrology (riverbed carve profiles), and carving workflows.

The original evaluation in `path_sdf_to_array` and `sdf_2d_polyline` performed brute-force distance evaluation:
- For every pixel $(i, j)$ in a $W \times H$ array, iterate across all $N$ path segments.
- Complexity: $\mathcal{O}(W \cdot H \cdot N)$.

For dense or detailed paths ($N \ge 100$) and high-resolution heightmaps ($512 \times 512$ or higher), this approach causes significant CPU bottlenecks.

HighMap provides an internal 2D Axis-Aligned Bounding Box tree (`hmap::BVH2D`). By indexing the segment set into a BVH, the nearest segment search per pixel is accelerated via branch-and-bound pruning.

## 2. Algorithmic Approach

### Segment Indexing
Path segments $[p_k, p_{k+1}]$ are packaged into `Segment2D` structs storing directional vector $d = p_{k+1} - p_k$ and inverse squared length $\frac{1}{\|d\|^2}$.
A `BVH2D` tree is built over the segments in $\mathcal{O}(N \log N)$ time.

### Branch-and-Bound Distance Query
For a given query point $p = (x, y)$:
1. Traverse the BVH using a fixed-capacity traversal stack.
2. At each node, compute the minimum squared distance from $p$ to the node's axis-aligned bounding box (`AABB2D`):
   $$d_{\text{box}}^2(p, \text{box}) = \max(0, \text{box.min}_x - p_x, p_x - \text{box.max}_x)^2 + \max(0, \text{box.min}_y - p_y, p_y - \text{box.max}_y)^2$$
3. If $d_{\text{box}}^2 \ge d_{\text{best}}^2$, the node and all its descendants are pruned.
4. For internal nodes, visit the closer child first to tighten $d_{\text{best}}^2$ quickly.
5. In leaf nodes, evaluate exact point-to-segment distance squared and update $d_{\text{best}}^2$.

Average per-query complexity drops to $\mathcal{O}(\log N)$, reducing overall runtime to $\mathcal{O}(W \cdot H \log N)$.

### Support for Domain Warping and Path Closure
- **Domain Warping**: Per-pixel noise displacements (`p_noise_x`, `p_noise_y`) warp coordinate $p = (x + dx, y + dy)$ before querying the tree. The tree does not require reconstruction.
- **Closed Paths**: When `path.is_closed() == true`, the closing segment is included in the BVH, and an orientation/winding sign test determines interior vs. exterior distances.

## 3. Benchmark Measurements

Google Benchmark was executed comparing the baseline brute-force approach against the BVH approach on a CPU with 6 physical cores (12 threads).

### Execution Time Comparison ($512 \times 512$ Grid)

| Segment Count ($N$) | Baseline Brute Force (ms) | BVH2D Accelerated (ms) | Speedup Factor |
|---|---|---|---|
| 10 | 39.8 | 6.71 | ~5.9x |
| 50 | 205.0 | 12.20 | ~16.8x |
| 200 | 734.0 | 20.20 | ~36.3x |
| 1000 | 2823.0 | 37.60 | ~75.1x |

### Execution Time Comparison ($256 \times 256$ Grid)

| Segment Count ($N$) | Baseline Brute Force (ms) | BVH2D Accelerated (ms) | Speedup Factor |
|---|---|---|---|
| 10 | 9.59 | 1.42 | ~6.8x |
| 50 | 51.00 | 4.37 | ~11.7x |
| 200 | 191.00 | 4.84 | ~39.5x |
| 1000 | 917.00 | 10.30 | ~89.0x |

## 4. API Design and Integration

The functionality is structured to preserve backward compatibility while providing the brute-force implementation for reference and validation:

1. **Explicit Implementations**:
   - `hmap::sdf_2d_polyline_bvh` and `hmap::path_sdf_to_array_bvh`
   - `hmap::sdf_2d_polyline_brute_force` and `hmap::path_sdf_to_array_brute_force`

2. **Default Wrapper Functions**:
   - `hmap::sdf_2d_polyline(path, shape, bbox, p_noise_x, p_noise_y, use_bvh = true)`
   - `hmap::path_sdf_to_array(path, shape, bbox_array, p_noise_x, p_noise_y, use_bvh = true)`

Existing call sites retain default `use_bvh = true` behavior without signature breakage.
