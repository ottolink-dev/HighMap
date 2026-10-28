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

## 5. Verification of Bézier Polyline SDF (`sdf_2d_polyline_bezier`)

During verification of `sdf_2d_polyline_bezier`, numerical issues were discovered causing corrupted/broken results (entire arrays returning `NaN` or `inf`).

### Root Cause Analysis

1. **Division by Zero on Colinear/Degenerate Segments**:
   The distance to a quadratic Bézier curve $B(t) = (1-t)^2 P_0 + 2t(1-t)P_1 + t^2 P_2$ involves the second-difference vector:
   $$b = P_0 - 2P_1 + P_2$$
   When the control points are colinear (or $P_1$ is at the midpoint of $P_0$ and $P_2$), $b \approx 0 \implies \text{dot}(b, b) = 0$.
   The original code computed:
   ```cpp
   float kk = 1.0f / dot(b, b);
   ```
   without checking for zero. This produced `inf`, which cascaded through Cardano's cubic solver into `NaN` for all grid queries.
   **Fix**: When $\text{dot}(b, b) < 10^{-10}$, the curve degenerates to a straight line segment $[P_0, P_2]$. The solver directly evaluates the point-to-line-segment distance.

2. **Domain Violations in Trigonometric Cardano Branch**:
   For the 3-real-roots case ($h < 0$ in the cubic discriminant), the solution evaluates:
   ```cpp
   float v = std::acos(q / (p * zz * 2.0f)) / 3.0f;
   ```
   Due to float32 precision limits, $\frac{q}{2 p z}$ occasionally evaluated to $1.0000002$ or $-1.0000002$, triggering domain errors in `std::acos` and returning `NaN`.
   **Fix**: Clamped the ratio to $[-1.0f, 1.0f]$ via `std::clamp`.

3. **Topology and Point Count Handling**:
   Quadratic Bézier paths are indexed as triplets $(P_{2k}, P_{2k+1}, P_{2k+2})$. Paths with an even number of points previously appended a synthetic point that was frequently colinear with the last segment, triggering the division by zero. Proper segment pairing was enforced.

### BVH2D Acceleration for Bézier Curves

Since quadratic Bézier curves satisfy the convex hull property:
$$B(t) \in \text{ConvexHull}(P_0, P_1, P_2) \quad \forall t \in [0, 1]$$
A tight and conservative bounding box for each quadratic curve is simply:
$$\text{AABB} = [\min(P_0, P_1, P_2), \max(P_0, P_1, P_2)]$$

A specialized `BVH2DBezier` tree was implemented using this bounding box representation. Per-pixel queries prune subtrees whose bounding box distance exceeds the current closest squared distance.

The API was structured consistently with linear polyline SDF:
- `hmap::sdf_2d_polyline_bezier_bvh`
- `hmap::sdf_2d_polyline_bezier_brute_force`
- `hmap::sdf_2d_polyline_bezier(path, shape, bbox, p_noise_x, p_noise_y, use_bvh = true)`
