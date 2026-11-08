# Architecture Overview

This document describes the architectural layout of the HighMap codebase, its module hierarchy, and core design principles.

## 1. Module Hierarchy

The codebase is organized into layered modules with unidirectional dependencies. Higher layers depend on foundational storage and mathematical abstractions, while core containers remain independent of domain-specific algorithms.

```text
+---------------------------------------------------------------+
|                      Asset & Scene Export                     |
|  (export/, image, scene, mesh generators, 3D file formats)    |
+---------------------------------------------------------------+
|                      Domain Applications                      |
|  (erosion/, hydrology/, flora/, rocks/, carving, roads)       |
+---------------------------------------------------------------+
|                      Array Processing                         |
|  (filters, transform, morphology, blending, synthesis)        |
+---------------------------------------------------------------+
|                      Core Data Structures                     |
|  (Array, HeightMap, VirtualArray, Point, Cloud, Path, Graph)  |
+---------------------------------------------------------------+
|                      Foundational Math                        |
|  (algebra, functions, random, geometry, coordinate frames)    |
+---------------------------------------------------------------+
```

## 2. Layers

### Foundational Math
Provides mathematical primitives and analytical functions:
- Vector and matrix types (using GLM).
- Coordinate systems and transformation frames (`coord_frame.hpp`).
- Bivariate continuous functions (`functions.hpp`, `algebra.hpp`).
- Deterministic random sampling and noise kernels (`random.hpp`, `kernels.hpp`).

### Core Data Structures
Defines data containers without coupling to specific simulation or processing logic:
- `Array`: 2D contiguous row-major floating-point buffer.
- `HeightMap`: Tiled multi-buffer representation for concurrent block computation.
- `VirtualArray`: Out-of-core tiled raster supporting memory-mapped disk storage and distributed tile processing with halo exchange.
- Geometric primitives: `Point`, `Edge`, `Cloud`, `Path`, `Graph`.
- Entity scattering container: `ScatterField` and `ScatterItem`.

### Array Processing
Implements functional and in-place transformations operating on raster data:
- `filters.hpp`: Spatial convolution, bilateral filtering, domain transformations.
- `morphology.hpp`: Dilation, erosion, opening, closing, distance transforms.
- `synthesis.hpp` & `blending.hpp`: Procedural surface construction, layer combination operators.
- `transform.hpp`: Resampling, rotations, cropping, coordinate remapping.

### Domain Simulations
Encapsulates specialized physical, biological, and geomorphological models:
- `erosion/`: Hydraulic erosion, sediment transport, and thermal talus diffusion.
- `hydrology/`: Watershed segmentation, depression filling, flow accumulation, and river graph routing.
- `flora/`: Botanical descriptors (`Species`), multi-species competition, canopy dynamics (`Forest`).
- `rocks/`: Geomorphological rock distributions, talus sorting, and Pareto power-law radius sampling (`RockField`).
- `carving.hpp` & `roads.hpp`: Terrain cutting, path-guided grading, and slope alignment.

### Export & I/O
Converts internal representations into exchange formats:
- Raster images: 8-bit, 16-bit, and 32-bit floating-point images via OpenCV, OpenEXR, and BigTIFF writers.
- 3D representations: Triangular meshes (`TerrainTriMesh`), heightfield point clouds, USD / glTF scenes.

## 3. Design Principles

### Contiguous Storage and Explicit Memory Control
Core buffers (`Array`) rely on dense, contiguous allocations in row-major order. Iteration order throughout the codebase prioritizes cache coherence ($y$ outer loop, $x$ inner loop).

### Separation of Storage and Algorithm
Algorithms accept references or views to data buffers rather than binding to specific container subclasses. For large terrains, algorithms execute through tiled abstractions (`HeightMap` or `VirtualArray`) with explicit boundary halo policies.

### Stateless Processing Functions
Image processing, filtering, and synthesis algorithms are implemented primarily as free functions operating on inputs and producing outputs or modifying referenced buffers in-place, avoiding hidden state.

### Deterministic Seeding
All stochastic procedures (noise generators, scattering passes, Monte Carlo simulations) accept explicit integer seeds to ensure reproducibility across runs.

## 4. Subsystem Deep Dives

Detailed architectural descriptions of specific subsystems:

- [Scatter Architecture](scatter_field.md): Detailed specification of `ScatterField`, the unified instance model, and specializations `Forest` and `RockField`.
- [Tiled Communication & Overlap Synchronization](tiled_communication.md): Collective communication interface (`ICommunicator`), sequential cooperative scheduler, and halo buffer synchronization.

## 5. Investigations and Benchmarks

Explorations, trial-and-error evaluations, and algorithmic benchmarks:

- [FFT vs. GPU Convolution & Morphology](investigations/fft_convolution_and_morphology.md): Experimental evaluation of Fourier-domain convolution and soft-max morphology against spatial GPU implementations.
