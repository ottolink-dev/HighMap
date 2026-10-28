/*
 @licstart  The following is the entire license notice for the JavaScript code in this file.

 The MIT License (MIT)

 Copyright (C) 1997-2020 by Dimitri van Heesch

 Permission is hereby granted, free of charge, to any person obtaining a copy of this software
 and associated documentation files (the "Software"), to deal in the Software without restriction,
 including without limitation the rights to use, copy, modify, merge, publish, distribute,
 sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all copies or
 substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
 BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

 @licend  The above is the entire license notice for the JavaScript code in this file
*/
var NAVTREE =
[
  [ "HighMap library (C++)", "index.html", [
    [ "Architecture Overview", "md_docs_2content_2dev_2architecture.html", [
      [ "1. Module Hierarchy", "md_docs_2content_2dev_2architecture.html#autotoc_md2", null ],
      [ "2. Layers", "md_docs_2content_2dev_2architecture.html#autotoc_md3", [
        [ "Foundational Math", "md_docs_2content_2dev_2architecture.html#autotoc_md4", null ],
        [ "Core Data Structures", "md_docs_2content_2dev_2architecture.html#autotoc_md5", null ],
        [ "Array Processing", "md_docs_2content_2dev_2architecture.html#autotoc_md6", null ],
        [ "Domain Simulations", "md_docs_2content_2dev_2architecture.html#autotoc_md7", null ],
        [ "Export & I/O", "md_docs_2content_2dev_2architecture.html#autotoc_md8", null ]
      ] ],
      [ "3. Design Principles", "md_docs_2content_2dev_2architecture.html#autotoc_md9", [
        [ "Contiguous Storage and Explicit Memory Control", "md_docs_2content_2dev_2architecture.html#autotoc_md10", null ],
        [ "Separation of Storage and Algorithm", "md_docs_2content_2dev_2architecture.html#autotoc_md11", null ],
        [ "Stateless Processing Functions", "md_docs_2content_2dev_2architecture.html#autotoc_md12", null ],
        [ "Deterministic Seeding", "md_docs_2content_2dev_2architecture.html#autotoc_md13", null ]
      ] ],
      [ "4. Subsystem Deep Dives", "md_docs_2content_2dev_2architecture.html#autotoc_md14", null ],
      [ "5. Investigations and Benchmarks", "md_docs_2content_2dev_2architecture.html#autotoc_md15", null ]
    ] ],
    [ "BVH2D Spatial Indexing for Path SDF Computation", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html", [
      [ "1. Context and Motivation", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md17", null ],
      [ "2. Algorithmic Approach", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md18", [
        [ "Segment Indexing", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md19", null ],
        [ "Branch-and-Bound Distance Query", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md20", null ],
        [ "Support for Domain Warping and Path Closure", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md21", null ]
      ] ],
      [ "3. Benchmark Measurements", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md22", [
        [ "Execution Time Comparison ($512 \\times 512$ Grid)", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md23", null ],
        [ "Execution Time Comparison ($256 \\times 256$ Grid)", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md24", null ]
      ] ],
      [ "4. API Design and Integration", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md25", null ],
      [ "5. Verification of Bézier Polyline SDF (<tt>sdf_2d_polyline_bezier</tt>)", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md26", [
        [ "Root Cause Analysis", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md27", null ],
        [ "BVH2D Acceleration for Bézier Curves", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md28", null ]
      ] ],
      [ "6. GPU (OpenCL) Assessment and Validation", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md29", [
        [ "Performance Analysis: GPU Brute Force vs CPU BVH2D", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md30", null ],
        [ "Rationale: Keeping the GPU Kernel Brute Force", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md31", null ],
        [ "Numerical Stability and Parity Fixes on GPU", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md32", null ]
      ] ]
    ] ],
    [ "FFT vs. GPU Convolution and Morphological Expansion", "md_docs_2content_2dev_2investigations_2fft__convolution__and__morphology.html", [
      [ "1. Context and Motivation", "md_docs_2content_2dev_2investigations_2fft__convolution__and__morphology.html#autotoc_md34", null ]
    ] ],
    [ "Scatter Architecture: ScatterField, Forest, and RockField", "md_docs_2content_2dev_2scatter__field.html", [
      [ "1. Class and Dependency Diagram", "md_docs_2content_2dev_2scatter__field.html#autotoc_md41", null ],
      [ "2. Design Principles", "md_docs_2content_2dev_2scatter__field.html#autotoc_md42", [
        [ "Uniform Instance Representation", "md_docs_2content_2dev_2scatter__field.html#autotoc_md43", null ],
        [ "Inheritance and Reuse", "md_docs_2content_2dev_2scatter__field.html#autotoc_md44", null ],
        [ "Descriptor Separation", "md_docs_2content_2dev_2scatter__field.html#autotoc_md45", null ]
      ] ],
      [ "3. Subsystem Connections", "md_docs_2content_2dev_2scatter__field.html#autotoc_md46", [
        [ "Flora Subsystem (<tt>Forest</tt>)", "md_docs_2content_2dev_2scatter__field.html#autotoc_md47", null ],
        [ "Rocks Subsystem (<tt>RockField</tt>)", "md_docs_2content_2dev_2scatter__field.html#autotoc_md48", null ]
      ] ],
      [ "4. Converting to Raster Fields", "md_docs_2content_2dev_2scatter__field.html#autotoc_md49", null ]
    ] ],
    [ "Core Concepts", "md_docs_2content_2getting__started_2core__concepts.html", [
      [ "1. Elevation Representations", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md51", [
        [ "Array", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md52", null ]
      ] ],
      [ "2. Bivariate Functions", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md53", null ],
      [ "3. Geometric Primitives", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md54", null ],
      [ "4. Scattering", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md55", null ]
    ] ],
    [ "Quickstart", "md_docs_2content_2getting__started_2quickstart.html", [
      [ "1. Example", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md57", null ],
      [ "2. CMake Integration", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md58", null ],
      [ "3. Basic Operations", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md59", null ],
      [ "4. Next Steps", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md60", null ]
    ] ],
    [ "HighMap Documentation", "md_docs_2index.html", [
      [ "Documentation", "md_docs_2index.html#autotoc_md62", [
        [ "Getting Started", "md_docs_2index.html#autotoc_md63", null ],
        [ "Developer Documentation", "md_docs_2index.html#autotoc_md64", null ]
      ] ],
      [ "Installation", "md_docs_2index.html#autotoc_md65", null ]
    ] ],
    [ "Optional Tracy CPU profiling", "md_docs_2profiling.html", null ],
    [ "Bibliography", "citelist.html", null ],
    [ "Namespaces", "namespaces.html", [
      [ "Namespace List", "namespaces.html", "namespaces_dup" ],
      [ "Namespace Members", "namespacemembers.html", [
        [ "All", "namespacemembers.html", "namespacemembers_dup" ],
        [ "Functions", "namespacemembers_func.html", "namespacemembers_func" ],
        [ "Variables", "namespacemembers_vars.html", null ],
        [ "Typedefs", "namespacemembers_type.html", null ],
        [ "Enumerations", "namespacemembers_enum.html", null ],
        [ "Enumerator", "namespacemembers_eval.html", "namespacemembers_eval" ]
      ] ]
    ] ],
    [ "Classes", "annotated.html", [
      [ "Class List", "annotated.html", "annotated_dup" ],
      [ "Class Index", "classes.html", null ],
      [ "Class Hierarchy", "hierarchy.html", "hierarchy" ],
      [ "Class Members", "functions.html", [
        [ "All", "functions.html", "functions_dup" ],
        [ "Functions", "functions_func.html", "functions_func" ],
        [ "Variables", "functions_vars.html", "functions_vars" ],
        [ "Typedefs", "functions_type.html", null ],
        [ "Enumerations", "functions_enum.html", null ],
        [ "Related Symbols", "functions_rela.html", null ]
      ] ]
    ] ],
    [ "Files", "files.html", [
      [ "File List", "files.html", "files_dup" ],
      [ "File Members", "globals.html", [
        [ "All", "globals.html", "globals_dup" ],
        [ "Functions", "globals_func.html", "globals_func" ],
        [ "Variables", "globals_vars.html", null ]
      ] ]
    ] ],
    [ "Examples", "examples.html", "examples" ]
  ] ]
];

var NAVTREEINDEX =
[
"add__noise_8cpp.html",
"classhmap_1_1Array.html#a0240a827bf1e8e52f4448746d244ea31",
"classhmap_1_1BigTiffWriter.html#ae4cc55e3fa4015a5bb35205d137cba4d",
"classhmap_1_1CoordFrame.html#aebce63606ddca495d41de463cafb33bd",
"classhmap_1_1FieldFunction.html#acfaece8d19fe130f277072fb82b3b673",
"classhmap_1_1KDTree.html#a957d4e8b3845faf98b017d6fa0eee2b5",
"classhmap_1_1PdfSampler.html#ae36e2260565931cbb6a9f6b1c8357fd0",
"classhmap_1_1ScatterField.html#ad4986df1aa3dff53929c88f0d5b26d33",
"classhmap_1_1TileStorage.html#a2feeb317dfb3fa8057fcd4daf4202fd8",
"coherent__noise_8hpp.html#a10ace74a842f9d2dc1f15dde91263a40",
"deposition_8cpp.html#a26ed284f1df761074c2a49a2fbbd582e",
"export__splatmap_8cpp.html",
"filters_8hpp.html#ae6c94ff5fe464c75e5ef05b05d5ccacf",
"functions_w.html",
"hydraulic__mise_8cpp.html#a9f4e20a91b26459133b83992e9d8d381",
"inverse__sampler__2d_8hpp_source.html",
"math_2array_8hpp.html#abc58609556b107fdb1825a6c423efdce",
"namespacehmap.html#a01df0f07ab6d0583d0794d7c03abc1eb",
"namespacehmap.html#a32464922c39713910ebeab632e0e8fd1",
"namespacehmap.html#a69877ef43645c767daa0735fde59821c",
"namespacehmap.html#aad6a39e70bf2ff70f8367e663d3027cf",
"namespacehmap.html#ad6fc7fb24960fa6a45bc657ae2d89f80",
"namespacehmap_1_1gpu.html#a5acf81d4e6a4f4950a9c5f0dd936e75d",
"namespacemembers.html",
"poi_8hpp.html#ade234d7bb6c23c67317077ab7a9082caa729ef3b0198a30fc445b80eed5bc745c",
"range_8cpp.html#a55574a9d91e2e7e20de954956402ba1b",
"selector_8hpp.html#a60843740c9f6153042e0466b696239fd",
"structhmap_1_1ColorAdjust.html#abf71052fb9aa59345fff4155d9b33c0d",
"structhmap_1_1OpenEXRWriter_1_1Impl.html#a8c7de7055e5343eb5ab2b52e7bc879e8",
"structhmap_1_1TerrainTriMesh_1_1BoundingBox.html#a5e5137be28d3192dc989d549fc957639",
"structhmap_1_1gpu_1_1McDonaldParams_1_1PhysicalParams.html#a6957f3c429e254be4d008b04377c9f8c",
"test__distance__transform_8cpp.html#a15c97e0ceaccca23fb324bdebbc881e2",
"test__kd__tree_8cpp.html",
"test__remap_8cpp.html#af3a227c7111847b66b5b17c0f244d0c2",
"texture__transform_8cpp.html#a7d27d0392fb6850a212bee00cb58700c",
"virtual__array_8hpp.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';