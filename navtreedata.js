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
    [ "Hydraulic Erosion Models Comparison", "md_docs_2content_2dev_2hydraulic__comparison.html", [
      [ "1. Visual Comparison", "md_docs_2content_2dev_2hydraulic__comparison.html#autotoc_md17", [
        [ "Erosion Difference Maps", "md_docs_2content_2dev_2hydraulic__comparison.html#autotoc_md18", null ]
      ] ],
      [ "2. Qualitative Characteristics", "md_docs_2content_2dev_2hydraulic__comparison.html#autotoc_md19", null ],
      [ "3. Algorithm Deep Dives", "md_docs_2content_2dev_2hydraulic__comparison.html#autotoc_md20", null ]
    ] ],
    [ "Multiscale Implicit Stream-Power Erosion (MISE)", "md_docs_2content_2dev_2hydraulic__mise.html", [
      [ "1. Overview", "md_docs_2content_2dev_2hydraulic__mise.html#autotoc_md22", null ],
      [ "2. Theoretical Formulation and Kinematics", "md_docs_2content_2dev_2hydraulic__mise.html#autotoc_md23", [
        [ "Drainage Graph and Depression Filling", "md_docs_2content_2dev_2hydraulic__mise.html#autotoc_md24", null ]
      ] ]
    ] ],
    [ "Hydraulic Particle Erosion", "md_docs_2content_2dev_2hydraulic__particle.html", [
      [ "1. Overview", "md_docs_2content_2dev_2hydraulic__particle.html#autotoc_md35", null ],
      [ "2. Particle Kinematics and Trajectories", "md_docs_2content_2dev_2hydraulic__particle.html#autotoc_md36", [
        [ "Surface Interpolation and Gradient Evaluation", "md_docs_2content_2dev_2hydraulic__particle.html#autotoc_md37", null ],
        [ "Trajectory Integration and Momentum Blending", "md_docs_2content_2dev_2hydraulic__particle.html#autotoc_md38", null ],
        [ "Velocity and Energy Conservation", "md_docs_2content_2dev_2hydraulic__particle.html#autotoc_md39", null ]
      ] ],
      [ "3. Sediment Transport Dynamics", "md_docs_2content_2dev_2hydraulic__particle.html#autotoc_md40", [
        [ "Deposition", "md_docs_2content_2dev_2hydraulic__particle.html#autotoc_md41", null ],
        [ "Erosion", "md_docs_2content_2dev_2hydraulic__particle.html#autotoc_md42", null ]
      ] ],
      [ "4. Bank Collapse and Thermal Relaxation", "md_docs_2content_2dev_2hydraulic__particle.html#autotoc_md43", null ],
      [ "5. Multi-Pass and Multi-Scale Architectures", "md_docs_2content_2dev_2hydraulic__particle.html#autotoc_md44", [
        [ "Multi-Pass Iterations", "md_docs_2content_2dev_2hydraulic__particle.html#autotoc_md45", null ],
        [ "Multi-Scale Cascade (<tt>hydraulic_particle_multiscale</tt>)", "md_docs_2content_2dev_2hydraulic__particle.html#autotoc_md46", null ]
      ] ],
      [ "6. C++ API Reference", "md_docs_2content_2dev_2hydraulic__particle.html#autotoc_md47", null ]
    ] ],
    [ "Hydraulic Stream Log Erosion", "md_docs_2content_2dev_2hydraulic__stream__log.html", [
      [ "1. Overview", "md_docs_2content_2dev_2hydraulic__stream__log.html#autotoc_md49", null ],
      [ "2. Mathematical Model and Kinematics", "md_docs_2content_2dev_2hydraulic__stream__log.html#autotoc_md50", [
        [ "Continuous Flow Accumulation ($D_\\infty$)", "md_docs_2content_2dev_2hydraulic__stream__log.html#autotoc_md51", null ],
        [ "Logarithmic Compression", "md_docs_2content_2dev_2hydraulic__stream__log.html#autotoc_md52", null ],
        [ "Gradient Modulation and Stream Power", "md_docs_2content_2dev_2hydraulic__stream__log.html#autotoc_md53", null ],
        [ "Ridge and Peak Preservation", "md_docs_2content_2dev_2hydraulic__stream__log.html#autotoc_md54", null ],
        [ "Surface Incision", "md_docs_2content_2dev_2hydraulic__stream__log.html#autotoc_md55", null ]
      ] ],
      [ "3. Deposition Mechanics", "md_docs_2content_2dev_2hydraulic__stream__log.html#autotoc_md56", null ],
      [ "4. Multiscale Architecture (<tt>hydraulic_stream_log_multiscale</tt>)", "md_docs_2content_2dev_2hydraulic__stream__log.html#autotoc_md57", [
        [ "Multiscale Algorithm Steps", "md_docs_2content_2dev_2hydraulic__stream__log.html#autotoc_md58", null ]
      ] ],
      [ "5. C++ API Reference", "md_docs_2content_2dev_2hydraulic__stream__log.html#autotoc_md59", null ]
    ] ],
    [ "BVH2D Spatial Indexing for Path SDF Computation", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html", [
      [ "1. Context and Motivation", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md61", null ],
      [ "2. Algorithmic Approach", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md62", [
        [ "Segment Indexing", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md63", null ],
        [ "Branch-and-Bound Distance Query", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md64", null ],
        [ "Support for Domain Warping and Path Closure", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md65", null ]
      ] ],
      [ "3. Benchmark Measurements", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md66", [
        [ "Execution Time Comparison ($512 \\times 512$ Grid)", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md67", null ],
        [ "Execution Time Comparison ($256 \\times 256$ Grid)", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md68", null ]
      ] ],
      [ "4. API Design and Integration", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md69", null ],
      [ "5. Verification of Bézier Polyline SDF (<tt>sdf_2d_polyline_bezier</tt>)", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md70", [
        [ "Root Cause Analysis", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md71", null ],
        [ "BVH2D Acceleration for Bézier Curves", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md72", null ]
      ] ],
      [ "6. GPU (OpenCL) Assessment and Validation", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md73", [
        [ "Performance Analysis: GPU Brute Force vs CPU BVH2D", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md74", null ],
        [ "Rationale: Keeping the GPU Kernel Brute Force", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md75", null ],
        [ "Numerical Stability and Parity Fixes on GPU", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md76", null ]
      ] ]
    ] ],
    [ "FFT vs. GPU Convolution and Morphological Expansion", "md_docs_2content_2dev_2investigations_2fft__convolution__and__morphology.html", [
      [ "1. Context and Motivation", "md_docs_2content_2dev_2investigations_2fft__convolution__and__morphology.html#autotoc_md78", null ]
    ] ],
    [ "Scatter Architecture: ScatterField, Forest, and RockField", "md_docs_2content_2dev_2scatter__field.html", [
      [ "1. Class and Dependency Diagram", "md_docs_2content_2dev_2scatter__field.html#autotoc_md85", null ],
      [ "2. Design Principles", "md_docs_2content_2dev_2scatter__field.html#autotoc_md86", [
        [ "Uniform Instance Representation", "md_docs_2content_2dev_2scatter__field.html#autotoc_md87", null ],
        [ "Inheritance and Reuse", "md_docs_2content_2dev_2scatter__field.html#autotoc_md88", null ],
        [ "Descriptor Separation", "md_docs_2content_2dev_2scatter__field.html#autotoc_md89", null ]
      ] ],
      [ "3. Subsystem Connections", "md_docs_2content_2dev_2scatter__field.html#autotoc_md90", [
        [ "Flora Subsystem (<tt>Forest</tt>)", "md_docs_2content_2dev_2scatter__field.html#autotoc_md91", null ],
        [ "Rocks Subsystem (<tt>RockField</tt>)", "md_docs_2content_2dev_2scatter__field.html#autotoc_md92", null ]
      ] ],
      [ "4. Converting to Raster Fields", "md_docs_2content_2dev_2scatter__field.html#autotoc_md93", null ]
    ] ],
    [ "Thermal Weathering and Talus Erosion Algorithms", "md_docs_2content_2dev_2thermal__erosion.html", [
      [ "1. Physical Foundations", "md_docs_2content_2dev_2thermal__erosion.html#autotoc_md95", null ]
    ] ],
    [ "Tiled Communication and Overlap Synchronization", "md_docs_2content_2dev_2tiled__communication.html", [
      [ "1. Overview and Purpose", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md110", null ],
      [ "2. Communicator Interface (<tt>ICommunicator</tt>)", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md111", [
        [ "Supported Operations", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md112", null ],
        [ "Overlap Synchronization Modes (<tt>SyncOperation</tt>)", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md113", null ]
      ] ],
      [ "3. Execution Engines", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md114", [
        [ "Single-Array Mode (<tt>SingleArrayComm</tt>)", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md115", null ],
        [ "Multithreaded Distributed Mode (<tt>ThreadedTileComm</tt>)", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md116", null ],
        [ "Out-of-Core Sequential Mode (<tt>SequentialTileComm</tt>)", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md117", null ]
      ] ],
      [ "4. How Overlap Buffering Works", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md118", null ],
      [ "5. Iteration Batching in Domain Algorithms", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md119", [
        [ "Batching Loop Pattern", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md120", null ]
      ] ],
      [ "6. Integration and Usage Example", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md121", null ]
    ] ],
    [ "Core Concepts", "md_docs_2content_2getting__started_2core__concepts.html", [
      [ "1. Elevation Representations", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md123", [
        [ "Array", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md124", null ]
      ] ],
      [ "2. Bivariate Functions", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md125", null ],
      [ "3. Geometric Primitives", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md126", null ],
      [ "4. Talus vs. Slope", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md127", [
        [ "Talus (Grid-Space Elevation Difference)", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md128", null ],
        [ "Slope (Metric Surface Gradient)", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md129", null ],
        [ "Relationship: talus = slope / nx", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md130", null ]
      ] ],
      [ "5. Scattering", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md131", null ]
    ] ],
    [ "Quickstart", "md_docs_2content_2getting__started_2quickstart.html", [
      [ "1. Example", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md133", null ],
      [ "2. CMake Integration", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md134", null ],
      [ "3. Basic Operations", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md135", null ],
      [ "4. Next Steps", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md136", null ]
    ] ],
    [ "HighMap Documentation", "md_docs_2index.html", [
      [ "Documentation", "md_docs_2index.html#autotoc_md138", [
        [ "Getting Started", "md_docs_2index.html#autotoc_md139", null ],
        [ "Developer Documentation", "md_docs_2index.html#autotoc_md140", null ]
      ] ],
      [ "Installation", "md_docs_2index.html#autotoc_md141", null ]
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
"classJaggedTest.html#a1f456dd76286876c4f22c4255a6074b6",
"classhmap_1_1BigTiffWriter.html#a82f4cd33c936daaf0da9c419971c7152",
"classhmap_1_1CoordFrame.html#a8d816754189a8244628bd4c6eefaf09c",
"classhmap_1_1Forest.html",
"classhmap_1_1KDTree.html#a28e014b1d02cf4a168ec74be0cb28ab1",
"classhmap_1_1PerlinHalfFunction.html#a1009eb75ebcbd52423fd37eb8a2de011",
"classhmap_1_1ScatterField.html#aeb89e13a2a22fc06a08bc02e0fdcaa61",
"classhmap_1_1Texture.html#a7e68630c4908e71e9cfb91e3e891d472",
"classhmap_1_1WorleyDoubleFunction.html#af26f43647e20fe6fc3ff88ed4c1e1274",
"core_8cpp.html#a0a3ae79b31c4eebc54ae686d5ba558a3",
"distance__transform_8cpp.html#abcb9259cf91bfd25e92a72cdf18c285e",
"filters_8hpp.html#a63d764200510382a7732cfaa750df4a4",
"functions_8hpp.html#ab8197c4c7a85dacb616fee7a7d0df195a803cca1f65b9596dbde7c54039773a39",
"grow__forest__soft__core_8cpp.html#a12b36f2c693eff09dda2732199468fe8",
"interpolate2d_8hpp.html#a304d2b2e755c4a1fea14ac8dc6d44082",
"local__metrics__gpu_8cpp.html#a36a5f96566c4074dd1032c63bf95ee5e",
"md_docs_2profiling.html",
"namespacehmap.html#a1ba9d4ce7afe3309f8ae8e617f7154eeadb9a04ab3508704543b5dd43ae495fe7",
"namespacehmap.html#a4c42905cbb910cb4017ef3e026a31d59",
"namespacehmap.html#a8a7a083e583afe812239d1fe6aa479cd",
"namespacehmap.html#ac14815cf025b7430d31d08125a5c6f0d",
"namespacehmap.html#aef71b544e3f4e5321cb3d41a7366e0f3",
"namespacehmap_1_1gpu.html#ab0737b29664736928a86f0755493895fa3fb8064f9087e5af79c85dce2c655112",
"opencv__wrapper_8cpp.html#a518f05118e1b9704160a7c3276a2f483",
"primitives_2functions_8hpp.html#a02c2283ef2c4fedf45dcbed47a86f853",
"reverse__midpoint_8cpp.html#af5440364a346b297246885e95793ef45",
"solve__poisson__gauss__seidel_8cpp.html#a0bcc84aa462aee73e1a3a52729f21ea6",
"structhmap_1_1ImageWriterConfig.html",
"structhmap_1_1RockSeedingOptions.html#abd5f9bec8368c3669a6803a84a8d3720",
"structhmap_1_1TileKeyHash.html#a7231b3026e5cbd723cba9e9ee184f500",
"terrain__features_8hpp_source.html",
"test__erosion_8cpp.html#a8b502ba0361a7963fa843302373b2e6f",
"test__math_8cpp.html",
"test__rocks_8cpp.html#affe75027effe93d9104900729a70eb3c",
"thermal__erosion_8hpp.html#a72c9d220589d2724b4e48e9c37ff55ad",
"virtual__texture_8inl.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';