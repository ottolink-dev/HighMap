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
        [ "Export &amp; I/O", "md_docs_2content_2dev_2architecture.html#autotoc_md8", null ]
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
        [ "Multi-Scale Cascade (<span class=\"tt\">hydraulic_particle_multiscale</span>)", "md_docs_2content_2dev_2hydraulic__particle.html#autotoc_md46", null ]
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
      [ "4. Multiscale Architecture (<span class=\"tt\">hydraulic_stream_log_multiscale</span>)", "md_docs_2content_2dev_2hydraulic__stream__log.html#autotoc_md57", [
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
      [ "5. Verification of Bézier Polyline SDF (<span class=\"tt\">sdf_2d_polyline_bezier</span>)", "md_docs_2content_2dev_2investigations_2bvh2d__path__sdf.html#autotoc_md70", [
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
        [ "Flora Subsystem (<span class=\"tt\">Forest</span>)", "md_docs_2content_2dev_2scatter__field.html#autotoc_md91", null ],
        [ "Rocks Subsystem (<span class=\"tt\">RockField</span>)", "md_docs_2content_2dev_2scatter__field.html#autotoc_md92", null ]
      ] ],
      [ "4. Converting to Raster Fields", "md_docs_2content_2dev_2scatter__field.html#autotoc_md93", null ]
    ] ],
    [ "Thermal Weathering and Talus Erosion Algorithms", "md_docs_2content_2dev_2thermal__erosion.html", [
      [ "1. Physical Foundations", "md_docs_2content_2dev_2thermal__erosion.html#autotoc_md95", null ]
    ] ],
    [ "Tiled Communication and Overlap Synchronization", "md_docs_2content_2dev_2tiled__communication.html", [
      [ "1. Overview and Purpose", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md110", null ],
      [ "2. Communicator Interface (<span class=\"tt\">ICommunicator</span>)", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md111", [
        [ "Supported Operations", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md112", null ],
        [ "Overlap Synchronization Modes (<span class=\"tt\">SyncOperation</span>)", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md113", null ]
      ] ],
      [ "3. Execution Engines", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md114", [
        [ "Single-Array Mode (<span class=\"tt\">SingleArrayComm</span>)", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md115", null ],
        [ "Multithreaded Distributed Mode (<span class=\"tt\">ThreadedTileComm</span>)", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md116", null ],
        [ "Out-of-Core Sequential Mode (<span class=\"tt\">SequentialTileComm</span>)", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md117", null ]
      ] ],
      [ "4. How Overlap Buffering Works", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md118", null ],
      [ "5. Iteration Batching in Domain Algorithms", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md119", [
        [ "Batching Loop Pattern", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md120", null ]
      ] ],
      [ "6. Integration and Usage Example", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md121", null ],
      [ "7. Buffer Shape Invariance and Nested Multiscale Pipelines", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md122", [
        [ "Incompatibility with Nested Multiscale Pyramids", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md123", null ],
        [ "Runtime Shape Validation", "md_docs_2content_2dev_2tiled__communication.html#autotoc_md124", null ]
      ] ]
    ] ],
    [ "Core Concepts", "md_docs_2content_2getting__started_2core__concepts.html", [
      [ "1. Elevation Representations", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md126", [
        [ "Array", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md127", null ]
      ] ],
      [ "2. Bivariate Functions", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md128", null ],
      [ "3. Geometric Primitives", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md129", null ],
      [ "4. Talus vs. Slope", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md130", [
        [ "Talus (Grid-Space Elevation Difference)", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md131", null ],
        [ "Slope (Metric Surface Gradient)", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md132", null ],
        [ "Relationship: talus = slope / nx", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md133", null ]
      ] ],
      [ "5. Scattering", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md134", null ],
      [ "6. Spatial Coordinates: UV, BBox, and CoordFrame", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md135", [
        [ "UV (<span class=\"tt\">[0, 1]^2</span> Normalized Domain)", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md136", null ],
        [ "BBox (<span class=\"tt\">BBox2D</span> Axis-Aligned Bounding Box)", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md137", null ],
        [ "CoordFrame (<span class=\"tt\">CoordFrame</span> Oriented Local Frame)", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md138", null ]
      ] ]
    ] ],
    [ "Quickstart", "md_docs_2content_2getting__started_2quickstart.html", [
      [ "1. Example", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md140", null ],
      [ "2. CMake Integration", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md141", null ],
      [ "3. Basic Operations", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md142", null ],
      [ "4. Next Steps", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md143", null ]
    ] ],
    [ "HighMap Documentation", "md_docs_2index.html", [
      [ "Documentation", "md_docs_2index.html#autotoc_md145", [
        [ "Getting Started", "md_docs_2index.html#autotoc_md146", null ],
        [ "Developer Documentation", "md_docs_2index.html#autotoc_md147", null ]
      ] ],
      [ "Installation", "md_docs_2index.html#autotoc_md148", null ]
    ] ],
    [ "Optional Tracy CPU profiling", "md_docs_2profiling.html", null ],
    [ "Deprecated List", "deprecated.html", null ],
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
        [ "Functions", "globals_func.html", "globals_func" ]
      ] ]
    ] ],
    [ "Examples", "examples.html", "examples" ]
  ] ]
];

var NAVTREEINDEX =
[
"add__noise_8cpp.html",
"classhmap_1_1Array.html#af586ddce153ac93df48a6d26609d7264",
"classhmap_1_1Cloud.html#aafa7e8dda6f65d37f30c3fa8414481c0",
"classhmap_1_1FieldFunction.html#a17613f0a81c1605e901170bf5224afbc",
"classhmap_1_1InverseSampler2D.html#ae33248c92a5bd80e35291f1e1cea16f8",
"classhmap_1_1PdfSampler.html",
"classhmap_1_1ScatterField.html#aa5cc94dd2526613f4d8f1a7db128e3c0",
"classhmap_1_1TerrainTriMesh.html#adab683cb7235023ccd24630f401f6755",
"classhmap_1_1WaveDuneFunction.html#ae329aa6a15d2bc57ee3a9dd54c9f1387",
"distance__transform_8cpp.html#a345246095a5b2601fa9cd2dacf851cbc",
"interpolate2d_8hpp.html",
"namespacehmap.html#a11c1600c2d67387c1b037fadcca9a069",
"namespacehmap.html#a3d23e9f95cfe4a8ceeecd6ea01d1d3e4",
"namespacehmap.html#a6b058e4341f421d058c4d0d408d30654",
"namespacehmap.html#ab0307c51c15ac3d4162974c970bc9ea9",
"namespacehmap.html#ad8482d5a5c0650c20ffe1dc86b5f8a61",
"namespacehmap_1_1gpu.html#a4ce73a2eef313a245bbe81f115695a6b",
"namespacehmap_1_1neighborhood.html#a6b1a4665bde6beaf3b96d40d027083a7",
"structhmap_1_1AABB2D.html#a02c15d457826615942f7bc826be9a88b",
"structhmap_1_1GaussianPush.html#a46e45b11a25e6e431318577e44f1a2e4",
"structhmap_1_1Recorder.html#a0ec174294138be0e1d677a7c33011959",
"structhmap_1_1TerrainTriMesh_1_1Neighbor.html#a1d8844b86573fc6e3e0037a6f9d4632c",
"structhmap_1_1gpu_1_1McDonaldParams_1_1PhysicalParams.html#a0bb1f5db548cbbb733e1d49e2c958b89",
"test__convolution__scaling_8cpp.html#a072356e877b58a1c0a9d10abde9680b4",
"test__grids_8cpp.html#a5e0ea3f9036f0624cc417d15b4b19013",
"test__point_8cpp.html#ae071b673bea366a8b822ebf6b1748b22",
"test__tiled__monolithic_8cpp.html#a252b2def3dc0b1630b406be8c3974ea2"
];

var SYNCONMSG = 'click to disable panel synchronization';
var SYNCOFFMSG = 'click to enable panel synchronization';
var LISTOFALLMEMBERS = 'List of all members';