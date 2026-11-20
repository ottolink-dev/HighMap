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
      [ "5. Scattering", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md134", null ]
    ] ],
    [ "Quickstart", "md_docs_2content_2getting__started_2quickstart.html", [
      [ "1. Example", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md136", null ],
      [ "2. CMake Integration", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md137", null ],
      [ "3. Basic Operations", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md138", null ],
      [ "4. Next Steps", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md139", null ]
    ] ],
    [ "HighMap Documentation", "md_docs_2index.html", [
      [ "Documentation", "md_docs_2index.html#autotoc_md141", [
        [ "Getting Started", "md_docs_2index.html#autotoc_md142", null ],
        [ "Developer Documentation", "md_docs_2index.html#autotoc_md143", null ]
      ] ],
      [ "Installation", "md_docs_2index.html#autotoc_md144", null ]
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
        [ "Functions", "globals_func.html", "globals_func" ]
      ] ]
    ] ],
    [ "Examples", "examples.html", "examples" ]
  ] ]
];

var NAVTREEINDEX =
[
"add__noise_8cpp.html",
"classhmap_1_1Array.html#afdba5b214c8378e78db1978188676703",
"classhmap_1_1Cloud.html#aacfbabf4294d4922d79dc9a94e501592",
"classhmap_1_1DrainageBasinCellBased.html#aeb36d7eba63e150534906f2e924e7be8",
"classhmap_1_1Graph.html#afe8d4e58ca035b43ed891556dd7cdd2c",
"classhmap_1_1Path.html#a0ef4a309bc181ea182c0b7a75c8fdae7",
"classhmap_1_1RockField.html#ab05996c590e4e7037dd3583da7c9359a",
"classhmap_1_1StepFunction.html#ad75592284cc760c917c42578878d99b9",
"classhmap_1_1ValueLinearNoiseFunction.html#adcd66fd0064d69239703672a996c7665",
"dir_28d1f6f5761f5c8b206221dff91b5473.html",
"globals.html",
"mountain__tibesti_8cpp.html",
"namespacehmap.html#a321149eda3417aa485a010c31267818ca8083c22f48c0732c103514d10b5fe7e6",
"namespacehmap.html#a660d632fc99c8af533a9d3b2326d8e7d",
"namespacehmap.html#aa7a5c39d07988b3dc0aaba918618ed93",
"namespacehmap.html#ad39c72e2b4c07d3013f7b0143c31137c",
"namespacehmap_1_1gpu.html#a28e454a987e0ad5d10821c91cc11136c",
"namespacehmap_1_1log.html#a2b618829e9c711302f726a3ec3f4ff49",
"selector__gpu_8cpp.html",
"structhmap_1_1InteractionMatrix.html#a7c486745d36fc7e001b88d6d0f3ae1fe",
"structhmap_1_1ScatterSeedingOptions.html",
"structhmap_1_1TileRegion.html#a9c49add411fe1d63bdc173f330d09958",
"test__band_8cpp.html",
"test__flora_8cpp.html",
"test__math_8cpp.html#aee56917134f4d3eedc45e35a894c3b05",
"test__scatter_8cpp.html#a3d92cac9d3c7ca43dae7c3a4db8a6e56",
"virtual__array_8inl.html#a54eec3c94e3409f3990876e68a581322"
];

var SYNCONMSG = 'click to disable panel synchronization';
var SYNCOFFMSG = 'click to enable panel synchronization';
var LISTOFALLMEMBERS = 'List of all members';