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
"classhmap_1_1Array.html#afb1e51a296eee1c7375afa500877b261",
"classhmap_1_1Cloud.html#aab450c9376e44def1d4e4465f0b1b691",
"classhmap_1_1DrainageBasinCellBased.html#ad576fdd4931b25ab57bec6eeca32204b",
"classhmap_1_1Graph.html#af9dd5d2b3ff699e1e064f4d189d74d55",
"classhmap_1_1Path.html#a0121bc7300b5a436eea60d3915103458",
"classhmap_1_1RockField.html#aa4616c39e3dd683da14cd2aff13bce80",
"classhmap_1_1StepFunction.html#aae15d0fef186bc37378c6dfc029dcbc7",
"classhmap_1_1ValueLinearNoiseFunction.html#aa0666aa17538459d7ca10513acac9203",
"dir_211501a338c3c7d58f3c9708c2a7240f.html",
"geometry_8hpp.html",
"morphology__gpu_8cpp.html",
"namespacehmap.html#a321149eda3417aa485a010c31267818ca0e9abbb30a51ad2951f6f0b527f6f7a1",
"namespacehmap.html#a62068967763814bfbf05a505993c544a",
"namespacehmap.html#aa3322e388ba870c9d6881082da2fd800",
"namespacehmap.html#ad160ef479ee398f252b4f12d02a47cb3",
"namespacehmap_1_1gpu.html#a13e5b18eb4ff1cf1f7357a9b29bf9368",
"namespacehmap_1_1gpu.html#aec782e6aba3f4cfbbd97aafa46cc700b",
"scatter__field_8hpp.html",
"structhmap_1_1ImageWriterConfig.html#aa40e252246ea357b6cc3c6f795651856",
"structhmap_1_1RockSimulationOptions.html#a0368f35abadef946b682e9b0762224c5",
"structhmap_1_1TileKeyHash.html#a7231b3026e5cbd723cba9e9ee184f500",
"terrain__tri__mesh_8hpp_source.html",
"test__erosion_8cpp.html#a96a4f2dbf67bc7e1100baf25ad164c77",
"test__logger_8cpp.html#acfaf7091403bf7bff33d37e9ed5fd5ef",
"test__rocks_8cpp.html#a38ea00803f849b75b5f278092a67f8f6",
"texture_8hpp.html"
];

var SYNCONMSG = 'click to disable panel synchronization';
var SYNCOFFMSG = 'click to enable panel synchronization';
var LISTOFALLMEMBERS = 'List of all members';