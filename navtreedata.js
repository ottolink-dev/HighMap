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
    [ "Core Concepts", "md_docs_2content_2getting__started_2core__concepts.html", [
      [ "1. Elevation Representations", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md110", [
        [ "Array", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md111", null ]
      ] ],
      [ "2. Bivariate Functions", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md112", null ],
      [ "3. Geometric Primitives", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md113", null ],
      [ "4. Talus vs. Slope", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md114", [
        [ "Talus (Grid-Space Elevation Difference)", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md115", null ],
        [ "Slope (Metric Surface Gradient)", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md116", null ],
        [ "Relationship: talus = slope / nx", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md117", null ]
      ] ],
      [ "5. Scattering", "md_docs_2content_2getting__started_2core__concepts.html#autotoc_md118", null ]
    ] ],
    [ "Quickstart", "md_docs_2content_2getting__started_2quickstart.html", [
      [ "1. Example", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md120", null ],
      [ "2. CMake Integration", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md121", null ],
      [ "3. Basic Operations", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md122", null ],
      [ "4. Next Steps", "md_docs_2content_2getting__started_2quickstart.html#autotoc_md123", null ]
    ] ],
    [ "HighMap Documentation", "md_docs_2index.html", [
      [ "Documentation", "md_docs_2index.html#autotoc_md125", [
        [ "Getting Started", "md_docs_2index.html#autotoc_md126", null ],
        [ "Developer Documentation", "md_docs_2index.html#autotoc_md127", null ]
      ] ],
      [ "Installation", "md_docs_2index.html#autotoc_md128", null ]
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
"classhmap_1_1Array.html#a000acd474d036adfa4e5e334139788ce",
"classhmap_1_1BigTiffWriter.html#ab76555869ed2459371d5187374ebb557",
"classhmap_1_1CoordFrame.html#ab1c8fa8ed0244f4a54759253b629255c",
"classhmap_1_1FieldFunction.html#abd0ed9592d4ff58306a41b1d118a9ae7",
"classhmap_1_1KDTree.html#a6a69c840c8d2421d96fa453063b135e2",
"classhmap_1_1PdfSampler.html",
"classhmap_1_1ScatterField.html#ac41168040f8d2d54d29e27b8fd9c3b3d",
"classhmap_1_1TileStorage.html",
"coherent__noise_8hpp.html#a041f105aa505b9f94c9d526e5d0978cc",
"dendry__array__control__function_8hpp_source.html",
"export__points__to__ply_8cpp.html#ae01b5471d28a8f4a9104e2ffea4484a9",
"filters_8hpp.html#ae475b4f5defd09d3e566f44dd2c21575",
"functions_vars_z.html",
"hydraulic__mise_8cpp.html#a7a5dcbf1ee87dd3c79068bfff522a6ed",
"interpolate__heightmap__virtual__array_8cpp.html#ace2076a7808458f99aeaeec2bc404d6b",
"math_2array_8hpp.html#aad00066161ccf5ba4d6f93ab48c10876",
"morphology_8hpp.html#af80bcff234bf1482cf26896b9d21d270",
"namespacehmap.html#a2e0228ee67d8c47cd32a3c9d9c6cb9d4",
"namespacehmap.html#a6108df0909730083cb6dd7cf4a217c1b",
"namespacehmap.html#aa4468d1e57ee67ba0be8fb1ab931045a",
"namespacehmap.html#ad1eb8e81ef31a0e0bc7a5bd093c3aefc",
"namespacehmap_1_1gpu.html#a26ed284f1df761074c2a49a2fbbd582e",
"namespacehmap_1_1log.html#a2b618829e9c711302f726a3ec3f4ff49",
"path__functions_8cpp.html#aca95a27cb4796288e6ac2e7bd265c611",
"profiles_8hpp.html#a25a4e4775f3c0918ed75f0503f5ed0a6",
"sdf__2d__polyline_8cpp.html#aee7ff7068cdfd831c2ad4891a38a817c",
"structhmap_1_1AABB2D.html#a9cdb78ad36d25875cb0ba86e20efd765",
"structhmap_1_1MiseParams.html#a3cb61732527fd8119798ef8370e96667",
"structhmap_1_1SpeciesLookup.html",
"structhmap_1_1VirtualArray.html#ae77e0c68919f9a4aa0ac291300d1d0f9",
"test__cloud_8cpp.html#a2635d53220155d41c0cdf04c61a5555f",
"test__gradient_8cpp.html#a20ea8c3252331571d45ccda6f3308653",
"test__point_8cpp.html#a61268b7c41e61f8da119dfa430717bde",
"test__validation_8cpp.html#a742dd3e5278dc15eeae1f659bbca1fb9",
"valley__fill_8cpp.html"
];

var SYNCONMSG = 'click to disable panel synchronisation';
var SYNCOFFMSG = 'click to enable panel synchronisation';