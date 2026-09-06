cmake_minimum_required(VERSION 3.22)

# -----------------------------------------------------------------------------
# Optional: force clang Usage: cmake -B build -DCMAKE_C_COMPILER=clang
# -DCMAKE_CXX_COMPILER=clang++
# -----------------------------------------------------------------------------

project(highmap-root LANGUAGES C CXX)

# -----------------------------------------------------------------------------
# Options
# -----------------------------------------------------------------------------

option(HIGHMAP_ENABLE_OPENCL "" ON)
option(HIGHMAP_ENABLE_LOGS "Enable logging" ON)
option(HIGHMAP_ENABLE_EXAMPLES "" ON)
option(HIGHMAP_ENABLE_TESTS "" ON)
option(HIGHMAP_ENABLE_BENCHMARKS "" ON)
option(HIGHMAP_ENABLE_PROFILE_BUILD "Enable build time profiling (-ftime-trace)" OFF)

# -----------------------------------------------------------------------------
# C++
# -----------------------------------------------------------------------------

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# -----------------------------------------------------------------------------
# Compiler flags
# -----------------------------------------------------------------------------

add_definitions(-D_USE_MATH_DEFINES)

if(MSVC)

  message(STATUS "Compiler: MSVC")

  add_definitions(-DHMAP_MSVC)

  # do NOT pass /Od here: add_compile_options lands after the per-config flags
  # on the MSVC command line, where it overrides Release's /O2 (last flag wins),
  # silently disabling optimization in every Release build
  add_compile_options(/W4 "$<$<NOT:$<CONFIG:Debug>>:/fp:fast>")

  add_compile_definitions(M_PI=3.14159265358979323846)

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")

  message(STATUS "Compiler: Clang")

  add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    $<$<CONFIG:Debug>:-Og>
    $<$<NOT:$<CONFIG:Debug>>:-Ofast>
    $<$<NOT:$<CONFIG:Debug>>:-ffast-math>
    -pthread
    -fPIC
    -Wno-deprecated-enum-enum-conversion
    -Wno-header-guard
    -Wno-gnu-zero-variadic-macro-arguments)

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

  message(STATUS "Compiler: GCC")

  add_compile_options(
    -Wall
    -Wextra
    $<$<CONFIG:Debug>:-Og>
    $<$<NOT:$<CONFIG:Debug>>:-Ofast>
    $<$<NOT:$<CONFIG:Debug>>:-ffast-math>
    -pthread
    -fPIC
    -Wno-free-nonheap-object
    -Wno-deprecated-enum-enum-conversion)

endif()

# -----------------------------------------------------------------------------
# Profiling
# -----------------------------------------------------------------------------

if(HIGHMAP_ENABLE_PROFILE_BUILD)
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag("-ftime-trace" HIGHMAP_HAS_FTIME_TRACE)
  if(HIGHMAP_HAS_FTIME_TRACE)
    message(STATUS "HIGHMAP_ENABLE_PROFILE_BUILD: -ftime-trace enabled")
    add_compile_options(-ftime-trace)
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message(STATUS "HIGHMAP_ENABLE_PROFILE_BUILD: -ftime-report enabled")
    add_compile_options(-ftime-report)
  endif()
endif()

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# -----------------------------------------------------------------------------
# Dependencies
# -----------------------------------------------------------------------------

find_package(GSL REQUIRED COMPONENTS gsl gslcblas)
find_package(assimp REQUIRED)
find_package(glm REQUIRED)
find_package(OpenMP REQUIRED)
find_package(OpenCV REQUIRED)
find_package(OpenCL REQUIRED)

include_directories(${OpenCV_INCLUDE_DIRS})

# OpenMP
add_compile_options(${OpenMP_CXX_FLAGS})

# OpenCL
add_compile_definitions(CL_HPP_MINIMUM_OPENCL_VERSION=120
                        CL_HPP_TARGET_OPENCL_VERSION=120)

# Logging
if(HIGHMAP_ENABLE_LOGS)
  add_compile_definitions(HIGHMAP_ENABLE_LOGS=1)
else()
  add_compile_definitions(HIGHMAP_ENABLE_LOGS=0)
endif()

# -----------------------------------------------------------------------------
# Status output
# -----------------------------------------------------------------------------

message(STATUS "┌─ HighMap configuration ─────────────")
message(STATUS "│ HIGHMAP_ENABLE_OPENCL:              ${HIGHMAP_ENABLE_OPENCL}")
message(STATUS "│ HIGHMAP_ENABLE_LOGS:                ${HIGHMAP_ENABLE_LOGS}")
message(
  STATUS "│ HIGHMAP_ENABLE_EXAMPLES:            ${HIGHMAP_ENABLE_EXAMPLES}")
message(STATUS "│ HIGHMAP_ENABLE_TESTS:               ${HIGHMAP_ENABLE_TESTS}")
message(
  STATUS "│ HIGHMAP_ENABLE_BENCHMARKS:          ${HIGHMAP_ENABLE_BENCHMARKS}")
message(
  STATUS "│ HIGHMAP_ENABLE_PROFILE_BUILD:       ${HIGHMAP_ENABLE_PROFILE_BUILD}")
message(STATUS "│ CMAKE_CXX_STANDARD:                 ${CMAKE_CXX_STANDARD}")
message(STATUS "│ CMAKE_BUILD_TYPE:                   ${CMAKE_BUILD_TYPE}")
message(STATUS "│ CMAKE_CXX_COMPILER_ID:              ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "└─────────────────────────────────────")

# -----------------------------------------------------------------------------
# External libs
# -----------------------------------------------------------------------------

add_subdirectory(external)

# -----------------------------------------------------------------------------
# Library
# -----------------------------------------------------------------------------

add_subdirectory(HighMap)

# -----------------------------------------------------------------------------
# Examples
# -----------------------------------------------------------------------------

if(HIGHMAP_ENABLE_EXAMPLES)
  add_subdirectory(examples)
endif()

# -----------------------------------------------------------------------------
# Tests
# -----------------------------------------------------------------------------

if(HIGHMAP_ENABLE_TESTS)
  message(STATUS "HighMap tests enabled")
  add_subdirectory(tests)
endif()

# -----------------------------------------------------------------------------
# Benchmarks
# -----------------------------------------------------------------------------

if(HIGHMAP_ENABLE_BENCHMARKS)
  message(STATUS "HighMap benchmarks enabled")
  add_subdirectory(benchmarks)
endif()
