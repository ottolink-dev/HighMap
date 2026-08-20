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
option(HIGHMAP_ENABLE_EXAMPLES "" ON)
option(HIGHMAP_ENABLE_TESTS "" ON)
option(HIGHMAP_ENABLE_BENCHMARKS "" ON)

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

# MATCHES so that AppleClang (macOS default compiler) is covered too:
# STREQUAL "Clang" silently skipped this branch on macOS, which left Apple
# Silicon builds without any optimization flags (same bug class as the MSVC
# /Od issue fixed in 7237678)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "^(Clang|AppleClang)$")

  message(STATUS "Compiler: Clang (${CMAKE_CXX_COMPILER_ID})")

  add_compile_options(
    -Wall
    -Wextra
    -Wpedantic
    -pthread
    -fPIC
    # do NOT apply optimization flags for Debug builds; mirrors the MSVC
    # branch ($<$<NOT:$<CONFIG:Debug>>:...>); on single-config generators
    # without an explicit build type the flags apply by default; note that
    # -Ofast is deprecated in recent clang, "-O3 -ffast-math" is equivalent
    $<$<NOT:$<CONFIG:Debug>>:-O3>
    $<$<NOT:$<CONFIG:Debug>>:-ffast-math>
    -Wno-deprecated-enum-enum-conversion
    -Wno-header-guard
    -Wno-gnu-zero-variadic-macro-arguments)

  add_compile_definitions(LOG_LEVEL=3)

elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

  message(STATUS "Compiler: GCC")

  add_compile_options(
    -Wall
    -Wextra
    -O3
    -Ofast
    -ffast-math
    -pthread
    -fPIC
    -Wno-free-nonheap-object
    -Wno-deprecated-enum-enum-conversion)

  add_compile_definitions(LOG_LEVEL=3)

endif()

# -----------------------------------------------------------------------------
# Architecture tuning
# -----------------------------------------------------------------------------

# Apple Silicon: target the M1 baseline so the compiler can exploit the full
# NEON/FMA instruction set available on every Apple Silicon chip (M1 and up)
if(APPLE AND CMAKE_SYSTEM_PROCESSOR MATCHES "arm64")

  message(STATUS "CPU tuning: Apple Silicon (-mcpu=apple-m1)")

  add_compile_options($<$<NOT:$<CONFIG:Debug>>:-mcpu=apple-m1>)

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

# Homebrew assimp may embed an absolute libz.tbd path pointing at a
# CommandLineTools SDK that no longer exists; rebind ZLIB::ZLIB (and the
# assimp target link libraries) to the zlib from the active SDK when that
# happens (macOS only)
if(APPLE)
  find_library(
    HMAP_SDK_ZLIB z
    PATHS "${CMAKE_OSX_SYSROOT}/usr/lib"
    NO_DEFAULT_PATH)

  if(HMAP_SDK_ZLIB AND TARGET ZLIB::ZLIB)
    get_target_property(_zlib_loc ZLIB::ZLIB IMPORTED_LOCATION)

    if(_zlib_loc AND NOT EXISTS "${_zlib_loc}")
      set_target_properties(ZLIB::ZLIB PROPERTIES IMPORTED_LOCATION
                                                  "${HMAP_SDK_ZLIB}")
    endif()
  endif()

  if(HMAP_SDK_ZLIB AND TARGET assimp::assimp)
    get_target_property(_assimp_libs assimp::assimp INTERFACE_LINK_LIBRARIES)
    set(_assimp_libs_fixed "")

    foreach(_lib ${_assimp_libs})
      if(_lib MATCHES "CommandLineTools/.*libz\\.tbd" AND NOT EXISTS
                                                          "${_lib}")
        list(APPEND _assimp_libs_fixed "${HMAP_SDK_ZLIB}")
      else()
        list(APPEND _assimp_libs_fixed "${_lib}")
      endif()
    endforeach()

    set_target_properties(assimp::assimp PROPERTIES INTERFACE_LINK_LIBRARIES
                                                    "${_assimp_libs_fixed}")
  endif()
endif()
find_package(glm REQUIRED)
find_package(OpenCV REQUIRED)
find_package(OpenCL REQUIRED)

# On macOS the OpenCL C/C++ headers come from keg-only Homebrew packages
# (brew install opencl-headers opencl-clhpp-headers) which are not linked
# into the default include path: point the compiler at them
if(APPLE)
  if(NOT DEFINED HOMEBREW_PREFIX)
    execute_process(
      COMMAND brew --prefix
      OUTPUT_VARIABLE HOMEBREW_PREFIX
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET)
  endif()

  foreach(_keg opencl-headers opencl-clhpp-headers)
    if(EXISTS "${HOMEBREW_PREFIX}/opt/${_keg}/include")
      include_directories("${HOMEBREW_PREFIX}/opt/${_keg}/include")
    endif()
  endforeach()
endif()

include_directories(${OpenCV_INCLUDE_DIRS})

# OpenMP
#
# Apple Clang does not ship OpenMP; on macOS it is provided by the Homebrew
# `libomp` keg (brew install libomp). Point FindOpenMP at it when present.
if(APPLE AND CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  if(NOT DEFINED HOMEBREW_PREFIX)
    execute_process(
      COMMAND brew --prefix
      OUTPUT_VARIABLE HOMEBREW_PREFIX
      OUTPUT_STRIP_TRAILING_WHITESPACE
      ERROR_QUIET)
  endif()

  if(EXISTS "${HOMEBREW_PREFIX}/opt/libomp")
    set(OpenMP_C_FLAGS
        "-Xpreprocessor -fopenmp -I${HOMEBREW_PREFIX}/opt/libomp/include")
    set(OpenMP_CXX_FLAGS
        "-Xpreprocessor -fopenmp -I${HOMEBREW_PREFIX}/opt/libomp/include")
    set(OpenMP_C_LIB_NAMES libomp)
    set(OpenMP_CXX_LIB_NAMES libomp)
    set(OpenMP_libomp_LIBRARY
        "${HOMEBREW_PREFIX}/opt/libomp/lib/libomp.dylib")
  endif()
endif()

find_package(OpenMP QUIET)

if(NOT OpenMP_CXX_FOUND)
  if(APPLE)
    # Only a few kernels use OpenMP (flow accumulation, drainage basins,
    # local metrics); tile-based threading keeps working without it
    message(
      WARNING
        "OpenMP not found: `brew install libomp` to enable it. Building "
        "without OpenMP support.")
    # without -fopenmp, clang ignores `#pragma omp` directives: silence the
    # associated warnings
    add_compile_options(-Wno-unknown-pragmas)
  else()
    message(FATAL_ERROR "OpenMP is required on this platform")
  endif()
endif()

# OpenMP
if(OpenMP_CXX_FOUND)
  # OpenMP_CXX_FLAGS is a space-separated string: split it so that
  # add_compile_options receives individual tokens (matters for the macOS
  # libomp hint "-Xpreprocessor -fopenmp -I...")
  separate_arguments(_openmp_cxx_flags UNIX_COMMAND "${OpenMP_CXX_FLAGS}")
  add_compile_options(${_openmp_cxx_flags})
endif()

# OpenCL
add_compile_definitions(CL_HPP_MINIMUM_OPENCL_VERSION=120
                        CL_HPP_TARGET_OPENCL_VERSION=120)

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
