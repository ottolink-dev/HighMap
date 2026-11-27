cmake_minimum_required(VERSION 3.22)

project(
  highmap
  VERSION 0.0.1
  LANGUAGES CXX
  DESCRIPTION
    "C++ library to generate two-dimensional terrain heightmaps for software rendering or video games."
)

file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

if(NOT HIGHMAP_ENABLE_OPENCV)
  list(FILTER SOURCES EXCLUDE REGEX ".*/(opencv_wrapper|opencv_writer)\\.cpp$")
endif()

if(NOT HIGHMAP_ENABLE_LIGHTUSD)
  list(FILTER SOURCES EXCLUDE REGEX ".*/export_usd\\.cpp$")
endif()

add_library(${PROJECT_NAME} STATIC ${SOURCES})

if(HIGHMAP_ENABLE_TRACY)
  find_package(Tracy CONFIG REQUIRED)
  # Public because the tile dispatcher is instantiated by library consumers.
  target_link_libraries(${PROJECT_NAME} PUBLIC Tracy::TracyClient)
  target_compile_definitions(${PROJECT_NAME} PUBLIC HIGHMAP_ENABLE_TRACY)
endif()

target_include_directories(
  ${PROJECT_NAME}
  PUBLIC ${PROJECT_SOURCE_DIR}/include
  PRIVATE ${PROJECT_SOURCE_DIR}/include/internal)

set_target_properties(${PROJECT_NAME} PROPERTIES VERSION ${PROJECT_VERSION})

target_compile_definitions(
  ${PROJECT_NAME} INTERFACE GLM_FORCE_SWIZZLE GLM_FORCE_INLINE
                            GLM_ENABLE_EXPERIMENTAL _USE_MATH_DEFINES)

target_link_libraries(
  ${PROJECT_NAME}
  PUBLIC
    assimp::assimp
    delaunator::delaunator
    glm::glm
    GSL::gsl
    GSL::gslcblas
    dkm::dkm
    hmm::hmm
    libnpy::libnpy
    nanoflann::nanoflann
    mixbox
    NoiseLib
    terrain-descriptors
    OpenMP::OpenMP_CXX
    $<$<BOOL:${HIGHMAP_ENABLE_OPENCV}>:${OpenCV_LIBS}>
    point_sampler
    clwrapper
    ${OpenCL_LIBRARIES}
    $<$<BOOL:${HIGHMAP_ENABLE_LIGHTUSD}>:lightusd_static>
    OpenEXR::OpenEXR
    TIFF::TIFF
    nlohmann_json::nlohmann_json
    $<IF:$<TARGET_EXISTS:zstd::libzstd_shared>,zstd::libzstd_shared,zstd::libzstd_static>
  PRIVATE nn-c::nn-c)

if(HIGHMAP_ENABLE_PCH)
  message(STATUS "Precompiled Headers (PCH) enabled for target ${PROJECT_NAME}")
  target_precompile_headers(
    ${PROJECT_NAME}
    PRIVATE
    <vector>
    <string>
    <memory>
    <cmath>
    <algorithm>
    <functional>
    <map>
    <unordered_map>
    <iostream>
    <glm/glm.hpp>)
endif()
