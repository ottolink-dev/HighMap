cmake_minimum_required(VERSION 3.22)

project(
  highmap
  VERSION 0.0.1
  LANGUAGES CXX
  DESCRIPTION
    "C++ library to generate two-dimensional terrain heightmaps for software rendering or video games."
)

file(GLOB_RECURSE SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_library(${PROJECT_NAME} STATIC ${SOURCES})

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
  assimp::assimp
  delaunator::delaunator
  glm::glm
  GSL::gsl
  GSL::gslcblas
  dkm::dkm
  FastNoiseLite::FastNoiseLite
  hmm::hmm
  libnpy::libnpy
  macro-logger::macro-logger
  nanoflann::nanoflann
  nn-c::nn-c
  NoiseLib
  terrain-descriptors
  ${OpenCV_LIBS}
  point_sampler
  clwrapper
  ${OpenCL_LIBRARIES})

# OpenMP is optional on macOS (see root CMakeLists)
if(OpenMP_CXX_FOUND)
  target_link_libraries(${PROJECT_NAME} OpenMP::OpenMP_CXX)
endif()
