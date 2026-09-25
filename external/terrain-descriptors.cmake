include(FetchContent)

FetchContent_Declare(
  terrain-descriptors
  GIT_REPOSITORY https://github.com/otto-link/terrain-descriptors.git
  GIT_TAG main
  GIT_SHALLOW TRUE)

FetchContent_GetProperties(terrain-descriptors)
if(NOT terrain-descriptors_POPULATED)
  FetchContent_Populate(terrain-descriptors)
  add_subdirectory(${terrain-descriptors_SOURCE_DIR}/lib
                   ${terrain-descriptors_BINARY_DIR}/lib)
endif()
