include(FetchContent)

FetchContent_Declare(
  Noise
  GIT_REPOSITORY https://github.com/otto-link/Noise.git
  GIT_TAG master
  GIT_SHALLOW TRUE)

FetchContent_GetProperties(Noise)
if(NOT noise_POPULATED)
  FetchContent_Populate(Noise)
  add_subdirectory(${noise_SOURCE_DIR}/NoiseLib ${noise_BINARY_DIR}/NoiseLib)
endif()
