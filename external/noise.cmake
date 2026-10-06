include(FetchContent)

FetchContent_Declare(
  Noise
  GIT_REPOSITORY https://github.com/otto-link/Noise.git
  GIT_TAG master
  GIT_SHALLOW TRUE
  SOURCE_SUBDIR NoiseLib)

FetchContent_MakeAvailable(Noise)
