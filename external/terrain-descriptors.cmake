include(FetchContent)

FetchContent_Declare(
  terrain-descriptors
  GIT_REPOSITORY https://github.com/otto-link/terrain-descriptors.git
  GIT_TAG main
  GIT_SHALLOW TRUE
  SOURCE_SUBDIR lib)

FetchContent_MakeAvailable(terrain-descriptors)
