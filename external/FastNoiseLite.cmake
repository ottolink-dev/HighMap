include(FetchContent)

FetchContent_Declare(
  FastNoiseLite
  GIT_REPOSITORY https://github.com/Auburn/FastNoiseLite.git
  GIT_TAG v1.1.1
  GIT_SHALLOW TRUE)

FetchContent_GetProperties(FastNoiseLite)
if(NOT fastnoiselite_POPULATED)
  FetchContent_Populate(FastNoiseLite)
  add_library(FastNoiseLite INTERFACE)
  add_library(FastNoiseLite::FastNoiseLite ALIAS FastNoiseLite)
  target_include_directories(FastNoiseLite
                             INTERFACE ${fastnoiselite_SOURCE_DIR}/Cpp)
endif()
