include(FetchContent)

FetchContent_Declare(
  libnpy
  GIT_REPOSITORY https://github.com/llohse/libnpy.git
  GIT_TAG v1.0.1
  GIT_SHALLOW TRUE)

FetchContent_GetProperties(libnpy)
if(NOT libnpy_POPULATED)
  FetchContent_Populate(libnpy)
  add_library(libnpy INTERFACE)
  add_library(libnpy::libnpy ALIAS libnpy)
  target_include_directories(libnpy INTERFACE ${libnpy_SOURCE_DIR}/include)
endif()
