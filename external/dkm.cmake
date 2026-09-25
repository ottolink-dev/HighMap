include(FetchContent)

FetchContent_Declare(
  dkm
  GIT_REPOSITORY https://github.com/genbattle/dkm.git
  GIT_TAG master
  GIT_SHALLOW TRUE)

FetchContent_GetProperties(dkm)
if(NOT dkm_POPULATED)
  FetchContent_Populate(dkm)
  add_library(dkm INTERFACE)
  add_library(dkm::dkm ALIAS dkm)
  target_include_directories(dkm INTERFACE ${dkm_SOURCE_DIR}/include)
endif()
