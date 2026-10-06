include(FetchContent)

FetchContent_Declare(
  dkm
  GIT_REPOSITORY https://github.com/genbattle/dkm.git
  GIT_TAG master
  GIT_SHALLOW TRUE
  SOURCE_SUBDIR include)

FetchContent_MakeAvailable(dkm)

if(NOT TARGET dkm)
  add_library(dkm INTERFACE)
  add_library(dkm::dkm ALIAS dkm)
  target_include_directories(dkm INTERFACE ${dkm_SOURCE_DIR}/include)
endif()
