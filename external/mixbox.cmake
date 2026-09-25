include(FetchContent)

FetchContent_Declare(
  mixbox
  GIT_REPOSITORY https://github.com/scrtwpns/mixbox.git
  GIT_TAG master
  GIT_SHALLOW TRUE)

FetchContent_GetProperties(mixbox)
if(NOT mixbox_POPULATED)
  FetchContent_Populate(mixbox)

  if(EXISTS ${mixbox_SOURCE_DIR}/cpp/mixbox.cpp)
    set(MIXBOX_DIR ${mixbox_SOURCE_DIR}/cpp)
  else()
    set(MIXBOX_DIR ${mixbox_SOURCE_DIR})
  endif()

  add_library(mixbox STATIC ${MIXBOX_DIR}/mixbox.cpp)
  add_library(mixbox::mixbox ALIAS mixbox)

  target_include_directories(mixbox PUBLIC ${MIXBOX_DIR})
endif()
