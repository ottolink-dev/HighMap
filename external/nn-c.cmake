include(FetchContent)

FetchContent_Declare(
  nn-c
  GIT_REPOSITORY https://github.com/sakov/nn-c.git
  GIT_TAG master
  GIT_SHALLOW TRUE)

FetchContent_GetProperties(nn-c)
if(NOT nn-c_POPULATED)
  FetchContent_Populate(nn-c)

  set(NNC_DIR ${nn-c_SOURCE_DIR}/nn)

  # Configure config.h (actually generate a dummy one)
  if(NOT EXISTS ${CMAKE_BINARY_DIR}/nn-c/config.h)
    file(WRITE ${CMAKE_BINARY_DIR}/nn-c/config.h "")
  endif()

  file(GLOB_RECURSE NNC_SRC ${NNC_DIR}/*.c)

  add_library(nn-c STATIC ${NNC_SRC})
  add_library(nn-c::nn-c ALIAS nn-c)

  target_compile_definitions(nn-c PRIVATE TRILIBRARY=1)
  target_compile_definitions(nn-c PRIVATE NO_TIMER=1)

  target_include_directories(nn-c PUBLIC ${NNC_DIR} ${CMAKE_BINARY_DIR}/nn-c)
endif()
