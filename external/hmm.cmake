include(FetchContent)

FetchContent_Declare(
  hmm
  GIT_REPOSITORY https://github.com/fogleman/hmm.git
  GIT_TAG master
  GIT_SHALLOW TRUE)

FetchContent_GetProperties(hmm)
if(NOT hmm_POPULATED)
  FetchContent_Populate(hmm)

  set(HMM_DIR ${hmm_SOURCE_DIR}/src)
  set(HMM_SRC ${HMM_DIR}/base.cpp ${HMM_DIR}/blur.cpp ${HMM_DIR}/heightmap.cpp
              ${HMM_DIR}/triangulator.cpp)

  add_library(hmm STATIC ${HMM_SRC})
  add_library(hmm::hmm ALIAS hmm)

  target_include_directories(hmm PUBLIC ${HMM_DIR})

  # Provide symlink or copy to support `#include "hmm/src/..."`
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/external_includes)
  if(NOT EXISTS ${CMAKE_BINARY_DIR}/external_includes/hmm)
    file(CREATE_LINK ${hmm_SOURCE_DIR} ${CMAKE_BINARY_DIR}/external_includes/hmm
         SYMBOLIC RESULT link_res)
    if(NOT link_res EQUAL 0)
      file(COPY ${hmm_SOURCE_DIR}/src
           DESTINATION ${CMAKE_BINARY_DIR}/external_includes/hmm)
    endif()
  endif()

  target_include_directories(hmm PUBLIC ${CMAKE_BINARY_DIR}/external_includes)
  target_link_libraries(hmm PUBLIC glm::glm)
endif()
