project(nn-c)

set(NNC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/nn-c/nn)

# Configure config.h (actually generate a dummy one)
if(NOT EXISTS ${CMAKE_BINARY_DIR}/nn-c/config.h)
  file(WRITE ${CMAKE_BINARY_DIR}/nn-c/config.h "")
endif()

file(GLOB_RECURSE NNC_SRC ${NNC_DIR}/*.c)

set(NNC_INCLUDE_DIR ${NNC_DIR})

add_library(${PROJECT_NAME} STATIC ${NNC_SRC})
add_library(${PROJECT_NAME}::${PROJECT_NAME} ALIAS ${PROJECT_NAME})

target_compile_definitions(${PROJECT_NAME} PRIVATE TRILIBRARY=1)
target_compile_definitions(${PROJECT_NAME} PRIVATE NO_TIMER=1)

target_include_directories(${PROJECT_NAME} PUBLIC ${NNC_INCLUDE_DIR}
                                                  ${CMAKE_BINARY_DIR}/nn-c)
