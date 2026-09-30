cmake_minimum_required(VERSION 3.22)
project(ex_hydraulic_fastscape)

add_executable(${PROJECT_NAME} ex_hydraulic_fastscape.cpp)
target_link_libraries(${PROJECT_NAME} PRIVATE highmap)
