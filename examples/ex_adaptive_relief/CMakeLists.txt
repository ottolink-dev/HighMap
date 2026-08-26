cmake_minimum_required(VERSION 3.22)

project(ex_adaptive_relief)

add_executable(${PROJECT_NAME} ex_adaptive_relief.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE highmap)
