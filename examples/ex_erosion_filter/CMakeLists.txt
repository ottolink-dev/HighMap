project(ex_erosion_filter)

add_executable(${PROJECT_NAME} ex_erosion_filter.cpp)

target_link_libraries(${PROJECT_NAME} PRIVATE highmap)
