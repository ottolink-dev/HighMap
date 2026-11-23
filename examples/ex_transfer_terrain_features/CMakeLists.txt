project(ex_transfer_terrain_features)

add_executable(${PROJECT_NAME} ex_transfer_terrain_features.cpp)
target_link_libraries(${PROJECT_NAME} highmap)
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                 ${CMAKE_BINARY_DIR}/bin)
