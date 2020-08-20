cmake_minimum_required(VERSION 3.1)
project(601_SubdivisionFields)

set(CMAKE_CXX_STANDARD 14)

add_executable(${PROJECT_NAME}_bin main.cpp)
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw tutorials)
