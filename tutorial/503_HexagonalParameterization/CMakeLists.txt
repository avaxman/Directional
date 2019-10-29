cmake_minimum_required(VERSION 2.8.12)
project(503_HexagonalParameterization)

set(LIBIGL_WITH_PNG CACHE BOOL ON)

add_executable(${PROJECT_NAME}_bin main.cpp)
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw igl::png tutorials)
