cmake_minimum_required(VERSION 3.1)
project(604_Advection)

set(CMAKE_CXX_STANDARD 14)

add_executable(${PROJECT_NAME}_bin main.cpp)
target_link_libraries(${PROJECT_NAME}_bin Directional igl::core igl::opengl_glfw tutorials)
