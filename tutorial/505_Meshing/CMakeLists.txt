cmake_minimum_required(VERSION 3.1)
project(505_Meshing)

add_executable(${PROJECT_NAME}_bin main.cpp)
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl_glfw  igl::cgal tutorials)
igl_copy_cgal_dll(${PROJECT_NAME}_bin)
