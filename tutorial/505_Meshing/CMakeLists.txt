cmake_minimum_required(VERSION 3.20)
project(505_Meshing)

add_executable(${PROJECT_NAME}_bin main.cpp)
target_link_libraries(${PROJECT_NAME}_bin PUBLIC tutorials polyscope)
### igl_copy_cgal_dll(${PROJECT_NAME}_bin)
