cmake_minimum_required(VERSION 3.20)
project(602_DiscreteExteriorCalculus)

add_executable(${PROJECT_NAME}_bin main.cpp)
target_link_libraries(${PROJECT_NAME}_bin PUBLIC tutorials polyscope)
