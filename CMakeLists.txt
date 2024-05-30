cmake_minimum_required(VERSION 3.20)
project(Directional)
set(CMAKE_CXX_STANDARD 17)

add_library(${PROJECT_NAME} INTERFACE include)
target_include_directories(${PROJECT_NAME} INTERFACE include)
