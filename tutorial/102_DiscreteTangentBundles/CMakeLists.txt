cmake_minimum_required(VERSION 3.16)
project(102_DiscreteTangentBundles)

set(CMAKE_CXX_STANDARD 14)

add_executable(${PROJECT_NAME}_bin main.cpp)
target_link_libraries(${PROJECT_NAME}_bin PUBLIC tutorials)
