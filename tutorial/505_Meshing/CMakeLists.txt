cmake_minimum_required(VERSION 3.20)
project(505_Meshing)

add_executable(${PROJECT_NAME}_bin main.cpp)

target_link_libraries(${PROJECT_NAME}_bin
    PUBLIC
        tutorials
        polyscope
)

# Link GMP only if enabled and found
if(USE_GMP AND GMP_LIBRARIES)
    target_compile_definitions(${PROJECT_NAME}_bin PUBLIC USE_GMP)
    target_include_directories(${PROJECT_NAME}_bin PUBLIC ${GMP_INCLUDE_DIR})
    target_link_libraries(${PROJECT_NAME}_bin PUBLIC ${GMP_LIBRARIES})
endif()
