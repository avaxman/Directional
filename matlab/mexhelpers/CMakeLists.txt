cmake_minimum_required(VERSION 3.1)
project(DirectionalMexHelpers)

### conditionally compile certain modules depending on libraries found on the system
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)

set(TEST_SHARED_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../tutorials/shared CACHE PATH "location of resources for testing")

### Adding libIGL and Directional: choose the path to your local copy
find_package(LIBIGL REQUIRED QUIET)
include(Directional)

find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY ENG_LIBRARY)
add_definitions(/DMATLAB_MEX_FILE) #define matlab macros

#Add the verbosity preprocessor definition
if(VERBOSE)
add_definitions(/DVERBOSE)
endif()

set(SOURCES
    FuncArgs.h
    FuncArgs.cpp
    mexHelpers.h
)

add_library(DirectionalMexHelpers ${SOURCES})
target_link_libraries(DirectionalMexHelpers PUBLIC Directional)
target_include_directories(DirectionalMexHelpers PRIVATE ${Matlab_INCLUDE_DIRS})
target_link_libraries(DirectionalMexHelpers PUBLIC igl::core ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY})