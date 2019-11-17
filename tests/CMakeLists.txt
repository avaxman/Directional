cmake_minimum_required(VERSION 3.1)
project(Directional_tests)

### conditionally compile certain modules depending on libraries found on the system
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)

### libIGL options:
option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_OPENGL      "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW"           ON)
option(LIBIGL_WITH_EMBREE      "Use Embree"         ON)

# Path to models to perform tests on
set(TEST_SHARED_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../tutorials/shared CACHE PATH "location of resources for testing")


# Add google test
set(GTEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../external/googletest)
set(GTEST_INCLUDE_DIR ${GTEST_DIR}/googletest/include)
add_subdirectory(${GTEST_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/gtest_bin)

### Adding libIGL and Directional: choose the path to your local copy
find_package(LIBIGL REQUIRED QUIET)
include(Directional)

# Enables testing
enable_testing()

# The libraries to link against.
set(LINK_LIBS igl::core igl::opengl_glfw gtest)

# Commutation tests
add_executable(CommutationsTests CommutationTests.cpp)
target_include_directories(CommutationsTests INTERFACE ${DIRECTIONAL_SOURCE_DIR} ${GTEST_INCLUDE_DIR})
target_link_libraries(CommutationsTests PUBLIC ${LINK_LIBS})
target_sources(CommutationsTests PRIVATE Helpers.h)
target_compile_definitions(CommutationsTests PUBLIC "-DTEST_SHARED_PATH=\"${TEST_SHARED_PATH}\"")
add_test(NAME CommutationsTest COMMAND CommutationsTests)

# Operator tests
add_executable(OperatorTests OperatorTests.cpp)
target_include_directories(OperatorTests INTERFACE ${DIRECTIONAL_SOURCE_DIR} ${GTEST_INCLUDE_DIR})
target_link_libraries(OperatorTests PUBLIC ${LINK_LIBS})
target_sources(OperatorTests PRIVATE Helpers.h DCELFixture.h)
target_compile_definitions(OperatorTests PUBLIC "-DTEST_SHARED_PATH=\"${TEST_SHARED_PATH}\"")
add_test(NAME OperatorTests COMMAND OperatorTests)

