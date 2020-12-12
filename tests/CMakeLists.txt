cmake_minimum_required(VERSION 3.11)
project(Directional_test)

find_package(Catch2 CONFIG REQUIRED)

if(MSVC)
    list(APPEND CMAKE_CXX_FLAGS "/MP")
endif()

add_executable(Directional_tests
    # Main file
    test_main.cpp
    # Separate test cases
        test_fem_operators.cpp
        test_shm_commutation.cpp
        test_includes.h
        )

target_compile_definitions(Directional_tests PRIVATE -DDIRECTIONAL_TEST_MODELS_DIR="${CMAKE_CURRENT_LIST_DIR}/models")

target_link_libraries(Directional_tests
        PRIVATE
        Directional
        Catch2::Catch2
        )

include(CTest)
include(Catch)
catch_discover_tests(Directional_tests)