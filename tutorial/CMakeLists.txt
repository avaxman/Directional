cmake_minimum_required(VERSION 3.20)
project(Directional_tutorials)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)  # Optional: use -std=c++20 instead of -std=gnu++20

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Build type")
endif()

option(DIRECTIONAL_VIEWER "Directional Viewer" ON)
option(USE_GMP "Enable GMP support if found" ON)

message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

set(DIRECTIONAL_ROOT "${CMAKE_CURRENT_LIST_DIR}/../")
set(DIRECTIONAL_SOURCE_DIR "${DIRECTIONAL_ROOT}/include")
set(DIRECTIONAL_EXTERNAL "${DIRECTIONAL_ROOT}/external")

message(STATUS "DIRECTIONAL_ROOT: ${DIRECTIONAL_ROOT}")

include_directories(${DIRECTIONAL_SOURCE_DIR})
include_directories(${DIRECTIONAL_EXTERNAL}/eigen)

# Polyscope
add_subdirectory(${DIRECTIONAL_EXTERNAL}/polyscope "polyscope")

# Tutorial data path and output path (create dirs if needed)
set(TUTORIAL_DATA_PATH "${CMAKE_CURRENT_SOURCE_DIR}/data" CACHE PATH "location of tutorial data resources")
set(TUTORIAL_OUTPUT_PATH "${CMAKE_CURRENT_SOURCE_DIR}/output" CACHE PATH "location of tutorial outputs")

file(MAKE_DIRECTORY "${TUTORIAL_DATA_PATH}")
file(MAKE_DIRECTORY "${TUTORIAL_OUTPUT_PATH}")

# Set runtime output dir and create
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")
file(MAKE_DIRECTORY "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")

# Interface library for tutorial targets
add_library(tutorials INTERFACE)

target_compile_definitions(tutorials INTERFACE
  TUTORIAL_DATA_PATH=\"${TUTORIAL_DATA_PATH}\"
  TUTORIAL_OUTPUT_PATH=\"${TUTORIAL_OUTPUT_PATH}\"
  GLM_ENABLE_EXPERIMENTAL
)


target_include_directories(tutorials INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${POLYSCOPE_SOURCE_DIR}/include
    ${POLYSCOPE_SOURCE_DIR}/deps/imgui/imgui
    ${POLYSCOPE_SOURCE_DIR}/deps/glm
)

# ===== GMP Handling (Cross-platform) =====
if(USE_GMP)
    message(STATUS "Attempting to find GMP...")

    find_path(GMP_INCLUDE_DIR
        NAMES gmp.h
        HINTS
            $ENV{GMP_DIR}/include
            ${CMAKE_PREFIX_PATH}/include
            /usr/include
            /usr/local/include
            "C:/Program Files/GMP/include"
    )

    find_library(GMP_LIBRARIES
        NAMES gmp libgmp
        HINTS
            $ENV{GMP_DIR}/lib
            ${CMAKE_PREFIX_PATH}/lib
            /usr/lib
            /usr/local/lib
            "C:/Program Files/GMP/lib"
    )

    if(GMP_INCLUDE_DIR AND GMP_LIBRARIES)
        message(STATUS "GMP found:")
        message(STATUS "  Include: ${GMP_INCLUDE_DIR}")
        message(STATUS "  Library: ${GMP_LIBRARIES}")
        add_compile_definitions(USE_GMP)
        include_directories(${GMP_INCLUDE_DIR})
        link_libraries(${GMP_LIBRARIES})
    else()
        message(WARNING "GMP not found. Disabling USE_GMP.")
        set(USE_GMP OFF CACHE BOOL "Enable GMP support if found" FORCE)
    endif()
endif()

if(USE_GMP)
    add_compile_definitions(USE_GMP_ENABLED)
endif()

# ===== Chapters =====
if(DIRECTIONAL_VIEWER)
  add_subdirectory("101_GlyphRendering")
  add_subdirectory("102_DiscreteTangentBundles")
  add_subdirectory("103_StreamlineTracing")
  add_subdirectory("104_FaceVertexEdgeData")
  add_subdirectory("105_Sparsity")
  add_subdirectory("106_PrincipalDirections")
  add_subdirectory("201_PrincipalMatching")
  add_subdirectory("202_Sampling")
  add_subdirectory("203_Combing")
  add_subdirectory("301_PowerFields")
  add_subdirectory("302_PolyVectors")
  add_subdirectory("303_GinzburgLandau")
  add_subdirectory("304_IntegrableFields")
  add_subdirectory("401_IndexPrescription")
  add_subdirectory("501_SeamlessIntegration")
  add_subdirectory("502_DifferentOrders")
  add_subdirectory("503_SeamsSingsRounding")
  add_subdirectory("504_LinearReductions")
  add_subdirectory("505_Meshing")
  add_subdirectory("601_FaceBasedFEM")
  add_subdirectory("602_DiscreteExteriorCalculus")
  add_subdirectory("603_HodgeDecomposition")
  add_subdirectory("604_HodgeBoundary")
  add_subdirectory("605_HarmonicFields")
endif()
