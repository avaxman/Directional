cmake_minimum_required(VERSION 3.2)

project(Directional_tutorials)

set(CMAKE_CXX_STANDARD 17)


### Choose which chapters to compile
option(DIRECTIONAL_VIEWER "Directional Viewer" ON)
option(DIRECTIONAL_MESHER "Directional Mesher" OFF)

message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

set(DIRECTIONAL_ROOT "${CMAKE_CURRENT_LIST_DIR}/../")
set(DIRECTIONAL_SOURCE_DIR "${DIRECTIONAL_ROOT}/include")
set(DIRECTIONAL_EXTERNAL "${DIRECTIONAL_ROOT}/external")

message(STATUS "DIRECTIONAL_ROOT: ${DIRECTIONAL_ROOT}")

include_directories(${DIRECTIONAL_SOURCE_DIR})
include_directories(${DIRECTIONAL_EXTERNAL}/eigen)
include_directories(${DIRECTIONAL_EXTERNAL}/glfw/include)

# define the path where all the input meshes are available
add_compile_definitions(DATA_PATH="${PROJECT_SOURCE_DIR}/data")

# make a bin folder to host all the executables
make_directory(${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

include(FetchContent)

# Fetch Polyscope
FetchContent_Declare(
    polyscope
    GIT_REPOSITORY https://github.com/nmwsharp/polyscope.git
    GIT_TAG        master
)
FetchContent_MakeAvailable(polyscope)


# Store location of tutorial data directory
set(TUTORIAL_DATA_PATH ${CMAKE_CURRENT_SOURCE_DIR}/data CACHE PATH "location of tutorial data resources")
add_library(tutorials INTERFACE)
target_compile_definitions(tutorials INTERFACE "-DTUTORIAL_DATA_PATH=\"${TUTORIAL_DATA_PATH}\"")
target_include_directories(tutorials INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} ${polyscope_SOURCE_DIR}/include)


# Chapter 1
if(DIRECTIONAL_VIEWER)
  add_subdirectory("101_GlyphRendering")
  add_subdirectory("102_DiscreteTangentBundles")
  add_subdirectory("103_PickingEditing")
  add_subdirectory("104_StreamlineTracing")
  add_subdirectory("105_FaceVertexEdgeData")
  add_subdirectory("106_Sparsity")
  add_subdirectory("201_PrincipalMatching")
  add_subdirectory("202_Sampling")
  add_subdirectory("203_Combing")
  add_subdirectory("301_PowerFields")
  add_subdirectory("302_PolyVectors")
  add_subdirectory("303_PolyCurlReduction")
  add_subdirectory("304_ConjugateFields")
  add_subdirectory("401_IndexPrescription")
  add_subdirectory("501_SeamlessIntegration")
  add_subdirectory("502_DifferentOrders")
  add_subdirectory("503_SeamsSingsRounding")
  add_subdirectory("504_LinearReductions")

  add_subdirectory("601_SubdivisionFields")
endif()

if(DIRECTIONAL_MESHER)
  add_subdirectory("505_Meshing")
endif()
