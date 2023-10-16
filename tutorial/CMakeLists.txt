cmake_minimum_required(VERSION 3.16)
project(Directional_tutorials)
message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

### conditionally compile certain modules depending on libraries found on the system
list(PREPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../cmake)

### Choose which chapters to compile
option(TUTORIALS_CHAPTER1 "Compile chapter 1" ON)
option(TUTORIALS_CHAPTER2 "Compile chapter 2" ON)
option(TUTORIALS_CHAPTER3 "Compile chapter 3" ON)
option(TUTORIALS_CHAPTER4 "Compile chapter 4" ON)
option(TUTORIALS_CHAPTER5 "Compile chapter 5" ON)
option(TUTORIALS_CHAPTER6 "Compile chapter 6" ON)


### Adding libIGL and Directional: choose the path to your local copy
include(Directional)

### Output directories
if(MSVC)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
else()
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")
endif()


# Store location of tutorial/shared directory
set(TUTORIAL_SHARED_PATH ${CMAKE_CURRENT_SOURCE_DIR}/shared CACHE PATH "location of shared tutorial resources")
add_library(tutorials INTERFACE)
target_compile_definitions(tutorials INTERFACE "-DTUTORIAL_SHARED_PATH=\"${TUTORIAL_SHARED_PATH}\"")
target_include_directories(tutorials INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})


# Chapter 1
if(TUTORIALS_CHAPTER1)
  add_subdirectory("101_GlyphRendering")
  add_subdirectory("102_DiscreteTangentBundles")
  add_subdirectory("103_PickingEditing")
  add_subdirectory("104_StreamlineTracing")
  add_subdirectory("105_FaceVertexEdgeData")
  add_subdirectory("106_Sparsity")
endif()

# Chapter 2
if(TUTORIALS_CHAPTER2)
  add_subdirectory("201_PrincipalMatching")
  add_subdirectory("202_Sampling")
  add_subdirectory("203_Combing")
endif()

# Chapter 3
if(TUTORIALS_CHAPTER3)
  add_subdirectory("301_PowerFields")
  add_subdirectory("302_PolyVectors")
  add_subdirectory("303_PolyCurlReduction")
  add_subdirectory("304_ConjugateFields")
endif()

# Chapter 4
if(TUTORIALS_CHAPTER4)
  add_subdirectory("401_IndexPrescription")
endif()

# Chapter 5
if(TUTORIALS_CHAPTER5)
  add_subdirectory("501_SeamlessIntegration")
  add_subdirectory("502_DifferentOrders")
  add_subdirectory("503_SeamsSingsRounding")
  add_subdirectory("504_LinearReductions")
  add_subdirectory("505_Meshing")
endif()

# Chapter 6
if(TUTORIALS_CHAPTER6)
  add_subdirectory("601_SubdivisionFields")
endif()


