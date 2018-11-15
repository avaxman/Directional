cmake_minimum_required(VERSION 3.1)
project(Directional_tutorials)
message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

### conditionally compile certain modules depending on libraries found on the system
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)

### libIGL options:
option(LIBIGL_WITH_VIEWER      "Use OpenGL viewer"  ON)
option(LIBIGL_WITH_OPENGL      "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW "Use GLFW"           ON)
option(LIBIGL_WITH_EMBREE      "Use Embree"         ON)


### Adding libIGL and Directional: choose the path to your local copy
find_package(LIBIGL REQUIRED QUIET)
include(Directional)

### Output directories
if(MSVC)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_BINARY_DIR})
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_BINARY_DIR})
else()
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "../")
endif()

### Choose which chapters to compile
option(TUTORIALS_CHAPTER1 "Compile chapter 1" ON)
option(TUTORIALS_CHAPTER2 "Compile chapter 2" ON)
option(TUTORIALS_CHAPTER3 "Compile chapter 3" ON)
option(TUTORIALS_CHAPTER4 "Compile chapter 4" ON)
option(TUTORIALS_CHAPTER5 "Compile chapter 5" ON)


# Store location of tutorial/shared directory
set(TUTORIAL_SHARED_PATH ${CMAKE_CURRENT_SOURCE_DIR}/shared CACHE PATH "location of shared tutorial resources")
add_library(tutorials INTERFACE)
target_compile_definitions(tutorials INTERFACE "-DTUTORIAL_SHARED_PATH=\"${TUTORIAL_SHARED_PATH}\"")
target_include_directories(tutorials INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})


# Chapter 1
if(TUTORIALS_CHAPTER1)
  add_subdirectory("101_GlyphRendering")
  add_subdirectory("102_PickingEditing")
  add_subdirectory("103_StreamlineTracing")
  add_subdirectory("104_DynamicStreamlines")
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
  add_subdirectory("501_SeamlessParameterization")
  add_subdirectory("502_MixedIntegerParameterization")
endif()

# Chapter 6
if(TUTORIALS_CHAPTER5)
  add_subdirectory("601_DifferentialOperators")
  add_subdirectory("602_HodgeDecomposition")
  add_subdirectory("603_HarmonicBasis")
  add_subdirectory("604_MassMatrices")
endif()

# Chapter 7
if(TUTORIALS_CHAPTER5)
  add_subdirectory("701_OneForms")
  add_subdirectory("702_HarmonicFormsDesign")
  add_subdirectory("703_SEC")
endif()

# Chapter 8
if(TUTORIALS_CHAPTER5)
  add_subdirectory("801_OperatorDerivation")
endif()



