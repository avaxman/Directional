cmake_minimum_required(VERSION 3.20)

project(Directional_tutorials)

set(CMAKE_CXX_STANDARD 17)


### Choose which chapters to compile
option(DIRECTIONAL_VIEWER "Directional Viewer" ON)
option(DIRECTIONAL_MESHER "Directional Mesher" ON)

message(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

set(DIRECTIONAL_ROOT "${CMAKE_CURRENT_LIST_DIR}/../")
set(DIRECTIONAL_SOURCE_DIR "${DIRECTIONAL_ROOT}/include")
set(DIRECTIONAL_EXTERNAL "${DIRECTIONAL_ROOT}/external")

message(STATUS "DIRECTIONAL_ROOT: ${DIRECTIONAL_ROOT}")

include_directories(${DIRECTIONAL_SOURCE_DIR})
include_directories(${DIRECTIONAL_EXTERNAL}/eigen)
message("\n\n == CMAKE recursively building Polyscope\n")
#set(POLYSCOPE_SOURCE_DIR ${DIRECTIONAL_EXTERNAL}/polyscope)

# Create the build directory for Polyscope
#set(POLYSCOPE_BUILD_DIR ${CMAKE_BINARY_DIR}/polyscope-build)
#file(MAKE_DIRECTORY ${POLYSCOPE_BUILD_DIR})

# add_subdirectory(${DIRECTIONAL_EXTERNAL}/polyscope ${CMAKE_CURRENT_BINARY_DIR})
add_subdirectory(${DIRECTIONAL_EXTERNAL}/polyscope "polyscope")

# define the path where all the input meshes are available
add_compile_definitions(DATA_PATH="${PROJECT_SOURCE_DIR}/data")

# make a bin folder to host all the executables
make_directory(${PROJECT_SOURCE_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${PROJECT_SOURCE_DIR}/bin")

if (DIRECTIONAL_MESHER)
  #include(ExternalProject)

  # Set the download and install directories
  #set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

  #et(DEFAULT_GMP_VERSION "6.3.0")

  # Allow the user to specify the GMP version via command line or cache
  #set(GMP_VERSION "${DEFAULT_GMP_VERSION}" CACHE STRING "GMP version")

  # Define the URL based on the specified version
  #set(GMP_URL "https://ftp.gnu.org/gnu/gmp/gmp-6.3.0.tar.bz2")

  # ExternalProject_Add for GMP
  #ExternalProject_Add(
  #        gmp
  #        URL ${GMP_URL}
  #        PREFIX ${EXTERNAL_INSTALL_LOCATION}/gmp
  #        CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env CC=${CMAKE_C_COMPILER} CXX=${CMAKE_CXX_COMPILER} <SOURCE_DIR>/configure --prefix=<INSTALL_DIR>
  #        BUILD_COMMAND ${CMAKE_COMMAND} --build .
  ##        INSTALL_COMMAND ${CMAKE_COMMAND} --install .
  #        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
  #)

  #set(GMP_DIR /usr/local/opt/gmp)

  #find_package(GMP REQUIRED)

  set(GMP_INCLUDE_DIR /usr/local/opt/gmp/include)
  set(GMP_LIBRARIES /usr/local/opt/gmp/lib/libgmp.dylib)

  # Include GMP headers
  include_directories(${GMP_INCLUDE_DIR})

  # Set GMP library and include directories
  #set(GMP_INCLUDE_DIR ${EXTERNAL_INSTALL_LOCATION}/gmp/include)
  #set(GMP_LIB_DIR ${EXTERNAL_INSTALL_LOCATION}/gmp/lib)

  # Add GMP include directory to the include paths
  #include_directories(${GMP_INCLUDE_DIR})

  # Add GMP library directory to link directories
  link_directories(${GMP_LIB_DIR})
endif()


# Store location of tutorial data directory
set(TUTORIAL_DATA_PATH ${CMAKE_CURRENT_SOURCE_DIR}/data CACHE PATH "location of tutorial data resources")
add_library(tutorials INTERFACE)
target_compile_definitions(tutorials INTERFACE "-DTUTORIAL_DATA_PATH=\"${TUTORIAL_DATA_PATH}\"")
target_include_directories(tutorials INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} ${POLYSCOPE_SOURCE_DIR}/include)
target_include_directories(tutorials INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} ${POLYSCOPE_SOURCE_DIR}/deps/imgui/imgui)
target_include_directories(tutorials INTERFACE ${CMAKE_CURRENT_SOURCE_DIR} ${POLYSCOPE_SOURCE_DIR}/deps/glm)
target_compile_definitions(tutorials INTERFACE GLM_ENABLE_EXPERIMENTAL)

# Chapter 1
if(DIRECTIONAL_VIEWER)
  add_subdirectory("101_GlyphRendering")
  add_subdirectory("102_DiscreteTangentBundles")
  add_subdirectory("103_StreamlineTracing")
  add_subdirectory("104_FaceVertexEdgeData")
  add_subdirectory("105_Sparsity")
  add_subdirectory("201_PrincipalMatching")
  add_subdirectory("202_Sampling")
  add_subdirectory("203_Combing")
  add_subdirectory("301_PowerFields")
  add_subdirectory("302_PolyVectors")
  #add_subdirectory("303_PolyCurlReduction")
  #add_subdirectory("304_ConjugateFields")
  add_subdirectory("401_IndexPrescription")
  add_subdirectory("501_SeamlessIntegration")
  add_subdirectory("502_DifferentOrders")
  add_subdirectory("503_SeamsSingsRounding")
  add_subdirectory("504_LinearReductions")
  #add_subdirectory("601_SubdivisionFields")
endif()

if(DIRECTIONAL_MESHER)
  add_subdirectory("505_Meshing")
endif()
