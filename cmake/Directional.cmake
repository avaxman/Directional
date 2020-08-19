cmake_minimum_required(VERSION 3.1)



################################################################################

### Configuration

set(DIRECTIONAL_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")
set(DIRECTIONAL_SOURCE_DIR "${DIRECTIONAL_ROOT}/include")
set(DIRECTIONAL_EXTERNAL "${DIRECTIONAL_ROOT}/external")

set(LIBIGL_ROOT "${DIRECTIONAL_EXTERNAL}/libigl")
set(LIBIGL_SOURCE_DIR "${LIBIGL_ROOT}/include")
set(LIBIGL_EXTERNAL "${LIBIGL_ROOT}/external")

include_directories(${DIRECTIONAL_SOURCE_DIR})
