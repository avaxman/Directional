cmake_minimum_required(VERSION 3.1)



################################################################################

### Configuration

set(LIBDIRECTIONAL_ROOT "${CMAKE_CURRENT_LIST_DIR}/..")
set(LIBDIRECTIONAL_SOURCE_DIR "${LIBDIRECTIONAL_ROOT}/include")
set(LIBDIRECTIONAL_EXTERNAL "${LIBDIRECTIONAL_ROOT}/external")

set(LIBIGL_ROOT "${LIBDIRECTIONAL_EXTERNAL}/libigl")
set(LIBIGL_SOURCE_DIR "${LIBIGL_ROOT}/include")
set(LIBIGL_EXTERNAL "${LIBIGL_ROOT}/external")

include_directories(${LIBDIRECTIONAL_SOURCE_DIR})
