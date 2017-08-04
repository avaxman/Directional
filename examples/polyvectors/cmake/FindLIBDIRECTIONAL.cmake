# - Try to find the LIBDIRECTIONAL library
# Once done this will define
#
#  LIBDIRECTIONAL_FOUND - system has LIBDIRECTIONAL
#  LIBDIRECTIONAL_INCLUDE_DIR - **the** LIBDIRECTIONAL include directory
if(NOT LIBDIRECTIONAL_FOUND)

FIND_PATH(LIBDIRECTIONAL_INCLUDE_DIR directional/representative_to_raw.h
   ${PROJECT_SOURCE_DIR}/../../include
   ${PROJECT_SOURCE_DIR}/../include
   ${PROJECT_SOURCE_DIR}/include
   /usr/include
   /usr/local/include
)


if(LIBDIRECTIONAL_INCLUDE_DIR)
   set(LIBDIRECTIONAL_FOUND TRUE)
   set(LIBDIRECTIONAL_INCLUDE_DIRS ${LIBDIRECTIONAL_INCLUDE_DIR})
endif()

endif()
