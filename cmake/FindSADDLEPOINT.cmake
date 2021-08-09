# - Try to find the SaddlePoint library
# Once done this will define
#
#  SADDLEPOINT_FOUND - system has SaddlePoint
#  SADDLEPOINT_INCLUDE_DIR - **the** SaddlePoint include directory
#  SADDLEPOINT_INCLUDE_DIRS - SaddlePoint include directories
#  SADDLEPOINTL_SOURCES - the SaddlePoint source files
if(NOT SADDLEPOINT_FOUND)
message("hello")

FIND_PATH(SADDLEPOINT_INCLUDE_DIR SaddlePoint/LMSolver.h
   ${PROJECT_SOURCE_DIR}/../../../include
   ${PROJECT_SOURCE_DIR}/../../include
   ${PROJECT_SOURCE_DIR}/../include
   ${PROJECT_SOURCE_DIR}/include
   /usr/include
   /usr/local/include
)

if(SADDLEPOINT_INCLUDE_DIR)
   set(SADDLEPOINT_FOUND TRUE)
   set(SADDLEPOINT_INCLUDE_DIRS ${SADDLEPOINT_INCLUDE_DIR})
endif()

endif()
