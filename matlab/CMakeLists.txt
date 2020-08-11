cmake_minimum_required(VERSION 3.1)
project(Directional_MEX_Files)

### conditionally compile certain modules depending on libraries found on the system
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/../cmake)

set(TEST_SHARED_PATH ${CMAKE_CURRENT_SOURCE_DIR}/../tutorials/shared CACHE PATH "location of resources for testing")

### Adding libIGL and Directional: choose the path to your local copy
find_package(LIBIGL REQUIRED QUIET)
include(Directional)


set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake) # add FindMatlab module
set(VERBOSE Off CACHE BOOL "Mex files output extra information")
set(COPY_MEX_OUTPUT On CACHE BOOL "Flags that compiled MEX files should be copied")
set(MEX_OUTPUT_COPY_DIR ${CMAKE_SOURCE_DIR}/dist CACHE PATH "Mex files output extra information")
message(STATUS "Starting to look for matlab")

find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY ENG_LIBRARY)
add_definitions(/DMATLAB_MEX_FILE) #define matlab macros

#Add the verbosity preprocessor definition
if(VERBOSE)
add_definitions(/DVERBOSE)
endif()

#Disable this when compiling 64 bit, otherwise Matlab starts complaining
#add_definitions(/DMX_COMPAT_32)

function(CopyMexFiles Target FuncName)
    if(COPY_MEX_OUTPUT)
        add_custom_command(TARGET ${Target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${Target}> ${CMAKE_CURRENT_SOURCE_DIR}/${Target}.m ${MEX_OUTPUT_COPY_DIR})
    endif()
endfunction()

#Setup the mex target function

function(CreateMexTarget TargetName)
    message(STATUS "Mex Target ${TargetName}")
    # Set filename
    set(CPP_FILE_NAME ${TargetName}.cpp)
    #INCLUDE_DIRECTORIES(${MATLAB_INCLUDE_DIRs})
    add_library(${TargetName} SHARED ${CPP_FILE_NAME} ../mexHelpers.h ../Matlabdef.def)

    # Include Directional, IGL and necessary matlab stuff
    target_include_directories(${TargetName} PRIVATE ${DIRECTIONAL_SOURCE_DIR})
    target_include_directories(${TargetName} PRIVATE ${Matlab_INCLUDE_DIRS})
    target_link_libraries(${TargetName} igl::core ${Matlab_MEX_LIBRARY} ${Matlab_MX_LIBRARY})

    target_link_libraries(${TargetName} ${MATLAB_LIBRARIES})

    # 32-bit or 64-bit mex, windows or not
    if(WIN32)
        if (CMAKE_CL_64)
            SET_TARGET_PROPERTIES(${TargetName} PROPERTIES SUFFIX .mexw64)
        else()
            SET_TARGET_PROPERTIES(${TargetName} PROPERTIES SUFFIX .mexw32)
        endif()
    else()
        if (CMAKE_SIZEOF_VOID_P MATCHES "8")
            SET_TARGET_PROPERTIES(${TargetName} PROPERTIES SUFFIX .mexa64 PREFIX "")
        else()
            SET_TARGET_PROPERTIES(${TargetName} PROPERTIES SUFFIX .mexglx PREFIX "")
        endif()
    endif()

    CopyMexFiles(${TargetName} ${TargetName})
endfunction()


IF(MATLAB_FOUND)
    message(STATUS "MATLAB Found, MATLAB MEX will be compiled.")
    add_subdirectory(mexhelpers)
    add_subdirectory(mex)
ELSE(MATLAB_FOUND)
    MESSAGE("MATLAB not found...nothing will be built.")
ENDIF(MATLAB_FOUND)
