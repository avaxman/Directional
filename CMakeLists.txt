cmake_minimum_required(VERSION 3.13)
project(Directional
        LANGUAGES CXX
        VERSION 1.0
        )

# set path to custom cmake modules
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(LIBIGL REQUIRED)

set(HEADERS
    include/directional/angle_bound_frame_fields.h
    include/directional/block_diag.h
    include/directional/columndirectional_to_rawfield.h
    include/directional/combing.h
    include/directional/conjugate_frame_fields.h
    include/directional/ConjugateFFSolverData.h
    include/directional/curl_matching.h
    include/directional/cut_mesh_with_singularities.h
    include/directional/get_directional_subdivision_suite.h
    include/directional/dcel.h
    include/directional/dual_cycles.h
    include/directional/get_P_inverse.h
    include/directional/get_P.h
    include/directional/get_W_inverse.h
    include/directional/get_W.h
    include/directional/dynamic_visualization.h
    include/directional/edge_diamond_mesh.h
    include/directional/effort_to_indices.h
    include/directional/FEM_masses.h
    include/directional/get_U_inverse.h
    include/directional/get_U.h
    include/directional/FEM_suite.h
    include/directional/field_local_global_conversions.h
    include/directional/glyph_lines_raw.h
    include/directional/harmonic_basis.h
    include/directional/hodge_decomposition.h
    include/directional/index_prescription.h
    include/directional/line_boxes.h
    include/directional/line_cylinders.h
    include/directional/non_conforming_mesh.h
    include/directional/parameterize.h
    include/directional/point_spheres.h
    include/directional/polycurl_reduction.h
    include/directional/polyroots.h
    include/directional/polyvector_field.h
    include/directional/polyvector_to_raw.h
    include/directional/power_field.h
    include/directional/power_to_raw.h
    include/directional/power_to_representative.h
    include/directional/principal_matching.h
    include/directional/rawfield_to_columndirectional.h
    include/directional/read_matching.h
    include/directional/read_raw_field.h
    include/directional/read_singularities.h
    include/directional/representative_to_raw.h
    include/directional/rotation_to_raw.h
    include/directional/rotation_to_representative.h
    include/directional/seam_lines.h
    include/directional/setup_parameterization.h
    include/directional/singularity_spheres.h
    include/directional/streamlines.h
    include/directional/subdivide_field.h
    include/directional/tree.h
    include/directional/vertex_area_mesh.h
    include/directional/visualization_schemes.h
    include/directional/write_matching.h
    include/directional/write_raw_field.h
    include/directional/write_singularities.h
    include/directional/get_subdivision_matrix.h
)
set(SRC
    include/directional/angle_bound_frame_fields.cpp
    include/directional/conjugate_frame_fields.cpp
    include/directional/cut_mesh_with_singularities.cpp
    include/directional/field_local_global_conversions.cpp
    include/directional/polycurl_reduction.cpp
    include/directional/polyroots.cpp
    include/directional/streamlines.cpp
)

option(DIRECTIONAL_BUILD_STATIC "Build static Directional + libigl library" OFF)

if(DIRECTIONAL_BUILD_STATIC)
    add_library(Directional STATIC ${SRC} ${HEADERS})
    target_include_directories(Directional PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(Directional PUBLIC igl::core)

    
    #Particularly nice for Visual Studio users
    if(MSVC)
        set(_SRC_LIST ${SRC} ${HEADERS} ${PUBLIC_HEADERS})
        set(_TREE_LIST )
        foreach(EL IN LISTS _SRC_LIST)
            list(APPEND _TREE_LIST "${CMAKE_CURRENT_SOURCE_DIR}/${EL}")
        endforeach()
        source_group(${CMAKE_CURRENT_SOURCE_DIR}/include ${_TREE_LIST})	
    endif()
else()
    # Create an interface library that will just add the necessary include directory and potentially compiler flags/
    # preprocessor definitions
    add_library(Directional INTERFACE)

    target_sources(Directional INTERFACE ${SRC} ${HEADERS})

    if(MSVC)
        target_compile_definitions(Directional INTERFACE "_USE_MATH_DEFINES" )
    endif()
    target_include_directories(Directional INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_link_libraries(Directional INTERFACE igl::core)
endif()

# Allow user to specify whether or not to build tests and tutorial
option(DIRECTIONAL_BUILD_TESTS "Build Directional tests" OFF)
option(DIRECTIONAL_BUILD_TUTORIAL "Build Directional tutorial" OFF)

if(DIRECTIONAL_BUILD_TUTORIAL)
    add_subdirectory(tutorial)
endif()

# Add unit tests
if(DIRECTIONAL_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()