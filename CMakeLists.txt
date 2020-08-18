cmake_minimum_required(VERSION 3.1)
project(Directional)

### conditionally compile certain modules depending on libraries found on the system
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(TEST_SHARED_PATH ${CMAKE_CURRENT_SOURCE_DIR}/tutorials/shared CACHE PATH "location of resources for testing")

### Adding libIGL and Directional: choose the path to your local copy
find_package(LIBIGL REQUIRED QUIET)

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake) # add FindMatlab module
option(DIRECTIONAL_BUILD_MEX "Build Directional MEX files" ON)
option(DIRECTIONAL_BUILD_TUTORIALS "Build Directional tutorials" ON)
option(DIRECTIONAL_BUILD_TESTS "Build Directional tests" ON)
option(DIRECTIONAL_BUILD_STATIC_LIBRARY "Build Directional as static library" ON)

set(HEADERS
include/directional/Deprecated/cycle_curvature.h
include/directional/Deprecated/draw_cycles.h
include/directional/Deprecated/representative_to_rotation.h
include/directional/Subdivision/DCEL.h
include/directional/Subdivision/directional_subdivision.h
include/directional/Subdivision/DirectionalSubdivider.h
include/directional/Subdivision/EdgeData.h
include/directional/Subdivision/EdgeFieldSubdivider.h
include/directional/Subdivision_new/build_directional_subdivision_operators.h
include/directional/Subdivision/FaceFieldSubdivider.h
include/directional/Subdivision/HalfBoxSplineCoefficientProvider.h
include/directional/glyph_arrows_raw.h
include/directional/Subdivision/HalfCurlCoefficientProvider.h
include/directional/Subdivision/LoopCoefficientProvider.h
include/directional/Subdivision_new/Sc_directional_triplet_provider.h
include/directional/Subdivision_new/Se_directional_triplet_provider.h
include/directional/Subdivision_new/gamma_subdivision_matrix.h
include/directional/Subdivision_new/subdivide_field.h
include/directional/Subdivision_new/iterate_branched_rings.h
#include/directional/surface_lic.h
include/directional/Subdivision/OneFormCoefficientProvider.h
include/directional/Subdivision/OneFormSubdivider.h
include/directional/Subdivision/RangeHelpers.h
include/directional/Subdivision/SimpleTimer.h
include/directional/Subdivision/subdivision.h
include/directional/Subdivision/SubdivisionBuilder.h
include/directional/Subdivision/VertexFieldSubdivider.h
include/directional/Subdivision_new/build_shm_subdivision.h
include/directional/Subdivision_new/build_subdivision_operators.h
include/directional/Subdivision_new/hbspline_coefficients.h
include/directional/Subdivision_new/is_edgedata_consistent.h
include/directional/Subdivision_new/iterate_rings.h
include/directional/Subdivision_new/loop_coefficients.h
include/directional/Subdivision_new/quadrisect.h
include/directional/Subdivision_new/Sc_triplet_provider.h
include/directional/Subdivision_new/Se_triplet_provider.h
include/directional/Subdivision_new/Sf_triplet_provider.h
include/directional/Subdivision_new/shm_edge_topology.h
include/directional/Subdivision_new/shm_halfcurl_coefficients.h
include/directional/Subdivision_new/shm_oneform_coefficients.h
include/directional/Subdivision_new/Sv_triplet_provider.h
include/directional/angle_bound_frame_fields.cpp
include/directional/angle_bound_frame_fields.h
include/directional/block_diag.h
include/directional/combing.h
include/directional/conjugate_frame_fields.cpp
include/directional/conjugate_frame_fields.h
include/directional/ConjugateFFSolverData.h
include/directional/curl_matching.h
include/directional/cut_mesh_with_singularities.cpp
include/directional/cut_mesh_with_singularities.h
include/directional/dcel.h
#include/directional/dijkstra.cpp
#include/directional/dijkstra.h
include/directional/dir_assert.h
include/directional/DirectionalGamma_Suite.h
include/directional/dual_cycles.h
include/directional/dynamic_visualization.h
include/directional/edge_diamond_mesh.h
include/directional/effort_to_indices.h
include/directional/FEM_masses.h
include/directional/FEM_suite.h
include/directional/field_local_global_conversions.cpp
include/directional/field_local_global_conversions.h
include/directional/Gamma_suite.h
include/directional/glyph_lines_raw.h
include/directional/harmonic_basis.h
include/directional/hodge_decomposition.h
include/directional/index_prescription.h
include/directional/iterate_rings.h
include/directional/line_boxes.h
include/directional/line_cylinders.h
include/directional/non_conforming_mesh.h
include/directional/parameterize.h
include/directional/point_spheres.h
include/directional/polycurl_reduction.cpp
include/directional/polycurl_reduction.h
include/directional/polyroots.cpp
include/directional/polyroots.h
include/directional/polyvector_field.h
include/directional/polyvector_to_raw.h
include/directional/power_field.h
include/directional/power_to_raw.h
include/directional/power_to_representative.h
include/directional/principal_matching.h
include/directional/read_matching.h
include/directional/read_raw_field.h
include/directional/read_singularities.h
include/directional/representative_to_raw.h
include/directional/rotation_to_raw.h
include/directional/rotation_to_representative.h
include/directional/seam_lines.h
include/directional/setup_parameterization.h
include/directional/singularity_spheres.h
include/directional/streamlines.cpp
include/directional/Subdivision_new/subdivide_directionals.h
include/directional/streamlines.h
include/directional/tree.h
include/directional/vertex_area_mesh.h
include/directional/visualization_schemes.h
include/directional/write_raw_field.h
#include/directional/directional_types.h
include/directional/write_singularities.h
)
source_group(TREE ${CMAKE_CURRENT_SOURCE_DIR}/include FILES ${HEADERS})

# Interface for header only?
if(${DIRECTIONAL_BUILD_STATIC_LIBRARY})
    add_library(Directional STATIC ${HEADERS})
else()
    add_library(Directional INTERFACE)
    target_include_directories(Directional INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/include)
    target_sources(Directional  ${HEADERS})
endif()
target_link_libraries(Directional PUBLIC igl::core)
target_include_directories(Directional PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Conditionally add parts of Directional for compilation
if(${DIRECTIONAL_BUILD_TUTORIALS})
    add_subdirectory(tutorial)
endif()
if(${DIRECTIONAL_BUILD_TESTS})
    add_subdirectory(tests)
endif()
if(${DIRECTIONAL_BUILD_MEX})
    add_subdirectory(matlab)
endif()
# Build Python?