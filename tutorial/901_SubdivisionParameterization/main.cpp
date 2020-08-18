#include <iostream>
#include <fstream>
#include <igl/readOFF.h>
#include <igl/readOBJ.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/local_basis.h>
#include <igl/avg_edge_length.h>
#include <igl/is_border_vertex.h>
#include <igl/adjacency_list.h>
#include <igl/vertex_triangle_adjacency.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_topology.h>
#include <igl/jet.h>
#include <igl/barycenter.h>
#include <igl/false_barycentric_subdivision.h>
#include <directional/visualization_schemes.h>
#include <directional/glyph_lines_raw.h>
#include <directional/seam_lines.h>
#include <directional/read_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/singularity_spheres.h>
#include <igl/edge_flaps.h>
#include <directional/combing.h>
#include <directional/polycurl_reduction.h>
#include <directional/write_raw_field.h>
#include <directional/setup_parameterization.h>
#include <directional/parameterize.h>
#include <directional/subdivide_directionals.h>

#include "tutorial_shared_path.h"
#include "directional/cut_mesh_with_singularities.h"
#include <unordered_set>

using namespace std;

Eigen::VectorXi matchingOrig, matchingCF, combedMatchingOrig, combedMatchingCF;
Eigen::VectorXd effortOrig, effortCF, combedEffortOrig, combedEffortCF;
Eigen::MatrixXi FMesh_coarse, FField, FSings, FSeams;
Eigen::MatrixXi EV_coarse, EF_coarse, FE_coarse;
// Fine level matrices
Eigen::MatrixXi EV_fine, EF_fine, FMesh_fine;
Eigen::VectorXi matchingCF_fine;
Eigen::MatrixXd VMesh_fine, rawfield_fine;

// Separate matrices for different visualizations
Eigen::MatrixXd VMesh_coarse, VField, VSings, VSeams, barycenters;
Eigen::MatrixXd CMesh, CField, CSings, CSeams;
Eigen::MatrixXd rawFieldOrig, rawFieldCF;
Eigen::MatrixXd combedFieldOrig, combedFieldCF;
Eigen::VectorXd curlOrig, curlCF; // norm of curl per edge

// The igl viewer
igl::opengl::glfw::Viewer viewer;

Eigen::VectorXi singVerticesOrig, singVerticesCF;
Eigen::VectorXi singIndicesOrig, singIndicesCF;

// Fine directional curl operator
Eigen::SparseMatrix<double> C_directional_fine;
Eigen::VectorXd l2Curl_fine;


//for averaging curl to faces, for visualization
Eigen::SparseMatrix<double> AE2F;

double curlMax, curlMaxOrig;
int N;
int targetLevel = 1;
int iter = 0;

// The set of parameters for calculating the curl-free fields
directional::polycurl_reduction_parameters params;

// Solver data (needed for precomputation)
directional::PolyCurlReductionSolverData pcrdata;

//Viewing modes
typedef enum
{
    ORIGINAL_FIELD, 
    ORIGINAL_CURL, 
    OPTIMIZED_FIELD, 
    OPTIMIZED_CURL,
    FINE_FIELD,
    FINE_CURL,
    COARSE_PARAM, 
    FINE_PARAM 
} ViewingModes;

// The visualization meshes
enum DisplayMeshes
{
    Mesh = 0,
    Rawfield,
    Singularities,
    Seams,
    FineMesh,
    FineField,
    CoarseParam,
    FineParam,
    DisplayMeshCount
};

// Meshes to show for each view
std::vector<unordered_set<DisplayMeshes>> meshesPerView = {
    {Mesh, Rawfield, Singularities, Seams},
    {Mesh},
    {Mesh, Rawfield, Singularities, Seams},
    {Mesh},
    {FineMesh, FineField},
    {FineMesh},
    {CoarseParam},
    {FineParam}
};

// Active view mode
ViewingModes viewingMode = ORIGINAL_FIELD;

// Helper for resetting viewer data on the Viewer object.
void resetMeshData(igl::opengl::ViewerData& target, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& C, bool showFaces, bool showLines)
{
    target.clear();
    target.set_mesh(V, F);
    target.set_colors(C);
    target.show_faces = showFaces;
    target.show_lines = showLines;
}

// Parameterize a curlfree optimized field
void parameterize_cf(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXi& EV, const Eigen::MatrixXi& EF, const Eigen::MatrixXi& FE,
    const Eigen::VectorXi& singVertices, const Eigen::MatrixXd& rawField, const Eigen::VectorXi& matching, igl::opengl::ViewerData& output)
{
    Eigen::MatrixXd combedField, VMeshCut;
    Eigen::MatrixXi FMeshCut;
    Eigen::VectorXi combedMatching;
    directional::ParameterizationData pd;
    directional::cut_mesh_with_singularities(V, F, singVertices, pd.face2cut);
    directional::combing(V, F, EV, EF, FE, pd.face2cut, rawField, matching, combedField, combedMatching);
    //directional::curl_matching()

    std::cout << "Setting up parameterization" << std::endl;

    Eigen::MatrixXi symmFunc(4, 2);
    symmFunc << 1, 0,
        0, 1,
        -1, 0,
        0, -1;

    directional::setup_parameterization(symmFunc, V, F, EV, EF, FE, combedMatching, singVertices, pd, VMeshCut, FMeshCut);

    double lengthRatio = 0.01;
    bool isInteger = false;  //do not do translational seamless.
    std::cout << "Solving parameterization" << std::endl;
    Eigen::MatrixXd cutReducedUV, cutFullUV, cornerWholeUV;
    directional::parameterize(V, F, FE, combedField, lengthRatio, pd, VMeshCut, FMeshCut, isInteger, cutReducedUV, cutFullUV, cornerWholeUV);
    cutFullUV = cutFullUV.block(0, 0, cutFullUV.rows(), 2);
    std::cout << "Done!" << std::endl;

    output.clear();
    output.set_mesh(VMeshCut, FMeshCut);
    output.set_uv(cutFullUV);
    output.show_texture = true;
    output.show_lines = false;
}

// Compute curl statistics on directional
void computeDirectionalCurl(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXi& EV, const Eigen::MatrixXi& EF, const Eigen::MatrixXd& rawField, const Eigen::VectorXi& matching, 
    Eigen::VectorXd& curlColumn, double& maxSqNorm, Eigen::VectorXd& l2CurlOnFaces)
{
    Eigen::MatrixXi EI_local, SFE_local;
    directional::shm_edge_topology(F, EV, EF, EI_local, SFE_local);
    Eigen::SparseMatrix<double> C, columndirectional_to_g2, avgToFaces;
    directional::Matched_Curl(EF, SFE_local, EI_local, matching, F.rows(), N, C);
    directional::columndirectional_to_gamma2_matrix(V, F, EV, SFE_local, EF, N, columndirectional_to_g2);
    Eigen::VectorXd columnDirectionalFine;
    directional::rawfield_to_columndirectional(rawField, N, columnDirectionalFine);
    // Output
    curlColumn = C * columndirectional_to_g2 * columnDirectionalFine;

    // Compute l2 norm
    Eigen::MatrixXd curlPerEdge;
    directional::columnfunction_to_rawfunction(curlColumn, N, 1, curlPerEdge);
    auto sqNorm = curlPerEdge.rowwise().norm();
    maxSqNorm = sqNorm.maxCoeff();
    // Copy to faces
    directional::Edge_To_Face_Average(SFE_local, EF.rows(), avgToFaces);
    l2CurlOnFaces = avgToFaces * sqNorm;
}

// Compute subdivision of coarse optimized field and apply parameterization
void construct_fine_level_parameterization()
{
    {
        Eigen::VectorXd curlColumn, l2CurlOnFaces;
        double maxSqNorm;
        computeDirectionalCurl(VMesh_coarse, FMesh_coarse, EV_coarse, EF_coarse, combedFieldCF, combedMatchingCF, curlColumn, maxSqNorm, l2CurlOnFaces);
        
        std::cout << "Coarse max l2 curl before: " << maxSqNorm << std::endl;
    }
    directional::subdivide_directionals(VMesh_coarse, FMesh_coarse, EV_coarse, EF_coarse, combedFieldCF, combedMatchingCF, targetLevel, VMesh_fine, FMesh_fine,
        EV_fine, EF_fine, rawfield_fine, matchingCF_fine);

    {
        Eigen::VectorXd curlColumn;
        double maxSqNorm;
        computeDirectionalCurl(VMesh_fine, FMesh_fine, EV_fine, EF_fine, rawfield_fine, matchingCF_fine, curlColumn, maxSqNorm, l2Curl_fine);

        std::cout << "Fine max l2 curl: " << maxSqNorm << std::endl;
    }
    
    // Set the fine mesh
    resetMeshData(viewer.data_list[DisplayMeshes::FineMesh], VMesh_fine, FMesh_fine, directional::default_mesh_color(), false, false);

    // Reconstruct FE where FE(f,0) is now the edge that has vertices F(f,0) and F(f,1) (igl::edge_topology style...)
    Eigen::MatrixXi FE_fine(FMesh_fine.rows(), 3);
    for (int e = 0; e < EV_fine.rows(); ++e)
    {
        for (int s = 0; s < 2; ++s)
        {
            int targetF = EF_fine(e, s);
            for (int j = 0; j < 3; ++j)
            {
                const int otherV = FMesh_fine(targetF, (j + 2) % 3);
                if ( EV_fine(e, 0) != otherV && EV_fine(e, 1) != otherV)
                {
                    FE_fine(targetF, j) = e;
                    break;
                }
            }
        }
    }
    // Parameterize in the fine level and output to libigl viewer mesh
    // We can reuse singVerticesCF since the coarse level vertex IDs remained unchanged and we assume the singularities to be stationary under the subdivision
    parameterize_cf(VMesh_fine, FMesh_fine, EV_fine, EF_fine, FE_fine, singVerticesCF, rawfield_fine, matchingCF_fine, viewer.data_list[DisplayMeshes::FineParam]);
}

/*
 * Update meshes for view and make the appropriate meshes visible
 */
void update_view()
{
    switch(viewingMode)
    {
    case ORIGINAL_FIELD:
        {
        directional::glyph_lines_raw(VMesh_coarse, FMesh_coarse, (viewingMode == ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF),
            directional::indexed_glyph_colors((viewingMode == ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF)), VField, FField, CField, 2.0);
        // Reset the mesh data for the field
        resetMeshData(viewer.data_list[1], VField, FField, CField, true, false);

        //singularity mesh
        directional::singularity_spheres(VMesh_coarse, FMesh_coarse, N, (viewingMode == ORIGINAL_FIELD ? singVerticesOrig : singVerticesCF), (viewingMode == ORIGINAL_FIELD ? singIndicesOrig : singIndicesCF), VSings, FSings, CSings, 1.5);

          resetMeshData(viewer.data_list[DisplayMeshes::Singularities], VSings, FSings, CSings, true, false);

        //seam mesh
        directional::seam_lines(VMesh_coarse, FMesh_coarse, EV_coarse, (viewingMode == ORIGINAL_FIELD ? combedMatchingOrig : combedMatchingCF), VSeams, FSeams, CSeams, 2.0);

        resetMeshData(viewer.data_list[DisplayMeshes::Seams], VSeams, FSeams, CSeams, true, false);
        viewer.data_list[Mesh].set_colors(directional::default_mesh_color());
        break;
        }
    case ORIGINAL_CURL:
        {
        Eigen::VectorXd currCurl = AE2F * curlOrig;
        igl::jet(currCurl, 0.0, curlMaxOrig, CMesh);
        viewer.data_list[DisplayMeshes::Mesh].set_colors(CMesh);
        }
        break;
    case OPTIMIZED_FIELD:
    {
        directional::glyph_lines_raw(VMesh_coarse, FMesh_coarse, (viewingMode == ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF),
            directional::indexed_glyph_colors((viewingMode == ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF)), VField, FField, CField, 2.0);
        // Reset the mesh data for the field
        resetMeshData(viewer.data_list[1], VField, FField, CField, true, false);

        //singularity mesh
        directional::singularity_spheres(VMesh_coarse, FMesh_coarse, N, (viewingMode == ORIGINAL_FIELD ? singVerticesOrig : singVerticesCF), (viewingMode == ORIGINAL_FIELD ? singIndicesOrig : singIndicesCF), VSings, FSings, CSings, 1.5);

        resetMeshData(viewer.data_list[DisplayMeshes::Singularities], VSings, FSings, CSings, true, false);

        //seam mesh
        directional::seam_lines(VMesh_coarse, FMesh_coarse, EV_coarse, (viewingMode == ORIGINAL_FIELD ? combedMatchingOrig : combedMatchingCF), VSeams, FSeams, CSeams, 2.0);

        resetMeshData(viewer.data_list[DisplayMeshes::Seams], VSeams, FSeams, CSeams, true, false);
        viewer.data_list[Mesh].set_colors(directional::default_mesh_color());
        break;
    }
    case OPTIMIZED_CURL:
        {
        Eigen::VectorXd currCurl = AE2F * curlCF;
        igl::jet(currCurl, 0.0, curlMaxOrig, CMesh);
        viewer.data_list[DisplayMeshes::Mesh].set_colors(CMesh);
        }
        break;
    case FINE_FIELD:
        {
        directional::glyph_lines_raw(VMesh_fine, FMesh_fine, rawfield_fine,
            directional::indexed_glyph_colors(rawfield_fine), VField, FField, CField, 2.0);
        // Reset the mesh data for the field
        resetMeshData(viewer.data_list[FineField], VField, FField, CField, true, false);
        viewer.data_list[FineMesh].set_colors(directional::default_mesh_color());
        }
        break;
    case FINE_CURL:
        {
        Eigen::MatrixXd Cs;
        igl::jet(l2Curl_fine, 0.0, curlMaxOrig * 0.25, Cs);
        viewer.data_list[DisplayMeshes::FineMesh].set_colors(Cs);
        break;
        }
    case COARSE_PARAM:
        {
        parameterize_cf(VMesh_coarse, FMesh_coarse, EV_coarse, EF_coarse, FE_coarse, singVerticesCF, combedFieldCF, combedMatchingCF, viewer.data_list[CoarseParam]);
        break;
        }
    case FINE_PARAM:
        viewer.data_list[FineParam].show_faces = true;
        break;
        /*
            COARSE_PARAM,
            FINE_PARAM*/
    default:
        break;
    }
    for(DisplayMeshes m = Mesh; m < DisplayMeshCount; m = static_cast<DisplayMeshes>(m+1))
    {
        bool show = meshesPerView[viewingMode].find(m) != meshesPerView[viewingMode].end();
        viewer.data_list[m].show_faces = show;
    }
    
}

bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
    // Set the view mode
    if ((key >= '1') && (key <= '8'))
    {
        viewingMode = static_cast<ViewingModes>(key - '1');
    }
        


    //do a batch of iterations
    if (key == 'A')
    {
        printf("--Improving Curl--\n");
        // Batch of 5 iterations
        for (int bi = 0; bi < 60; ++bi)
        {

            printf("\n\n **** Batch %d ****\n", iter);
            // Solve
            directional::polycurl_reduction_solve(pcrdata, params, rawFieldCF, iter == 0);
            ++iter;
            // Adjust the smoothness weight
            params.wSmooth *= params.redFactor_wsmooth;
        }

        // Reconstruct matching that minimizes the curl
        directional::curl_matching(VMesh_coarse, FMesh_coarse, EV_coarse, EF_coarse, FE_coarse, rawFieldCF, matchingCF, effortCF, curlCF);
        // Get indices from the effort for the matching
        directional::effort_to_indices(VMesh_coarse, FMesh_coarse, EV_coarse, EF_coarse, effortCF, matchingCF, N, singVerticesCF, singIndicesCF);
        // Comb the field
        directional::combing(VMesh_coarse, FMesh_coarse, EV_coarse, EF_coarse, FE_coarse, rawFieldCF, matchingCF, combedFieldCF);
        // Apply curl matching on the combed fields
        directional::curl_matching(VMesh_coarse, FMesh_coarse, EV_coarse, EF_coarse, FE_coarse, combedFieldCF, combedMatchingCF, combedEffortCF, curlCF);
        curlMax = curlCF.maxCoeff();
        std::cout << "curlMax optimized: " << curlMax << std::endl;
    }
    if(key == 'S')
    {
        std::cout << "Subdividing and parameterizing " << std::endl;
        construct_fine_level_parameterization();
    }

    // Write the current coarse field to file.
    if (key == 'W') {
        if (directional::write_raw_field(TUTORIAL_SHARED_PATH "/cheburashka-subdivision-cf.rawfield", rawFieldCF))
            std::cout << "Saved raw field" << std::endl;
        else
            std::cout << "Unable to save raw field. " << std::endl;
    }

    //update_triangle_mesh();
    update_view();
    return false;
}

int main(int argc, char *argv[])
{
    auto keyAction = [](const std::string& key, const std::string& description)
    {
        std::cout << "  " << key << "      " << description << std::endl;
    };
    keyAction("A", "Optimize 60 batches for curl reduction.");
    keyAction("1", "Show original field.");
    keyAction("2", "Show L2 norm of original-field curl.");
    keyAction("3", "Show Curl-reduced field.");
    keyAction("4", "Show Curl of curl-reduced field.");
    keyAction("S", "Subdivide and parameterize.");
    keyAction("5", "Show subdivided curl-reduced field.");
    keyAction("6", "Show L2 norm of subdivided curl-reduced field curl.");
    keyAction("7", "Show coarse parameterization of curl-reduced field.");
    keyAction("8", "Show fine parameterization of subdivided curl-reduced field.");

    // Load a mesh
    igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka-subdivision.off", VMesh_coarse, FMesh_coarse);
    directional::read_raw_field(TUTORIAL_SHARED_PATH "/cheburashka-subdivision.rawfield", N, rawFieldOrig);
    igl::edge_topology(VMesh_coarse, FMesh_coarse, EV_coarse, FE_coarse, EF_coarse);

    // compute barycenters
    igl::barycenter(VMesh_coarse, FMesh_coarse, barycenters);

    // Setup initial mathcing and determine singularities
    Eigen::VectorXi prinIndices;
    directional::curl_matching(VMesh_coarse, FMesh_coarse, EV_coarse, EF_coarse, FE_coarse, rawFieldOrig, matchingOrig, effortOrig, curlOrig);
    directional::effort_to_indices(VMesh_coarse, FMesh_coarse, EV_coarse, EF_coarse, effortOrig, matchingOrig, N, singVerticesOrig, singIndicesOrig);
    curlMaxOrig = curlOrig.maxCoeff();
    curlMax = curlMaxOrig;
    std::cout << "curlMax original: " << curlMax << std::endl;

    directional::combing(VMesh_coarse, FMesh_coarse, EV_coarse, EF_coarse, FE_coarse, rawFieldOrig, matchingOrig, combedFieldOrig);
    directional::curl_matching(VMesh_coarse, FMesh_coarse, EV_coarse, EF_coarse, FE_coarse, combedFieldOrig, combedMatchingOrig, combedEffortOrig, curlOrig);

    //trivial constraints
    Eigen::VectorXi b; b.resize(1); b << 0;
    Eigen::MatrixXd bc; bc.resize(1, 6); bc << rawFieldOrig.row(0).head(6);
    Eigen::VectorXi blevel; blevel.resize(1); b << 1;
    // Precompute polycurl reduction data.
    directional::polycurl_reduction_precompute(VMesh_coarse, FMesh_coarse, b, bc, blevel, rawFieldOrig, pcrdata);

    rawFieldCF = rawFieldOrig;
    matchingCF = matchingOrig;
    effortCF = effortOrig;
    combedFieldCF = combedFieldOrig;
    combedMatchingCF = combedMatchingOrig;
    combedEffortCF = combedEffortOrig;
    curlCF = curlOrig;
    singVerticesCF = singVerticesOrig;
    singIndicesCF = singIndicesOrig;


    //Setup coarse mesh
    viewer.data().set_mesh(VMesh_coarse, FMesh_coarse);
    viewer.data().set_colors(directional::default_mesh_color());
    viewer.data().show_lines = false;
    viewer.data().set_face_based(true);

    // Add a mesh for every display mesh
    for(int i = 1; i < DisplayMeshes::DisplayMeshCount; ++i)
    {
        viewer.append_mesh();
    }

    // Update view
    update_view();

    viewer.selected_data_index = 0;

    //creating the AE2F operator
    std::vector<Eigen::Triplet<double> > AE2FTriplets;
    for (int i = 0; i < EF_coarse.rows(); i++) {
        AE2FTriplets.push_back(Eigen::Triplet<double>(EF_coarse(i, 0), i, 1.0));
        AE2FTriplets.push_back(Eigen::Triplet<double>(EF_coarse(i, 1), i, 1.0));
    }
    AE2F.resize(FMesh_coarse.rows(), EF_coarse.rows());
    AE2F.setFromTriplets(AE2FTriplets.begin(), AE2FTriplets.end());

    viewer.callback_key_down = &key_down;
    viewer.launch();

    return 0;
}
