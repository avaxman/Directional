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
#include <directional/Subdivision_new/subdivide_directionals.h>

#include "tutorial_shared_path.h"
#include "directional/cut_mesh_with_singularities.h"

using namespace std;

Eigen::VectorXi matchingOrig, matchingCF, combedMatchingOrig, combedMatchingCF;
Eigen::VectorXd effortOrig, effortCF, combedEffortOrig, combedEffortCF;
Eigen::MatrixXi FMesh, FField, FSings, FSeams;
Eigen::MatrixXi EV, EF, FE;
// Fine level matrices
Eigen::MatrixXi EV_fine, EF_fine, FMesh_fine;
Eigen::VectorXi matchingCF_fine;
Eigen::MatrixXd VMesh_fine, rawfield_fine;

// Separate matrices for different visualizations
Eigen::MatrixXd VMesh, VField, VSings, VSeams, barycenters;
Eigen::MatrixXd CMesh, CField, CSings, CSeams;
Eigen::MatrixXd rawFieldOrig, rawFieldCF;
Eigen::MatrixXd combedFieldOrig, combedFieldCF;
Eigen::VectorXd curlOrig, curlCF; // norm of curl per edge
igl::opengl::glfw::Viewer viewer;

Eigen::VectorXi singVerticesOrig, singVerticesCF;
Eigen::VectorXi singIndicesOrig, singIndicesCF;


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


typedef enum
{
    ORIGINAL_FIELD, 
    ORIGINAL_CURL, 
    OPTIMIZED_FIELD, 
    OPTIMIZED_CURL,
    FINE_FIELD,
    COARSE_PARAM, 
    FINE_PARAM, 
    FINE_CURL
} ViewingModes;

enum DisplayMeshes
{
    Mesh = 0,
    Rawfield,
    Singularities,
    Seams,
    CoarseParam,
    FineParam,
    FineMesh,
    DisplayMeshCount
};

// Active view mode
ViewingModes viewingMode = ORIGINAL_FIELD;


void update_triangle_mesh()
{
    if ((viewingMode == ORIGINAL_FIELD) || (viewingMode == OPTIMIZED_FIELD) || viewingMode==FINE_FIELD) {
        viewer.data_list[0].set_colors(directional::default_mesh_color());
    }
    else if (viewingMode == COARSE_PARAM)
    {
        viewer.data_list[CoarseParam].show_faces = true;
    }
    else if(viewingMode == FINE_PARAM)
    {
        viewer.data_list[FineParam].show_faces = true;
    }
    else if(viewingMode == FINE_FIELD)
    {
        viewer.data_list[FineMesh].show_faces = true;
    }
    else {  //curl viewing - currently averaged to the face
        Eigen::VectorXd currCurl = AE2F * (viewingMode == ORIGINAL_CURL ? curlOrig : curlCF);
        igl::jet(currCurl, 0.0, curlMaxOrig, CMesh);
        viewer.data_list[0].set_colors(CMesh);
    }
}

void resetMeshData(igl::opengl::ViewerData& target, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& C, bool showFaces, bool showLines)
{
    target.clear();
    target.set_mesh(V, F);
    target.set_colors(C);
    target.show_faces = showFaces;
    target.show_lines = showLines;
}

void parameterize(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXi& EV, const Eigen::MatrixXi& EF, const Eigen::MatrixXi& FE,
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

void construct_fine_level_parameterization()
{
    directional::subdivision::subdivide_directionals(VMesh, FMesh, combedFieldCF, combedMatchingCF, targetLevel, VMesh_fine, FMesh_fine, rawfield_fine, matchingCF_fine, EF_fine, EV_fine);


    std::cout << "Max F val " << FMesh_fine.maxCoeff() << ", V count : " << VMesh_fine.rows() << std::endl;
    Eigen::MatrixXi FE_fine(FMesh_fine.rows(), 3);
    for(int e = 0; e < EV_fine.rows(); ++e)
    {
        int targetF = EF_fine(e, 0);
        for(int j = 0; j < 3; ++j)
        {
            const int v = FMesh_fine(targetF, j);
            if(EV_fine(e,0) != v && EV_fine(e,1) != v)
            {
                FE_fine(targetF, j) = e;
                break;
            }
        }
        targetF = EF_fine(e, 1);
        for (int j = 0; j < 3; ++j)
        {
            const int v = FMesh_fine(targetF, j);
            if (EV_fine(e, 0) != v && EV_fine(e, 1) != v)
            {
                FE_fine(targetF, j) = e;
                break;
            }
        }
    }
    // Expand singularities
    Eigen::VectorXi sings  = singIndicesCF;
    sings.conservativeResize(VMesh_fine.rows());
    for (int i = VMesh.rows(); i < sings.size(); ++i) sings(i) = 0;

    return;
    // Parameterize in the fine level and output to libigl viewer mesh
    parameterize(VMesh_fine, FMesh_fine, EV_fine, EF_fine, FE_fine, sings, rawfield_fine, matchingCF_fine, viewer.data_list[DisplayMeshes::FineParam]);
}


void update_raw_field_mesh()
{
    using namespace std;
    using namespace Eigen;

    if ((viewingMode == ORIGINAL_CURL) || (viewingMode == OPTIMIZED_CURL)) {
        for (int i = 1; i < DisplayMeshCount; i++) {  //hide all other meshes
            viewer.data_list[i].show_faces = false;
            viewer.data_list[i].show_lines = false;
        }
    }
    else if((viewingMode == COARSE_PARAM || viewingMode==FINE_PARAM))
    {
        viewer.data_list[Rawfield].show_faces = false;
    }
    else if(viewingMode == FINE_FIELD)
    {
        Eigen::MatrixXd fieldV,fieldC;
        Eigen::MatrixXi fieldF;
        directional::glyph_lines_raw(VMesh_fine, FMesh_fine, rawfield_fine,
            directional::indexed_glyph_colors(rawfield_fine), fieldV, fieldF, fieldC, 2.0);
        // Reset the mesh data for the field
        resetMeshData(viewer.data_list[1], fieldV, fieldF, fieldC, true, false);
    }
    else {
        directional::glyph_lines_raw(VMesh, FMesh, (viewingMode == ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF),
            directional::indexed_glyph_colors((viewingMode == ORIGINAL_FIELD ? combedFieldOrig : combedFieldCF)), VField, FField, CField, 2.0);
        // Reset the mesh data for the field
        resetMeshData(viewer.data_list[1], VField, FField, CField, true, false);

        //singularity mesh
        directional::singularity_spheres(VMesh, FMesh, N, (viewingMode == ORIGINAL_FIELD ? singVerticesOrig : singVerticesCF), (viewingMode == ORIGINAL_FIELD ? singIndicesOrig : singIndicesCF), VSings, FSings, CSings, 1.5);

        resetMeshData(viewer.data_list[2], VSings, FSings, CSings, true, false);

        //seam mesh
        directional::seam_lines(VMesh, FMesh, EV, (viewingMode == ORIGINAL_FIELD ? combedMatchingOrig : combedMatchingCF), VSeams, FSeams, CSeams, 2.0);

        resetMeshData(viewer.data_list[2], VSeams, FSeams, CSeams, true, false);
    }

}


bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
    // Set the view mode
    if ((key >= '1') && (key <= '7'))
    {
        viewingMode = static_cast<ViewingModes>(key - '1');
    }
        


    //do a batch of iterations
    if (key == 'A')
    {
        printf("--Improving Curl--\n");
        // Batch of 5 iterations
        for (int bi = 0; bi < 5; ++bi)
        {

            printf("\n\n **** Batch %d ****\n", iter);
            // Solve
            directional::polycurl_reduction_solve(pcrdata, params, rawFieldCF, iter == 0);
            ++iter;
            // Adjust the smoothness weight
            params.wSmooth *= params.redFactor_wsmooth;
        }

        // Reconstruct matching that minimizes the curl
        directional::curl_matching(VMesh, FMesh, EV, EF, FE, rawFieldCF, matchingCF, effortCF, curlCF);
        // Get indices from the effort for the matching
        directional::effort_to_indices(VMesh, FMesh, EV, EF, effortCF, matchingCF, N, singVerticesCF, singIndicesCF);
        // Comb the field
        directional::combing(VMesh, FMesh, EV, EF, FE, rawFieldCF, matchingCF, combedFieldCF);
        // Apply curl matching on the combed fields
        directional::curl_matching(VMesh, FMesh, EV, EF, FE, combedFieldCF, combedMatchingCF, combedEffortCF, curlCF);
        curlMax = curlCF.maxCoeff();
        std::cout << "curlMax optimized: " << curlMax << std::endl;
    }
    else if(key == 'C')
    {
        parameterize(VMesh, FMesh, EV, EF, FE, singVerticesCF, combedFieldCF, combedMatchingCF, viewer.data_list[DisplayMeshes::CoarseParam]);
    }
    if(key == 'S')
    {
        std::cout << "Subdividing and parameterizing " << std::endl;
        construct_fine_level_parameterization();
    }

    // Write the current field to file.
    if (key == 'W') {
        if (directional::write_raw_field(TUTORIAL_SHARED_PATH "/cheburashka-cf.rawfield", rawFieldCF))
            std::cout << "Saved raw field" << std::endl;
        else
            std::cout << "Unable to save raw field. " << std::endl;
    }

    update_triangle_mesh();
    update_raw_field_mesh();
    return false;
}

int main(int argc, char *argv[])
{

    std::cout <<
        "  A      Optimize 5 batches for curl reduction." << std::endl <<
        "  1      Original field" << std::endl <<
        "  2      L2 norm of original-field curl" << std::endl <<
        "  3      Curl-reduced field" << std::endl <<
        "  4      Curl of curl-reduced field." << std::endl;

    // Load a mesh
    igl::readOFF(TUTORIAL_SHARED_PATH "/cheburashka-subdivision.off", VMesh, FMesh);
    directional::read_raw_field(TUTORIAL_SHARED_PATH "/cheburashka-subdivision.rawfield", N, rawFieldOrig);
    igl::edge_topology(VMesh, FMesh, EV, FE, EF);

    igl::barycenter(VMesh, FMesh, barycenters);

    Eigen::VectorXi prinIndices;
    directional::curl_matching(VMesh, FMesh, EV, EF, FE, rawFieldOrig, matchingOrig, effortOrig, curlOrig);
    directional::effort_to_indices(VMesh, FMesh, EV, EF, effortOrig, matchingOrig, N, singVerticesOrig, singIndicesOrig);
    curlMaxOrig = curlOrig.maxCoeff();
    curlMax = curlMaxOrig;
    std::cout << "curlMax original: " << curlMax << std::endl;

    directional::combing(VMesh, FMesh, EV, EF, FE, rawFieldOrig, matchingOrig, combedFieldOrig);
    directional::curl_matching(VMesh, FMesh, EV, EF, FE, combedFieldOrig, combedMatchingOrig, combedEffortOrig, curlOrig);

    //trivial constraints
    Eigen::VectorXi b; b.resize(1); b << 0;
    Eigen::MatrixXd bc; bc.resize(1, 6); bc << rawFieldOrig.row(0).head(6);
    Eigen::VectorXi blevel; blevel.resize(1); b << 1;
    directional::polycurl_reduction_precompute(VMesh, FMesh, b, bc, blevel, rawFieldOrig, pcrdata);

    rawFieldCF = rawFieldOrig;
    matchingCF = matchingOrig;
    effortCF = effortOrig;
    combedFieldCF = combedFieldOrig;
    combedMatchingCF = combedMatchingOrig;
    combedEffortCF = combedEffortOrig;
    curlCF = curlOrig;
    singVerticesCF = singVerticesOrig;
    singIndicesCF = singIndicesOrig;


    //triangle mesh setup
    viewer.data().set_mesh(VMesh, FMesh);
    viewer.data().set_colors(directional::default_mesh_color());
    viewer.data().show_lines = false;
    viewer.data().set_face_based(true);

    //apending and updating raw field mesh
    viewer.append_mesh();

    //singularity mesh
    viewer.append_mesh();

    //seam mesh
    viewer.append_mesh();

    // Coarse param mesh
    viewer.append_mesh();

    // Fine param mesh
    viewer.append_mesh();

    // FIne mesh
    viewer.append_mesh();

    update_triangle_mesh();
    update_raw_field_mesh();
    viewer.selected_data_index = 0;

    //creating the AE2F operator
    std::vector<Eigen::Triplet<double> > AE2FTriplets;
    for (int i = 0; i < EF.rows(); i++) {
        AE2FTriplets.push_back(Eigen::Triplet<double>(EF(i, 0), i, 1.0));
        AE2FTriplets.push_back(Eigen::Triplet<double>(EF(i, 1), i, 1.0));
    }
    AE2F.resize(FMesh.rows(), EF.rows());
    AE2F.setFromTriplets(AE2FTriplets.begin(), AE2FTriplets.end());

    viewer.callback_key_down = &key_down;
    viewer.launch();

    return 0;
}
