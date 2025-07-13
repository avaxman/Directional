#include <iostream>
#include <Eigen/Core>
#include <directional/readOBJ.h>
#include <directional/write_raw_field.h>
#include <directional/polyvector_to_raw.h>
#include <directional/raw_to_polyvector.h>
#include <directional/polyvector_iteration_functions.h>
#include <directional/PolyVectorData.h>
#include <directional/polyvector_field.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>
#include <directional/IntrinsicVertexTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/curl_matching.h>
#include <directional/extrinsic_intrinsic_matrices.h>

Eigen::VectorXi constVertices;
directional::TriMesh mesh;
directional::IntrinsicVertexTangentBundle vtb;
directional::CartesianField pvFieldConjugate, rawFieldConjugate,constraintsField, rawFieldOrig, pvFieldOrig;
Eigen::MatrixXd constVectors, constSources;
Eigen::VectorXd alignWeights;
double smoothWeight, roSyWeight, globalAlignWeight;
directional::DirectionalViewer viewer;
directional::PolyVectorData pvData;
std::vector<directional::PvIterationFunction> iterationFunctions;   //The iteration functions used in sequence
int N = 4;

void callbackFunc() {
    ImGui::PushItemWidth(100);
    if (ImGui::Button("Perform Iteration")) {
        directional::polyvector_iterate(pvData, pvFieldConjugate, iterationFunctions, 1);
        directional::polyvector_to_raw(pvFieldConjugate, rawFieldConjugate, N%2==0);
        directional::principal_matching(rawFieldConjugate);
        viewer.set_cartesian_field(rawFieldConjugate,"Conjugate Field", 0);
        viewer.set_field_color({107.0/255.0, 8.0/255.0, 125.0}, 0);
    }
    if (ImGui::Button("Perform Implicit Step")) {
        directional::polyvector_iterate(pvData, pvFieldConjugate, std::vector<directional::PvIterationFunction>(), 1);
        directional::polyvector_to_raw(pvFieldConjugate, rawFieldConjugate, N%2==0);
        directional::principal_matching(rawFieldConjugate);
        viewer.set_cartesian_field(rawFieldConjugate,"Conjugate Field", 0);
        viewer.set_field_color({107.0/255.0, 8.0/255.0, 125.0}, 0);
    }
    ImGui::PopItemWidth();
}


int main()
{
    // Load mesh
    directional::readOBJ(TUTORIAL_DATA_PATH "/inspired_mesh.obj", mesh);
    vtb.init(mesh);
    pvFieldConjugate.init(vtb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
    pvFieldOrig.init(vtb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
    
    //Soft constraints on principal directions
    constVertices.resize(2*mesh.V.rows());
    constVectors.resize(2*mesh.V.rows(),3);
    constSources.resize(2*mesh.V.rows(),3);
    alignWeights.resize(2*mesh.V.rows());
    Eigen::VectorXd confidence = mesh.GaussianCurvature.cwiseAbs();
    for (int i=0;i<mesh.V.rows();i++)
        if (mesh.isBoundaryVertex(i))
            confidence(i) = 0.0;  //not aligning to these vertices
    
    confidence = (confidence.array()-confidence.cwiseAbs().minCoeff())/(confidence.maxCoeff()-confidence.minCoeff());
    
    for (int i=0;i<mesh.V.rows();i++){
        constVertices(2*i) = i;
        constSources.row(2*i) = mesh.V.row(i);
        constVectors.row(2*i) = mesh.minVertexPrincipalDirections.row(i);
        alignWeights(2*i) = confidence(i);
        
        constSources.row(2*i+1) = mesh.V.row(i);
        constVertices(2*i+1) = i;
        constVectors.row(2*i+1) = mesh.maxVertexPrincipalDirections.row(i);
        alignWeights(2*i+1) = confidence(i);
    }
    
    smoothWeight = 1.0;
    roSyWeight = 1.0;
    globalAlignWeight = 0.01;
    //TODO: probably need a soft conjugate
    pvData.N = N;
    pvData.tb = &vtb;
    pvData.verbose = true;
    pvData.constSpaces = constVertices;
    pvData.constVectors = constVectors;
    pvData.wAlignment = globalAlignWeight*alignWeights;
    pvData.wSmooth = smoothWeight;
    pvData.wRoSy = roSyWeight;
    
    //Computing regular PolyuVector field without iterations
    directional::polyvector_field(pvData, pvFieldOrig);
    directional::polyvector_to_raw(pvFieldOrig, rawFieldOrig, N%2==0);
    directional::principal_matching(rawFieldOrig);
    
    //Iterating for a conjugate field
    pvData.iterationMode = true;
    pvData.initImplicitFactor = 0.5;
    pvData.implicitScheduler = 0.9;
    iterationFunctions.push_back(directional::conjugate);
    iterationFunctions.push_back(directional::hard_normalization);
    directional::polyvector_field(pvData, pvFieldConjugate);
    rawFieldConjugate = rawFieldOrig;
    
    //Visualization
    viewer.init();
    viewer.set_surface_mesh(mesh);
    viewer.highlight_vertices(constVertices);
    //viewer.set_raw_field(constSources, constVectors, "Constraints",  0, 0.3*mesh.avgEdgeLength);
    //viewer.set_field_color(directional::DirectionalViewer::default_vector_constraint_color());
    //viewer.set_cartesian_field(rawFieldOrig,"Original Field", 1);
    viewer.set_cartesian_field(rawFieldConjugate,"Conjugate Field", 0);
    
    Eigen::MatrixXd extField(mesh.V.rows(), 3*N);
    extField<<mesh.minVertexPrincipalDirections, mesh.maxVertexPrincipalDirections, -mesh.minVertexPrincipalDirections, -mesh.maxVertexPrincipalDirections;
    viewer.set_raw_field(mesh.V, extField, "Principal directions Field", 1, 0.3*mesh.avgEdgeLength);
    viewer.set_surface_vertex_data(confidence, "Confidence function");
    //viewer.set_surface_vertex_data(mesh.vertexPrincipalCurvatures.col(1), "Max curvature", 0);
    viewer.set_field_color({107.0/255.0, 8.0/255.0, 125.0/255.0}, 0);
    viewer.set_callback(callbackFunc);
    viewer.launch();
}
