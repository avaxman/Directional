#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
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
Eigen::MatrixXd constVectors;
double smoothWeight, roSyWeight, alignWeight;
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
        viewer.set_cartesian_field(rawFieldConjugate,"Conjugate Field", 2);
        viewer.set_field_color({107.0/255.0, 8.0/255.0, 125.0}, 2);
    }
    ImGui::PopItemWidth();
}






int main()
{
    // Load mesh
    directional::readOFF(TUTORIAL_DATA_PATH "/botanic-garden-bubble.off", mesh);
    vtb.init(mesh);
    pvFieldConjugate.init(vtb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
    pvFieldOrig.init(vtb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
    
    //discovering and constraining sharp edges
    std::vector<int> constVerticeslist;
    std::vector<Eigen::Vector3d> constVectorslist;
    for (int i=0;i<mesh.boundEdges.rows();i++){
        int vertexIndex = (mesh.EF(mesh.boundEdges(i),0)==-1 ? mesh.EV(mesh.boundEdges(i),1) : mesh.EV(mesh.boundEdges(i),0));
        constVerticeslist.push_back(vertexIndex);
        constVectorslist.push_back((mesh.V.row(mesh.EV(mesh.boundEdges(i),0))-mesh.V.row(mesh.EV(mesh.boundEdges(i),1))).normalized());
    }
    
    constVertices.resize(constVerticeslist.size());
    constVectors.resize(constVectorslist.size(),3);
    Eigen::MatrixXd constSources(constVectorslist.size(),3);  //For visualization
    for (int i=0;i<constVertices.size();i++){
        constVertices(i)=constVerticeslist[i];
        constVectors.row(i)=constVectorslist[i];
        constSources.row(i) = mesh.V.row(constVertices(i));
    }
    
    smoothWeight = 1.0;
    roSyWeight = 10.0;
    alignWeight = 0.001;
    
    
    pvData.N = N;
    pvData.tb = &vtb;
    pvData.verbose = true;
    pvData.constSpaces = constVertices;
    pvData.constVectors = constVectors;
    pvData.wAlignment = alignWeight*Eigen::VectorXd::Constant(constVertices.size(),1.0);
    pvData.wSmooth = smoothWeight;
    pvData.wRoSy = roSyWeight;
    
    //Computing regular PolyuVector field without iterations
    directional::polyvector_field(pvData, pvFieldOrig);
    directional::polyvector_to_raw(pvFieldOrig, rawFieldOrig, N%2==0);
    directional::principal_matching(rawFieldOrig);
    
    //Iterating for a conjugate field
    pvData.iterationMode = true;
    pvData.initImplicitFactor = 1.0;
    pvData.implicitScheduler = 0.9;
    iterationFunctions.push_back(directional::conjugate);
    directional::polyvector_field(pvData, pvFieldConjugate);
    rawFieldConjugate = rawFieldOrig;
    
    //Visualization
    viewer.init();
    viewer.set_surface_mesh(mesh,0);
    viewer.set_raw_field(constSources, constVectors, mesh.avgEdgeLength, "Constraints",  0);
    viewer.set_field_color(directional::DirectionalViewer::default_vector_constraint_color());
    viewer.set_cartesian_field(rawFieldOrig,"Original Field", 1);
    Eigen::MatrixXd extField(mesh.V.rows(), 3*N);
    extField<<mesh.minVertexPrincipalDirections, mesh.maxVertexPrincipalDirections, -mesh.minVertexPrincipalDirections, -mesh.maxVertexPrincipalDirections;
    viewer.set_cartesian_field(rawFieldConjugate,"Conjugate Field", 2);
    viewer.set_raw_field(mesh.V, extField,1.0, "Principal directions Field", 3);
    viewer.set_vertex_data(mesh.vertexPrincipalCurvatures.col(0), mesh.vertexPrincipalCurvatures.col(0).minCoeff(), mesh.vertexPrincipalCurvatures.col(0).maxCoeff(), "Min curvature", 0);
    viewer.set_vertex_data(mesh.vertexPrincipalCurvatures.col(1), mesh.vertexPrincipalCurvatures.col(1).minCoeff(), mesh.vertexPrincipalCurvatures.col(1).maxCoeff(), "Max curvature", 0);
    viewer.set_field_color({107.0/255.0, 8.0/255.0, 125.0/255.0}, 2);
    viewer.set_callback(callbackFunc);
    viewer.launch();
}
