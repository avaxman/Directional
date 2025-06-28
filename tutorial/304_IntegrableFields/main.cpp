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
#include <directional/PCFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/curl_matching.h>
#include <directional/extrinsic_intrinsic_matrices.h>

Eigen::VectorXi constFaces;
directional::TriMesh mesh;
directional::PCFaceTangentBundle ftb;
directional::CartesianField pvFieldCurlFree, rawFieldCurlFree,constraintsField, rawFieldOrig, pvFieldOrig;
Eigen::MatrixXd constVectors;
Eigen::VectorXd curlOrig,curlCF;

double smoothWeight, roSyWeight, alignWeight;

directional::DirectionalViewer viewer;

int N = 4;
typedef enum {CONSTRAINTS, HARD_PRESCRIPTION, SOFT_PRESCRIPTION} ViewingModes;
ViewingModes viewingMode=CONSTRAINTS;


int main()
{
    // Load mesh
    directional::readOFF(TUTORIAL_DATA_PATH "/cheburashka.off", mesh);
    ftb.init(mesh);
    pvFieldCurlFree.init(ftb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
    pvFieldOrig.init(ftb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
    
    //discovering and constraining sharp edges
    std::vector<int> constFaceslist;
    std::vector<Eigen::Vector3d> constVectorslist;
    for (int i=0;i<mesh.EF.rows();i++){
        if (mesh.faceNormals.row(mesh.EF(i,0)).dot(mesh.faceNormals.row(mesh.EF(i,1)))<0.1){
            constFaceslist.push_back(mesh.EF(i,0));
            constFaceslist.push_back(mesh.EF(i,1));
            constVectorslist.push_back((mesh.V.row(mesh.EV(i,0))-mesh.V.row(mesh.EV(i,1))).normalized());
            constVectorslist.push_back((mesh.V.row(mesh.EV(i,0))-mesh.V.row(mesh.EV(i,1))).normalized());
        }
    }
    
    constFaces.resize(constFaceslist.size());
    constVectors.resize(constVectorslist.size(),3);
    Eigen::MatrixXd constSources(constVectorslist.size(),3);  //For visualization
    for (int i=0;i<constFaces.size();i++){
        constFaces(i)=constFaceslist[i];
        constVectors.row(i)=constVectorslist[i];
        constSources.row(i) = mesh.barycenters.row(constFaces(i));
    }
    
    smoothWeight = 1.0;
    roSyWeight = 1.0;
    alignWeight = 1.0;
    
    directional::PolyVectorData pvData;
    pvData.N = N;
    pvData.tb = &ftb;
    pvData.verbose = true;
    pvData.constSpaces = constFaces;
    pvData.constVectors = constVectors;
    pvData.wAlignment = alignWeight*Eigen::VectorXd::Constant(constFaces.size(),1.0);
    pvData.wSmooth = smoothWeight;
    pvData.wRoSy = roSyWeight;
    
    //Computing regular PolyuVector field without iterations
    directional::polyvector_field(pvData, pvFieldOrig);
    directional::polyvector_to_raw(pvFieldOrig, rawFieldOrig, N%2==0);
    directional::curl_matching(rawFieldOrig, curlOrig);
    
    //Iterating for a curl-free field
    int numIterations = 25;
    pvData.iterationMode = true;
    std::vector<directional::PvIterationFunction> iterationFunctions;
    iterationFunctions.push_back(directional::soft_rosy);
    iterationFunctions.push_back(directional::curl_projection);
    directional::polyvector_field(pvData, pvFieldCurlFree);
    directional::polyvector_iterate(pvData, pvFieldCurlFree, iterationFunctions, numIterations);
    directional::polyvector_to_raw(pvFieldCurlFree, rawFieldCurlFree, N%2==0);
    directional::curl_matching(rawFieldCurlFree, curlCF);
    
    //Visualization
    viewer.init();
    viewer.set_surface_mesh(mesh,0);
    viewer.set_raw_field(constSources, constVectors, mesh.avgEdgeLength, "Constraints",  0);
    viewer.set_field_color(directional::DirectionalViewer::default_vector_constraint_color());
    viewer.highlight_faces(constFaces,"Const Faces", 0);
    viewer.toggle_field_highlight(true,0);
    viewer.set_cartesian_field(rawFieldOrig,"Original Field", 1);
    viewer.set_cartesian_field(rawFieldCurlFree,"Curl-free Field", 2);
    viewer.set_field_color({107.0/255.0, 8.0/255.0, 125.0}, 2);
    viewer.set_edge_data(curlOrig, curlOrig.cwiseAbs().minCoeff(), curlOrig.cwiseAbs().maxCoeff()/10.0, "Original Abs Curl", 0);  //to increase sensitivity
    viewer.set_edge_data(curlCF, curlOrig.cwiseAbs().minCoeff(), curlOrig.cwiseAbs().maxCoeff()/10.0, "Optimized Abs Curl", 0);
    viewer.launch();
}
