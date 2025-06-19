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
directional::CartesianField pvFieldConjugate, rawFieldConjugate,constraintsField, rawFieldOrig, pvFieldOrig;
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
    directional::readOFF(TUTORIAL_DATA_PATH "/botanic-garden-bubble.off", mesh);
    ftb.init(mesh);
    pvFieldConjugate.init(ftb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
    pvFieldOrig.init(ftb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
    
    //discovering and constraining sharp edges
    std::vector<int> constFaceslist;
    std::vector<Eigen::Vector3d> constVectorslist;
    for (int i=0;i<mesh.boundEdges.rows();i++){
        int faceIndex = (mesh.EF(mesh.boundEdges(i),0)==-1 ? mesh.EF(mesh.boundEdges(i),1) : mesh.EF(mesh.boundEdges(i),0));
        constFaceslist.push_back(faceIndex);
        constVectorslist.push_back((mesh.V.row(mesh.EV(mesh.boundEdges(i),0))-mesh.V.row(mesh.EV(mesh.boundEdges(i),1))).normalized());
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
    roSyWeight = 10.0;
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
    directional::principal_matching(rawFieldOrig);
    
    //Iterating for a conjugate field
    pvData.numIterations = 25;
    pvData.initImplicitFactor = 1.0;
    std::vector<directional::PvIterationFunction> iterationFunctions;
    iterationFunctions.push_back(directional::conjugate);
    //iterationFunctions.push_back(directional::curl_projection);
    directional::polyvector_field(pvData, pvFieldConjugate, iterationFunctions);
    directional::polyvector_to_raw(pvFieldConjugate, rawFieldConjugate, N%2==0);
    directional::principal_matching(rawFieldConjugate);
    
    //Visualization
    viewer.init();
    viewer.set_surface_mesh(mesh,0);
    viewer.set_raw_field(constSources, constVectors, mesh.avgEdgeLength, "Constraints",  0);
    viewer.set_field_color(directional::DirectionalViewer::default_vector_constraint_color());
    viewer.highlight_faces(constFaces,"Const Faces", 0);
    viewer.toggle_field_highlight(true,0);
    viewer.set_cartesian_field(rawFieldOrig,"Original Field", 1);
    //viewer.set_raw_field(mesh.barycenters, mesh.minFacePrincipalDirectionals, mesh.avgEdgeLength, "Min Curvature Field",  2);
    viewer.set_cartesian_field(rawFieldConjugate,"Conjugate Field", 2);
    viewer.set_field_color({107.0/255.0, 8.0/255.0, 125.0}, 2);
    //viewer.set_raw_field(mesh.barycenters, mesh.maxFacePrincipalDirectionals, mesh.avgEdgeLength, "Max Curvature Field",  3);
    //viewer.set_field_color({125.0, 107.0/255.0, 8.0/255.0}, 3);
    //viewer.set_cartesian_field(rawFieldCurlFree,"Curl-free Field", 2);
    //viewer.set_field_color({107.0/255.0, 8.0/255.0, 125.0}, 2);
    //viewer.set_edge_data(curlOrig, curlOrig.cwiseAbs().minCoeff(), curlOrig.cwiseAbs().maxCoeff(), "Original Abs Curl", 0);
    //viewer.set_edge_data(curlCF, curlOrig.cwiseAbs().minCoeff(), curlOrig.cwiseAbs().maxCoeff(), "Optimized Abs Curl", 0);
    viewer.launch();
}
