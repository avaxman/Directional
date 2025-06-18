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
directional::CartesianField pvFieldGL, rawFieldGL,constraintsField, rawFieldOrig, pvFieldOrig;
Eigen::MatrixXd constVectors;

double smoothWeight, roSyWeight, alignWeight;

directional::DirectionalViewer viewer;

int N = 4;
typedef enum {CONSTRAINTS, HARD_PRESCRIPTION, SOFT_PRESCRIPTION} ViewingModes;
ViewingModes viewingMode=CONSTRAINTS;


int main()
{
    // Load mesh
    directional::readOFF(TUTORIAL_DATA_PATH "/deckel_input.off", mesh);
    ftb.init(mesh);
    pvFieldGL.init(ftb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
    pvFieldOrig.init(ftb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
    
    //discovering and constraining sharp edges
    std::vector<int> constFaceslist;
    std::vector<Eigen::Vector3d> constVectorslist;
    for (int i=0;i<mesh.EF.rows();i++){
        int face1Index = mesh.EF(i,0);
        int face2Index = mesh.EF(i,1);
        if (mesh.faceNormals.row(face1Index).dot(mesh.faceNormals.row(face2Index))<0.5){
            constFaceslist.push_back(face1Index);
            constFaceslist.push_back(face2Index);
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
    pvData.initImplicitFactor = 5.0;
    pvData.implicitScheduler = 1.0;
    
    //Computing regular PolyVector field without iterations
    directional::polyvector_field(pvData, pvFieldOrig);
    directional::polyvector_to_raw(pvFieldOrig, rawFieldOrig, N%2==0);
    directional::principal_matching(rawFieldOrig);
    
    //Iterating for a curl-free field
    pvData.numIterations = 25;
    std::vector<directional::PvIterationFunction> iterationFunctions;
    iterationFunctions.push_back(directional::hard_rosy);
    directional::polyvector_field(pvData, pvFieldGL, iterationFunctions);
    directional::polyvector_to_raw(pvFieldGL, rawFieldGL, N%2==0);
    //std::cout<<"pvFieldGL.intField: "<<pvFieldGL.intField<<std::endl;
    directional::principal_matching(rawFieldGL);
    
    //Visualization
    viewer.init();
    viewer.set_surface_mesh(mesh,0);
    viewer.set_raw_field(constSources, constVectors, mesh.avgEdgeLength, "Constraints",  0);
    viewer.set_field_color(directional::DirectionalViewer::default_vector_constraint_color());
    viewer.highlight_faces(constFaces,"Const Faces", 0);
    viewer.toggle_field_highlight(true,0);
    viewer.set_cartesian_field(rawFieldOrig,"Original Field", 1);
    viewer.set_cartesian_field(rawFieldGL,"Ginzburg-Landau Field", 2);
    viewer.set_field_color({107.0/255.0, 8.0/255.0, 125.0}, 2);
    viewer.launch();
}
