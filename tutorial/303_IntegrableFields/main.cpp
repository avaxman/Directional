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
    for (int i=0;i<constFaces.size();i++){
        constFaces(i)=constFaceslist[i];
        constVectors.row(i)=constVectorslist[i];
    }

    //generating the viewing fields
    Eigen::MatrixXd rawFieldConstraints=Eigen::MatrixXd::Zero(mesh.F.rows(),(N+3)*3);
    Eigen::VectorXi posInFace=Eigen::VectorXi::Zero(mesh.F.rows());
    for (int i=0;i<constFaces.size();i++){
        rawFieldConstraints.block(constFaces(i),3*posInFace(constFaces(i)), 1,3)=constVectors.row(i);
        posInFace(constFaces(i))++;
    }

    //Just to show the other direction if N is even, since we are by default constraining it
    if (N%2==0)
        rawFieldConstraints.middleCols(rawFieldConstraints.cols()/2, rawFieldConstraints.cols()/2)=-rawFieldConstraints.middleCols(0, rawFieldConstraints.cols()/2);

    constraintsField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, N);
    constraintsField.set_extrinsic_field(rawFieldConstraints);

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
    
    //Iterating for a curl-free field
    pvData.numIterations = 25;
    std::vector<directional::PvIterationFunction> iterationFunctions;
    iterationFunctions.push_back(directional::soft_rosy);
    iterationFunctions.push_back(directional::curl_projection);
    directional::polyvector_field(pvData, pvFieldCurlFree, iterationFunctions);
    directional::polyvector_to_raw(pvFieldCurlFree, rawFieldCurlFree, N%2==0);
    
    //std::cout<<rawFieldOrig.matching<<std::endl;
   
    //testing extrinsic to intrinsic matrix
    /*Eigen::VectorXd extVec = rawFieldOrig.flatten(false);
    Eigen::VectorXd intVec = rawFieldOrig.flatten(true);
    Eigen::SparseMatrix<double> EI = directional::face_extrinsic_to_intrinsic_matrix_2D<double>(&mesh, N);
    std::cout<<"IE difference: "<<(EI*extVec - intVec).cwiseAbs().maxCoeff()<<std::endl;
    Eigen::SparseMatrix<double> IE = directional::face_intrinsic_to_extrinsic_matrix_2D<double>(&mesh, N);
    std::cout<<"IE difference: "<<(IE*intVec - extVec).cwiseAbs().maxCoeff()<<std::endl;*/

    //std::cout<<"Raw Field first: "<<rawFieldOrig.extField<<std::endl;
    //std::cout<<"Optimized Field first: "<<rawFieldSoft.extField.row(0)<<std::endl;

    directional::principal_matching(rawFieldOrig);
    //std:: cout<<"Max curl original: "<<curlOrig.cwiseAbs().maxCoeff()<<std::endl;
    directional::principal_matching(rawFieldCurlFree);
    //std:: cout<<"Max curl optimized: "<<curlCF.cwiseAbs().maxCoeff()<<std::endl;

      //Testing for curl
    /*for (int i=0;i<mesh.innerEdges.size();i++){
        std::cout<<"Edge i : "<<i<<std::endl;
        std::cout<<"Matching: "<<rawFieldSoft.matching(i)<<std::endl;
        Eigen::RowVector3d edgeVector = mesh.V.row(mesh.EV(i,1))-mesh.V.row(mesh.EV(i,0));
        for (int n=0;n<N;n++){
            Eigen::RowVector3d vecLeft = rawFieldSoft.extField.block(mesh.EF(i,0), 3*n, 1,3);
            Eigen::RowVector3d vecRight = rawFieldSoft.extField.block(mesh.EF(i,1), 3*((n+rawFieldSoft.matching(i))%N), 1,3);
            std::cout<<"dot product left: "<<vecLeft.dot(edgeVector)<<std::endl;
            std::cout<<"dot product right: "<<vecRight.dot(edgeVector)<<std::endl;
        }

    }*/


    //triangle mesh setup
    viewer.init();
    viewer.set_surface_mesh(mesh,0);
    viewer.set_cartesian_field(constraintsField,"Constraints",0);
    viewer.highlight_faces(constFaces,"Const Faces", 0);
    viewer.toggle_field_highlight(true,0);
    viewer.set_cartesian_field(rawFieldOrig,"Original Field", 1);
    //viewer.toggle_cartesian_field(false, 1);
    viewer.set_cartesian_field(rawFieldCurlFree,"Curl-free Field", 2);
    //viewer.set_edge_data(curlOrig, curlOrig.cwiseAbs().minCoeff(), curlOrig.cwiseAbs().maxCoeff(), "Original Abs Curl", 0);
    //viewer.set_edge_data(curlCF, curlCF.cwiseAbs().minCoeff(), curlCF.cwiseAbs().maxCoeff(), "Optimized Abs Curl", 0);
    viewer.launch();
}
