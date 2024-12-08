#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/write_raw_field.h>
#include <directional/polyvector_to_raw.h>
#include <directional/raw_to_polyvector.h>
#include <directional/polyvector_field.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/curl_matching.h>
#include <directional/extrinsic_intrinsic_matrices.h>

Eigen::VectorXi constFaces;
directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField pvFieldSoft, rawFieldSoft,constraintsField, rawFieldOrig, pvFieldOrig;
Eigen::MatrixXd constVectors;
Eigen::VectorXd curlOrig,curlCF;

double smoothWeight, roSyWeight, alignWeight;

directional::DirectionalViewer viewer;

int N = 1;
typedef enum {CONSTRAINTS, HARD_PRESCRIPTION, SOFT_PRESCRIPTION} ViewingModes;
ViewingModes viewingMode=CONSTRAINTS;


int main()
{
    // Load mesh
    directional::readOFF(TUTORIAL_DATA_PATH "/cheburashka.off", mesh);
    ftb.init(mesh);
    pvFieldSoft.init(ftb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
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

    directional::polyvector_field(ftb, constFaces, constVectors, smoothWeight, roSyWeight, alignWeight*Eigen::VectorXd::Constant(constFaces.size(),1.0), N, pvFieldOrig, false, false);
    directional::polyvector_field(ftb, constFaces, constVectors, smoothWeight, roSyWeight, alignWeight*Eigen::VectorXd::Constant(constFaces.size(),1.0), N, pvFieldSoft, true, true, 1);
    directional::polyvector_to_raw(pvFieldOrig, rawFieldOrig, N%2==0);
    //std::cout<<rawFieldOrig.matching<<std::endl;
    directional::polyvector_to_raw(pvFieldSoft, rawFieldSoft, N%2==0);

    //testing extrinsic to intrinsic matrix
    Eigen::VectorXd extVec = rawFieldOrig.flatten(false);
    Eigen::VectorXd intVec = rawFieldOrig.flatten(true);
    Eigen::SparseMatrix<double> EI = directional::face_extrinsic_to_intrinsic_matrix_2D<double>(&mesh, N);
    std::cout<<"IE difference: "<<(EI*extVec - intVec).cwiseAbs().maxCoeff()<<std::endl;
    Eigen::SparseMatrix<double> IE = directional::face_intrinsic_to_extrinsic_matrix_2D<double>(&mesh, N);
    std::cout<<"IE difference: "<<(IE*intVec - extVec).cwiseAbs().maxCoeff()<<std::endl;

    //std::cout<<"Raw Field first: "<<rawFieldOrig.extField<<std::endl;
    //std::cout<<"Optimized Field first: "<<rawFieldSoft.extField.row(0)<<std::endl;

    directional::principal_matching(rawFieldOrig);
    //std:: cout<<"Max curl original: "<<curlOrig.cwiseAbs().maxCoeff()<<std::endl;
    directional::principal_matching(rawFieldSoft);
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
    viewer.set_mesh(mesh,0);
    viewer.set_field(constraintsField,"Constraints",0, 0);
    viewer.highlight_faces(constFaces,0);
    viewer.set_field(rawFieldOrig,"Original Field", 0, 0);
    viewer.set_field(rawFieldSoft,"Curl-free Field", 0, 1);
    //viewer.set_edge_data(curlOrig, curlOrig.cwiseAbs().minCoeff(), curlOrig.cwiseAbs().maxCoeff(), "Original Abs Curl", 0);
    //viewer.set_edge_data(curlCF, curlCF.cwiseAbs().minCoeff(), curlCF.cwiseAbs().maxCoeff(), "Optimized Abs Curl", 0);
    viewer.launch();
}
