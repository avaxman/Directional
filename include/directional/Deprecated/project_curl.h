// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_PROJECT_CURL
#define DIRECTIONAL_PROJECT_CURL

#include <Eigen/Core>
#include <vector>
#include <set>
#include <directional/PCFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/curl_matrices.h>
#include <directional/sparse_block.h>
#include <directional/extrinsic_intrinsic_matrices.h>

namespace directional {

inline void gradient_descent(const Eigen::SparseMatrix<double>& A,
                             const Eigen::VectorXd& gradMasses,
                             const Eigen::VectorXd& b,
                             const Eigen::VectorXd& initx,
                             const double tol,
                             const int maxIterations,
                             Eigen::VectorXd& resultx)
{
    
    using namespace Eigen;
    using namespace std;
    
    resultx = initx;
    double denValue = b.norm();
    //cout<<"(A*initx-b).lpNorm<Infinity>(): "<<(A*initx-b).lpNorm<2>();
    for (int i=0;i<maxIterations;i++){
        VectorXd Ax=A*resultx;
        VectorXd Axb = Ax-b;
        double currValue = Axb.norm()/denValue;
        // cout<<"currValue: "<<currValue<<endl;
        if (currValue<tol)
            break;
        
        VectorXd gradient = gradMasses.array()*Axb.array();
        VectorXd Ag = A*gradient;
        
        double alpha = (Ag.dot(resultx))/(Ag.dot(gradient));
        resultx = resultx - alpha*gradient;
        //cout<<"(A*resultx-b).lpNorm<Infinity>(): "<<(A*resultx-b).lpNorm<Infinity>();
    }
    
}


//This only works with face based raw fields, with a given matching
//Reducing the curl of a field by solving for the closest raw field that is curl free.
//the optional objMatrix is in case we want to minimize an objective (default: closeness)
//The optional reducMatrix is if the field has reduced degrees of freedom (for instance, symmetry). We have that field = reducMatrix*trueDofField
//Currently hard constraints are ignored
inline void project_curl(const CartesianField& origField,
                         const Eigen::VectorXi& constFaces,   //these are only in case of hard constraints, otherwise leave empty (soft constraints should be baked into objMatrix
                         const Eigen::MatrixXd& constVectors,
                         CartesianField& curlFreeField,
                         const Eigen::SparseMatrix<double>& objMatrix=Eigen::SparseMatrix<double>(),
                         const Eigen::VectorXd& objRhs=Eigen::VectorXd(),
                         const Eigen::SparseMatrix<double>& reducMatrix = Eigen::SparseMatrix<double>()){
    
    using namespace Eigen;
    using namespace std;
    
    assert(origField.fieldType == fieldTypeEnum::RAW_FIELD && origField.tb->discTangType() ==discTangTypeEnum::FACE_SPACES && "project_curl(): field should be a face-based raw field!");
    
    VectorXd rawFieldVec = origField.flatten(true);
    
    VectorXi constinFace(constFaces.size());
    
    //Finding the closest const field in face
    for (int i=0;i<constFaces.size();i++){
        double maxDot=-32767.0;
        //cout<<"b.row(i): "<<b.row(i)<<endl;
        for (int j=0;j<origField.N;j++){
            RowVector3d currVector3D = origField.extField.block(constFaces(i), 3*j,1,3).normalized();
            //cout<<"currVector3D: "<<currVector3D<<endl;
            double currDot = currVector3D.dot(constVectors.row(i));
            if (currDot>maxDot){
                maxDot=currDot;
                constinFace(i)=j;
            }
        }
        //cout<<"maxDot: "<<maxDot<<endl;
    }
    
    SparseMatrix<double> E;
    VectorXd rhs;
    if (objMatrix.nonZeros()==0){
        E.resize(rawFieldVec.size(), rawFieldVec.size());
        E.setIdentity();
        rhs = rawFieldVec;
    } else {
        E = objMatrix;
        rhs = objRhs;
    }
    
    SparseMatrix<double> R;
    if (reducMatrix.nonZeros()==0){
        R.resize(rawFieldVec.size(), rawFieldVec.size());
        R.setIdentity();
    } else R = reducMatrix;
    
    //TODO: hard constraints
    
    SparseMatrix<double> C = directional::curl_matrix_2D<double>(*((directional::PCFaceTangentBundle*)(origField.tb))->mesh, origField.matching, true, origField.N);
    SparseMatrix<double> CR=C*R;
    SparseMatrix<double> CRt=CR.transpose();
    SparseMatrix<double> ER=E*R;
    
    //cout<<"(C*rawFieldVec).lpNorm<Infinity>() before: "<<(C*rawFieldVec).lpNorm<Infinity>()<<endl;
    
    Eigen::Matrix2i orderMat; orderMat<<0,1,2,3;
    std::vector<Eigen::SparseMatrix<double>> matVec;
    matVec.push_back(ER.transpose()*ER);
    matVec.push_back(CRt);
    matVec.push_back(CR);
    matVec.push_back(Eigen::SparseMatrix<double>(CR.rows(), CR.rows()));
    Eigen::SparseMatrix<double> A;
    directional::sparse_block(orderMat, matVec, A);
    Eigen::Vector<double, Eigen::Dynamic> b(ER.rows()+CR.rows());
    b.head(ER.rows()) = ER.transpose()*rhs;
    b.tail(CR.rows()).setZero();
    
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(A);
    assert("project_curl(): Decomposition failed!" && solver.info() == Eigen::Success);
    Eigen::Vector<double, Eigen::Dynamic> x = solver.solve(b);
    assert("project_curl(): Solver failed!" && solver.info() == Eigen::Success);
    
    
    VectorXd cfFieldVec = R*x.head(ER.rows());
    //cout<<"(C*rawFieldVec).lpNorm<Infinity>() after: "<<(C*cfFieldVec).lpNorm<Infinity>()<<endl;
    
    curlFreeField=origField;
    Eigen::SparseMatrix<double> IE = directional::face_intrinsic_to_extrinsic_matrix_2D<double>(*((directional::PCFaceTangentBundle*)(origField.tb))->mesh, origField.N);
    curlFreeField.set_extrinsic_field(IE*cfFieldVec);
    
    //checking with extrinsic matrix
    SparseMatrix<double> Cext = directional::curl_matrix_2D<double>(*((directional::PCFaceTangentBundle*)(origField.tb))->mesh, origField.matching, false, origField.N);
    VectorXd v = Cext * IE * cfFieldVec;
    
    Eigen::Index maxIndex;
    double maxAbsVal = v.cwiseAbs().maxCoeff(&maxIndex);
    
    /*cout << "Extrinsic curl max abs: " << maxAbsVal << endl;
    cout << "At index: " << maxIndex << endl;
    cout << "Actual value: " << v[maxIndex] << endl;
    cout << "2D value: " << (C*cfFieldVec)[maxIndex] <<endl;
    cout<<"Extrinsic curl: "<<(Cext*IE*cfFieldVec).lpNorm<Infinity>()<<endl;*/
    
}
};


#endif
