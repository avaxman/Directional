// This file is part of libdirectional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_PARAMETERIZE_H
#define DIRECTIONAL_PARAMETERIZE_H
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_topology.h>
#include <directional/tree.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <igl/min_quad_with_fixed.h>
#include <igl/matlab_format.h>
#include <iostream>
#include <fstream>

#include <Eigen/Core>
#include <queue>
#include <vector>
#include <cmath>


namespace directional
{
  // Reorders the vectors in a face (preserving CCW) so that the principal matching across most edges, except a small set (called a cut), is an identity, making it ready for cutting and parameterization.
  // Important: if the Raw field in not CCW ordered, the result is unpredictable.
  // Input:
  //  V:      #V x 3 vertex coordinates
  //  F:      #F x 3 face vertex indices
  //  EV:     #E x 2 edges to vertices indices
  //  EF:     #E x 2 edges to faces indices
  //  rawField: #F by 3*N  The directional field, assumed to be ordered CCW, and in xyzxyz raw format. The degree is inferred by the size.
  // Output:
  // matching: #E matching function, where vector k in EF(i,0) matches to vector (k+matching(k))%N in EF(i,1). In case of boundary, there is a -1. Expect most matching =0 due to the combing.
  //  effort: #E updated principal-matching efforts.
  IGL_INLINE void parameterize(const Eigen::MatrixXd& wholeV,
                               const Eigen::MatrixXi& wholeF,
                               const Eigen::MatrixXi& FE,
                               const Eigen::MatrixXd rawField,
                               const Eigen::VectorXd& edgeWeights,
                               const double edgeLength,
                               const Eigen::SparseMatrix<double> i2vtMat,
                               const Eigen::SparseMatrix<double> vt2cMat,
                               const Eigen::SparseMatrix<double> constraintMat,
                               const Eigen::MatrixXd& cutV,
                               const Eigen::MatrixXi& cutF,
                               Eigen::MatrixXd& cutUV)
  
  
  {
    using namespace Eigen;
    using namespace std;
    
    ofstream myfile;
    myfile.open ("parameterize.m");
    //TODO: in vertex space, not corner...
    int N = rawField.cols()/3;
    //constructing face differentials
    vector<Triplet<double>> d0Triplets;
    vector<Triplet<double>> M1Triplets;
    VectorXd gamma(3*N*wholeF.rows());
    for (int i=0;i<cutF.rows();i++){
      for (int j=0;j<3;j++){
        for (int k=0;k<N;k++){
          d0Triplets.push_back(Triplet<double>(3*N*i+N*j+k, N*cutF(i,j)+k, -1.0));
          d0Triplets.push_back(Triplet<double>(3*N*i+N*j+k, N*cutF(i,(j+1)%3)+k, 1.0));
          Vector3d edgeVector=(cutV.row(cutF(i,(j+1)%3))-cutV.row(cutF(i,j))).transpose();
          gamma(3*N*i+N*j+k)=(rawField.block(i, 3*k, 1,3)*edgeVector)(0,0)*edgeLength;
          M1Triplets.push_back(Triplet<double>(3*N*i+N*j+k, 3*N*i+N*j+k, edgeWeights(FE(i,j))));
        }
      }
    }
    SparseMatrix<double> d0(3*N*wholeF.rows(), N*cutV.rows());
    d0.setFromTriplets(d0Triplets.begin(), d0Triplets.end());
    SparseMatrix<double> M1(3*N*wholeF.rows(), 3*N*wholeF.rows());
    M1.setFromTriplets(M1Triplets.begin(), M1Triplets.end());
    
    SparseMatrix<double> d0T=d0.transpose();
    SparseMatrix<double> cutEtE=d0T*M1*d0;
    
    //setting the first variables to zero to resolve translation amiguity
    SparseMatrix<double> firstVertexZeroMat(i2vtMat.cols(), i2vtMat.cols()-N);
    vector<Triplet<double> > firstVertexZeroTriplets;
    for (int i=0;i<i2vtMat.cols()-N;i++)
      firstVertexZeroTriplets.push_back(Triplet<double>(i+N, i, 1.0));
    firstVertexZeroMat.setFromTriplets(firstVertexZeroTriplets.begin(), firstVertexZeroTriplets.end());
    
    
    myfile<<igl::matlab_format(firstVertexZeroMat,"firstVertexZeroMat")<<std::endl;
    
    //filtering out sign symmetry
    SparseMatrix<double> SymmMat(firstVertexZeroMat.cols(), firstVertexZeroMat.cols()/2);
    vector<Triplet<double>> SymmMatTriplets;
    for (int i=0;i<firstVertexZeroMat.cols();i+=N){
      for (int j=0;j<N/2;j++){
        SymmMatTriplets.push_back(Triplet<double>(i+j, i/2+j, 1.0));
        SymmMatTriplets.push_back(Triplet<double>(i+j+N/2, i/2+j, -1.0));
      }
    }
    
    SymmMat.setFromTriplets(SymmMatTriplets.begin(), SymmMatTriplets.end());
    myfile<<igl::matlab_format(SymmMat,"SymmMat")<<std::endl;
    myfile<<igl::matlab_format(cutEtE,"cutEtE")<<std::endl;
    myfile<<igl::matlab_format(constraintMat,"constraintMat")<<std::endl;
    myfile<<igl::matlab_format(vt2cMat,"vt2cMat")<<std::endl;
    myfile<<igl::matlab_format(vt2cMat,"i2vtMat")<<std::endl;
    
    SparseMatrix<double> EtE = (vt2cMat/**i2vtMat*/*firstVertexZeroMat*SymmMat).transpose()*cutEtE*(vt2cMat/**i2vtMat*/*firstVertexZeroMat*SymmMat);
    SparseMatrix<double> C = constraintMat*(/**i2vtMat*/firstVertexZeroMat*SymmMat);
    
    myfile<<igl::matlab_format(EtE,"EtE")<<std::endl;
  
    //reducing constraintMat
    SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > qrsolver;
    qrsolver.compute(C.transpose());
    int CRank=qrsolver.rank();
    cout<<"CRank: "<<CRank<<endl;
    
    //creating sliced permutation matrix
    VectorXi PIndices=qrsolver.colsPermutation().indices();
    cout<<"PIndices: "<<PIndices<<endl;
    
    vector<Triplet<double> > CTriplets;
    for (int k=0; k<C.outerSize(); ++k)
      for (SparseMatrix<double>::InnerIterator it(C,k); it; ++it)
      {
        for (int j=0;j<CRank;j++)
          if (it.row()==PIndices(j))
            CTriplets.push_back(Triplet<double>(j, it.col(), it.value()));
        
      }
    
    C.resize(CRank, C.cols());
    C.setFromTriplets(CTriplets.begin(), CTriplets.end());
     myfile<<igl::matlab_format(C,"C")<<std::endl;
    SparseMatrix<double> A(EtE.rows()+C.rows(),EtE.rows()+C.rows());
    
    vector<Triplet<double>> ATriplets;
    for (int k=0; k<EtE.outerSize(); ++k)
      for (SparseMatrix<double>::InnerIterator it(EtE,k); it; ++it)
        ATriplets.push_back(Triplet<double>(it.row(), it.col(), it.value()));
    
    for (int k=0; k<C.outerSize(); ++k){
      for (SparseMatrix<double>::InnerIterator it(C,k); it; ++it){
        ATriplets.push_back(Triplet<double>(it.row()+EtE.rows(), it.col(), it.value()));
        ATriplets.push_back(Triplet<double>(it.col(), it.row()+EtE.rows(), it.value()));
      }
    }
    
    A.setFromTriplets(ATriplets.begin(), ATriplets.end());
    
    VectorXd b=VectorXd::Zero(EtE.rows()+C.rows());
    b.segment(0,EtE.rows())=(vt2cMat/**i2vtMat*/*firstVertexZeroMat*SymmMat).transpose()*d0T*M1*gamma;
    //cout<<"gamma: "<<gamma<<endl;
    
    myfile<<igl::matlab_format(A,"A")<<std::endl;
    myfile<<igl::matlab_format(b,"b")<<std::endl;
    myfile.close();
    
    SimplicialLDLT<SparseMatrix<double> > ldltsolver;
    cout<<"Computing A..."<<endl;
    VectorXd x;
    ldltsolver.compute(A);
    if(ldltsolver.info()!=Success) {
      cout<<"LDLT failed, trying LU"<<endl;
      SparseLU<SparseMatrix<double> > lusolver;
      lusolver.compute(A);
      if(lusolver.info()!=Success) {
        cout<<"LU failed as well!"<<endl;
        return;
      }
      x = lusolver.solve(b);
    } else{
      cout<<"Computing A done!"<<endl;
      x = ldltsolver.solve(b);
      if(ldltsolver.info()!=Success) {
        cout<<"Solving failed!!!"<<endl;
        return;
      }
    }
    
    cout<<"C*x.head(C.cols()): "<<C*x.head(C.cols())<<endl;
    
    //cout<<"d0*(vt2cMat*i2vtMat*firstVertexZeroMat*SymmMat)*x.head(SymmMat.cols())-gamma: "<<gamma<<endl;
    
    
    //the results are packets of N functions for each vertex, and need to be allocated for corners
    VectorXd cutUVVec=vt2cMat*/*i2vtMat**/firstVertexZeroMat*SymmMat*x.head(EtE.cols());
    cutUV.conservativeResize(cutV.rows(),N/2);
    for (int i=0;i<cutV.rows();i++)
      cutUV.row(i)<<cutUVVec.segment(N*i,N/2).transpose();
  }
  
  
}

#endif


