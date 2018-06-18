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
#include <igl/matlab/MatlabWorkspace.h>

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
                               const Eigen::SparseMatrix<double> vt2cMat,
                               const Eigen::SparseMatrix<double> constraintMat,
                               Eigen::MatrixXd& cornerUV)
  

  {
    using namespace Eigen;
    using namespace std;
    
    //TODO: in vertex space, not corner...
    int N = rawField.cols()/3;
    //constructing face differentials
    vector<Triplet<double>> d0Triplets;
    vector<Triplet<double>> M1Triplets;
    VectorXd gamma(3*N*wholeF.rows());
    for (int i=0;i<wholeF.rows();i++){
      for (int j=0;j<3;j++){
        for (int k=0;k<N;k++){
          d0Triplets.push_back(Triplet<double>(3*N*i+N*j+k, 3*N*i+N*j+k, -1.0));
          d0Triplets.push_back(Triplet<double>(3*N*i+N*j+k, 3*N*i+N*(j+1)%3+k, 1.0));
          Vector3d edgeVector=(wholeV.row(wholeF(i,(j+1)%3))-wholeV.row(wholeF(i,j))).transpose();
          gamma(3*N*i+N*j+k)=(rawField.block(i, 3*k, 1,3)*edgeVector)(0,0);
          M1Triplets.push_back(Triplet<double>(3*N*i+N*j+k, 3*N*i+N*j+k, edgeWeights(FE(i,j))));
        }
      }
    }
    SparseMatrix<double> d0(3*N*wholeF.rows(), 3*N*wholeF.rows());
    d0.setFromTriplets(d0Triplets.begin(), d0Triplets.end());
    SparseMatrix<double> M1(3*N*wholeF.rows(), 3*N*wholeF.rows());
    M1.setFromTriplets(M1Triplets.begin(), M1Triplets.end());
    
    SparseMatrix<double> d0T=d0.transpose();
    SparseMatrix<double> vt2cMatTranspose=vt2cMat.transpose();

    SparseMatrix<double> removeFirstVertexMat(vt2cMat.cols()/2, (vt2cMat.cols()-N)/2);
    vector<Triplet<double> > removeFirstVertexMatTriplets;
    for (int i=0;i<(vt2cMat.cols()-N)/2;i++)
      removeFirstVertexMatTriplets.push_back(Triplet<double>(i+N/2, i, 1.0));
    removeFirstVertexMat.setFromTriplets(removeFirstVertexMatTriplets.begin(), removeFirstVertexMatTriplets.end());
    SparseMatrix<double> EtE=d0T*M1*d0;
    
    //filtering out symmetry - only for N=4!
    SparseMatrix<double> SymmMat(vt2cMat.cols(), vt2cMat.cols()/2);
    vector<Triplet<double>> SymmMatTriplets;
    for (int i=0;i<vt2cMat.cols();i+=4){
      SymmMatTriplets.push_back(Triplet<double>(i, i/2, 1.0));
      SymmMatTriplets.push_back(Triplet<double>(i+1, i/2+1, 1.0));
      SymmMatTriplets.push_back(Triplet<double>(i+2, i/2, -1.0));
      SymmMatTriplets.push_back(Triplet<double>(i+3, i/2+1, -1.0));
    }
    
    SymmMat.setFromTriplets(SymmMatTriplets.begin(), SymmMatTriplets.end());
    
    EtE = removeFirstVertexMat.transpose()*SymmMat.transpose()*vt2cMat.transpose()*EtE*vt2cMat*SymmMat*removeFirstVertexMat;
    SparseMatrix<double> C = constraintMat*SymmMat*removeFirstVertexMat;
    
    SparseMatrix<double> A(EtE.rows()+C.rows(),EtE.rows()+C.rows());
    
    vector<Triplet<double>> ATriplets;
    for (int k=0; k<EtE.outerSize(); ++k)
      for (SparseMatrix<double>::InnerIterator it(EtE,k); it; ++it)
        ATriplets.push_back(Triplet<double>(it.row(), it.col(), it.value()));
    
    for (int k=0; k<constraintMat.outerSize(); ++k){
      for (SparseMatrix<double>::InnerIterator it(constraintMat,k); it; ++it){
        ATriplets.push_back(Triplet<double>(it.row()+EtE.rows(), it.col(), it.value()));
        ATriplets.push_back(Triplet<double>(it.col(), it.row()+EtE.rows(), it.value()));
      }
    }
    
    VectorXd b=VectorXd::Zero(EtE.rows()+C.rows());
    b.segment(0,EtE.rows())=removeFirstVertexMat.transpose()*SymmMat.transpose()*vt2cMatTranspose*d0T*M1*gamma;
    
    igl::matlab::MatlabWorkspace mw;
    mw.save(A,"A");
    mw.save_index(b,"b");
    mw.write("sphere.mat");
    
    SimplicialLDLT<SparseMatrix<double> > solver;
    cout<<"Computing A..."<<endl;
    solver.compute(A);
    if(solver.info()!=Success) {
      cout<<"Compute failed!!!"<<endl;
      return;
    }
     cout<<"Computing A done!"<<endl;
    VectorXd x = solver.solve(b);
    if(solver.info()!=Success) {
      cout<<"Solving failed!!!"<<endl;
      return;
    }
    

    //the results are packets of N functions for each vertex, and need to be allocated for corners
    VectorXd cornerUVVec=VectorXd::Zero(N*3*wholeF.rows(),1);
    cornerUVVec.tail(3*N*wholeF.rows()-N)=(vt2cMat*SymmMat*x);
  }
}




#endif


