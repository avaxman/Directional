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
                               const Eigen::VectorXi& cut2wholeIndices,
                               const Eigen::VectorXi& edge2TransitionIndices,
                               const Eigen::SparseMatrix<double> vt2cMat,
                               const Eigen::SparseMatrix<double> constraintMat,
                               Eigen::MatrixXd& cutUV)
  

  {
    using namespace Eigen;
    using namespace std;
    
    //TODO: in vertex sapce, not corner...
    int N = rawField.cols()/3;
    //constructing face differentials
    Triplet<double> d0Triplets;
    Triplet<double> M1Triplets;
    VectorXd gamma(3*N*wholeF.rows());;
    for (int i=0;i<wholeF.rows();i++){
      for (int j=0;j<3;j++){
        for (int k=0;k<N;k++){
          d0Triplets.push_back(Triplet<double>(3*N*i+N*j+k, 3*N*i+N*j+k, -1.0));
          d0Triplets.push_back(Triplet<double>(3*N*i+N*j+k, 3*N*i+N*(j+1)%3+k, 1.0));
          Rowvector3d edgeVector=wholeV.row(F(i,(j+1)%3))-wholeV.row(F(i,j));
          gamma(3*N*i+N*j+k)=rawField.block(i, 3*k, 1,3).dot(edgeVector);
          M1Triplets.push_back(Triplet<double>(3*N*i+N*j+k, 3*N*i+N*j+k, edgeWeights(FE(i,j))));
        }
      }
    }
    SparseMatrix<double> d0(3*N*wholeF.rows(), 3*N*wholeF.rows());
    d0.setFromTriplets(d0Triplets.begin(), d0Triplets.end());
    SparseMatrix<double> M1(3*N*wholeF.rows(), 3*N*wholeF.rows());
    M1.setFromTriplets(M1Triplets.begin(), M1Triplets.end());
    
    SparseMatrix<double> d0T=d0.transpose();
    EtE=d0T*M1*d0;
    
    SparseMatrix A(3*N*wholeF.rows()+constraintMat.rows());
    
    vector<Triplet<double>> ATriplets;
    for (int k=0; k<EtE.outerSize(); ++k)
      for (SparseMatrix<double>::InnerIterator it(EtE,k); it; ++it)
        ATriplets.push_back(Triplet<double>(it.row(), it.col(), it.value()));
    
    for (int k=0; k<constraintMat.outerSize(); ++k){
      for (SparseMatrix<double>::InnerIterator constraintMat(EtE,k); it; ++it){
        ATriplets.push_back(Triplet<double>(it.row()+3*N*wholeF.rows(), it.col(), it.value()));
        ATriplets.push_back(Triplet<double>(it.col(), it.row()+3*N*wholeF.rows(), it.value()));
      }
    }
    
    VectorXd b=VectorXd::Zero(3*N*wholeF.rows()+constraintMat.rows());
    b.segment(0,3*N*wholeF.rows())=d0T*gamma;
    
    VectorXd x = A.ldlt().solve(b);
    
    
  }
}




#endif


