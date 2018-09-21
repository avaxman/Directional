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
  // Creates a parameterization of (currently supported) (u,v) functions from a directional field by solving the Poisson equation, with custom edge weights
  // Input:
  //  wholeV:      #V x 3 vertex coordinates of the original mesh
  //  wholeF:      #F x 3 face vertex indices of the original mesh
  //  FE:          #F x 3 faces to edges indices
  //  rawField:    #F by 3*N  The directional field, assumed to be ordered CCW, and in xyzxyz raw format. The degree is inferred by the size. (currently only supporting sign-symmetric 4-fields)
  //  edgeWeights  #E x 3  weight per edge (for the Poisson system)
  //  edgeLength   #edgeLength of quad mesh (scaling the gradient)
  //  vt2cMat:     #V+#T (translational jumps)  x #cutV - a map between whole vertex values and the translational jump to the vertices of the cut mesh
  //  constraintMat: matrix of constraints around singularities and nodes in the cut graph
  //  cutV:        #cV x 3 vertices of the cut mesh
  //  cutF:      #F x 3 faces of the cut mesh
  // integerRound;   which variables (from #V+#T) are rounded iteratively to integers. for each "x" entry that means that the [4*x,4*x+4] entries of vt will be integer
  //  Output:
  // cutUV:        #cV x 2 (u,v) coordinates per cut vertex
  IGL_INLINE void parameterize(const Eigen::MatrixXd& wholeV,
                               const Eigen::MatrixXi& wholeF,
                               const Eigen::MatrixXi& FE,
                               const Eigen::MatrixXd rawField,
                               const Eigen::VectorXd& edgeWeights,
                               const double edgeLength,
                               //const Eigen::SparseMatrix<double> i2vtMat,
                               const Eigen::SparseMatrix<double> vt2cMat,
                               const Eigen::SparseMatrix<double> symmMat,
                               const Eigen::SparseMatrix<double> constraintMat,
                               const Eigen::MatrixXd& cutV,
                               const Eigen::MatrixXi& cutF,
                               const Eigen::VectorXi& integerVars,
                               Eigen::MatrixXd& cutUV)
  
  
  {
    using namespace Eigen;
    using namespace std;
    
    //TODO: in vertex space, not corner...
    int N = rawField.cols()/3;
    int numVars = symmMat.cols();
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
    
    //the variables that should be fixed in the end
    VectorXi fixedMask(numVars);
    for (int i=0;i<N/2;i++)
      fixedMask(i)=1;  //first vertex is always (0,0)
    for (int i=0;i<integerVars.size();i++)
      fixedMask(i)=1;
    
    //the variables that were already fixed in the next iteration
    VectorXi alreadyFixed(numVars);
    for (int i=0;i<N/2;i++)
      alreadyFixed(i)=1;  //first vertex is always (0,0)
    
    //the values for the fixed variables (size is as all variables)
    VectorXd fixedValues(numVars);
    fixedValues.setZero();
    
    SparseMatrix<double> Efull=d0*vt2cMat*symmMat;
    do{
      //the non-fixed variables to all variables
      SparseMatrix<double> var2AllMat(numVars, numVars-alreadyFixed.sum());
      int varCounter=0;
      vector<Triplet<double> > var2AllTriplets;
      for (int i=0;i<numVars;i++)
        if (!fixedMask(i))
          var2AllTriplets.push_back(Triplet<double>(i, varCounter++, 1.0));
      var2AllMat.setFromTriplets(var2AllTriplets.begin(), var2AllTriplets.end());
      
      SparseMatrix<double> Epart = Efull*var2AllMat;
      VectorXd torhs = Efull*fixedValues;
      SparseMatrix<double> EtE = Epart.transpose()*M1*Epart;
      SparseMatrix<double> C = constraintMat*(symmMat*var2AllMat);
      
      //reducing constraintMat
      SparseQR<SparseMatrix<double>, COLAMDOrdering<int> > qrsolver;
      qrsolver.compute(C.transpose());
      int CRank=qrsolver.rank();
      //cout<<"CRank: "<<CRank<<endl;
      
      //creating sliced permutation matrix
      VectorXi PIndices=qrsolver.colsPermutation().indices();
      //cout<<"PIndices: "<<PIndices<<endl;
      
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
      //myfile<<igl::matlab_format(C,"C")<<std::endl;
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
      b.segment(0,EtE.rows())=Epart.transpose()*M1*(gamma+torhs);
      //cout<<"gamma: "<<gamma<<endl;
      
      SimplicialLDLT<SparseMatrix<double> > ldltsolver;
      //cout<<"Computing A..."<<endl;
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
      
      cout<<"(C*x.head(C.cols())).lpNorm<Infinity>(): "<<(C*x.head(C.cols())).lpNorm<Infinity>()<<endl;
      cout<<"(d0*(vt2cMat*firstVertexZeroMat*SymmMat)*x.head(C.cols())-gamma).lpNorm<Infinity>(): "<<(d0*(vt2cMat*firstVertexZeroMat*SymmMat)*x.head(C.cols())-gamma).lpNorm<Infinity>()<<endl;
    }while(!sum(alreadyFixed-fixedMask)!=0);
    
   
    //the results are packets of N functions for each vertex, and need to be allocated for corners
    VectorXd cutUVVec=vt2cMat*firstVertexZeroMat*SymmMat*x.head(EtE.cols());
    cutUV.conservativeResize(cutV.rows(),N/2);
    for (int i=0;i<cutV.rows();i++)
      cutUV.row(i)<<cutUVVec.segment(N*i,N/2).transpose();
  }
  
  
}

#endif


