// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_HARMONIC_BASIS_H
#define DIRECTIONAL_HARMONIC_BASIS_H

#include <queue>
#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <igl/igl_inline.h>
#include <igl/diag.h>
#include <igl/local_basis.h>
#include <igl/edge_topology.h>
#include <directional/FEM_masses.h>
#include <directional/FEM_suite.h>
#include <directional/dual_cycles.h>
#include <igl/euler_characteristic.h>
#include <igl/per_face_normals.h>
#include <igl/boundary_loop.h>


namespace directional
{
  
  
  // Creating non-conforming mid-edge mesh, where the faces are between the midedges of each original face. This is generally only for visualization
  // Input:
  //  VMesh:      #V x 3 conforming mesh vertices
  //  FMesh:      #F x 3 conforming mesh faces
  //  EV:         #E x 2 edges to vertices indices
  //  EF:         #E x 2 edges to faces indices
  //  FE:         #F x 3 faces to edges indices
  // Output:
  //  Gv:    #3f x V Conforming gradient matrix, returning vector of xyzxyz per face gradient vectors
  //  Ge:    #3f x V Non-conforming gradient of the same style, but for mid-edge functions
  //  J:    #3F x 3F rotation operator [Nx] per face
  //  C:  Curl operator which is basically (JGe)^T * Mchi
  //  C:  Divergence operator which is basically Gv^T * Mchi
  
  IGL_INLINE void harmonic_basis(const Eigen::MatrixXd& V,
                                 const Eigen::MatrixXi& F,
                                 const Eigen::MatrixXi& EV,
                                 const Eigen::MatrixXi& FE,
                                 const Eigen::MatrixXi& EF,
                                 
                                 std::vector<Eigen::MatrixXd>& harmFields)
  {
    
    using namespace Eigen;
    using namespace std;
    
    SparseMatrix<double> Gv, Ge, J, Mv, Mchi, Mf, Me, C, D;
    VectorXd MvVec, MeVec, MfVec, MchiVec;
    
    directional::FEM_suite(V, F, EV, FE, EF, Gv, Ge, J, C, D);
    directional::FEM_masses(V, F, EV, FE, EF, MvVec, MeVec, MfVec, MchiVec);
    
    igl::diag(MvVec,Mv);
    igl::diag(MeVec,Me);
    igl::diag(MfVec,Mf);
    igl::diag(MchiVec,Mchi);
    
    Eigen::SparseMatrix<double> basisCycles;
    Eigen::VectorXd cycleCurvature;
    Eigen::VectorXi vertex2cycle;
    Eigen::VectorXi innerEdges;
    std::vector<std::vector<int>> boundaryLoops;
    
    dual_cycles(V,F,EV,EF,basisCycles,cycleCurvature,vertex2cycle,innerEdges);
    int eulerChar = igl::euler_characteristic(V,F);
    igl::boundary_loop(F, boundaryLoops);
    int numBoundaries=boundaryLoops.size();
    int numGenerators=2-numBoundaries-eulerChar;
    assert(numBoundaries==0 && "Currently not working with boundaries!");
    
    
     SparseMatrix<double> Lv = D*Gv;   //Gv^T * Mchi * Gv
    for (int cycle=basisCycles.rows()-numGenerators;cycle<basisCycles.rows();cycle++){
      SparseVector<double> singleCycle = basisCycles.row(cycle);
      
      //VectorXi bmask=VectorXi::Zero(V.rows());
      //VectorXd bcall=VectorXd::Zero(V.rows());
      VectorXd candidateFunc=VectorXd::Zero(V.rows());
      
      VectorXi cycleFaces=VectorXi::Zero(F.rows());
      for (SparseVector<double>::InnerIterator it(singleCycle); it; ++it)
      {
        candidateFunc(EV(it.index(),0))=(it.value() > 0 ? 1.0 : 0.0);
        candidateFunc(EV(it.index(),1))=(it.value() > 0 ? 0.0 : 1.0);
        cycleFaces(EF(it.index(),0))=1;
        cycleFaces(EF(it.index(),1))=1;
      }
      
      VectorXd candidateFieldVec = Gv*candidateFunc;
      for(int i=0;i<F.rows();i++)
        if (!cycleFaces(i))
          candidateFieldVec.segment(3*i,3).setZero();
      
      //solving for exact part
      VectorXd exactFunc;
      igl::min_quad_with_fixed_data<double> mqwfExact;
      // Linear term is 0
      VectorXd B = -D*candidateFieldVec;
      VectorXd Beq;
      SparseMatrix<double> Aeq;
      Eigen::VectorXi b(1); b(0)=0;
      Eigen::VectorXd bc(1); bc(0)=0;

      igl::min_quad_with_fixed_precompute(Lv,b,Aeq,true,mqwfExact);
      igl::min_quad_with_fixed_solve(mqwfExact,B,bc,Beq,exactFunc);
      
      //FIltering exact part
      VectorXd harmFieldVec = candidateFieldVec-Gv*exactFunc;
      harmFieldVec=harmFieldVec/harmFieldVec.norm();
      
      std::cout<<"harmFieldVec.norm(): "<<harmFieldVec.norm()<<std::endl;
      
      //sanity check:
      std::cout<<"(D*harmFieldVec).lpNorm<Infinity>(): "<<(D*harmFieldVec).lpNorm<Infinity>()<<std::endl;
      std::cout<<"(C*harmFieldVec).lpNorm<Infinity>(): "<<(C*harmFieldVec).lpNorm<Infinity>()<<std::endl;
      
      harmFields.push_back(Eigen::MatrixXd(F.rows(),3));
      for (int i=0;i<F.rows();i++)
        for (int j=0;j<3;j++)
          harmFields[harmFields.size()-1](i,j)=harmFieldVec(3*i+j);

    }

  }
}

#endif


