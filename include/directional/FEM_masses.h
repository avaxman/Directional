// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_FEM_MASSES_H
#define DIRECTIONAL_FEM_MASSES_H

#include <queue>
#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <igl/igl_inline.h>
#include <igl/local_basis.h>
#include <igl/edge_topology.h>
#include <igl/doublearea.h>
#include <igl/massmatrix.h>


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
  //  MvVec:    #V the Voronoi areas of each vertex
  //  MeVec:    #E the diamond areas of each edge
  //  MfVec:    #F triangle areas
  //  FchiVec:  #3F triangle ares suitable for face-based field (3 times each triangel area)
  
  IGL_INLINE void FEM_masses(const Eigen::MatrixXd& V,
                             const Eigen::MatrixXi& F,
                             const Eigen::MatrixXi& EV,
                             const Eigen::MatrixXi& FE,
                             const Eigen::MatrixXi& EF,
                             Eigen::VectorXd& MvVec,
                             Eigen::VectorXd& MeVec,
                             Eigen::VectorXd& MfVec,
                             Eigen::VectorXd& MchiVec)
  {
    
    using namespace Eigen;
    using namespace std;
    
    VectorXd dblA;
    igl::doublearea(V,F,dblA);
    
    MfVec = dblA*0.5;
    MchiVec.conservativeResize(F.rows()*3);
    MvVec=VectorXd::Zero(V.rows());
    MeVec=VectorXd::Zero(EV.rows());
    for (int i=0;i<F.rows();i++){
      for (int j=0;j<3;j++){
        MchiVec(3*i+j)=MfVec(i);
        MvVec(F(i,j))+=MfVec(i)/3.0;
        MeVec(FE(i,j))+=MfVec(i)/3.0;
      }
    }
  }
}

#endif


