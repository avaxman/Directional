// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_MIDEDGE_MESH_H
#define DIRECTIONAL_MIDEDGE_MESH_H

#include <queue>
#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <igl/igl_inline.h>
#include <igl/local_basis.h>
#include <igl/edge_topology.h>
#include <igl/parula.h>


namespace directional
{
  
  
  // Creating non-conforming mid-edge mesh, where the faces are between the midedges of each original face, and colored by a mid-edge function. This is generally only for visualization
  // Input:
  //  VMesh:      #V x 3 conforming mesh vertices
  //  FMesh:      #F x 3 conforming mesh faces
  //  EV:         #E x 2 edges to vertices indices
  //  EF:         #E x 2 edges to faces indices
  //  FE:         #F x 3 faces to edges indices
  // Output:
  //  VMidEdge:   #E x 3 the vertices of the midedge mesh.
  //  FMidEdge:   #F x 3 the faces of the midedge mesh.
  
  IGL_INLINE void non_conforming_mesh(const Eigen::MatrixXd& V,
                                      const Eigen::MatrixXi& F,
                                      const Eigen::MatrixXi& EV,
                                      const Eigen::MatrixXi& FE,
                                      const Eigen::MatrixXi& EF,
                                      const Eigen::VectorXd& midEdgeFunc,
                                      Eigen::MatrixXd& VNCMesh,
                                      Eigen::MatrixXi& FNCMesh,
                                      Eigen::MatrixXd& CNCMesh)
  {
    
    using namespace Eigen;
    using namespace std;
    
    VNCMesh.conservativeResize(6*F.rows(),3);
    FNCMesh.conservativeResize(4*F.rows(),3);
    
    Eigen::VectorXd funcScalars(6*F.rows());
    
    for (int i=0;i<F.rows();i++){
      VNCMesh.block(6*i,0,6,3)<<V.row(F(i,0)),
      (V.row(F(i,0))+V.row(F(i,1)))/2.0,
      V.row(F(i,1)),
      (V.row(F(i,1))+V.row(F(i,2)))/2.0,
      V.row(F(i,2)),
      (V.row(F(i,2))+V.row(F(i,0)))/2.0,
      
      FNCMesh.block(4*i,0,4,3)<<6*i,6*i+1, 6*i+5,
      6*i+2,6*i+3, 6*i+1,
      6*i+4,6*i+5, 6*i+3,
      6*i+1,6*i+3, 6*i+5;
      
      funcScalars.segment(6*i,6)<<midEdgeFunc(FE(i,0))+midEdgeFunc(FE(i,2))-midEdgeFunc(FE(i,1)),
      midEdgeFunc(FE(i,0)),
      midEdgeFunc(FE(i,1))+midEdgeFunc(FE(i,0))-midEdgeFunc(FE(i,2)),
      midEdgeFunc(FE(i,1)),
      midEdgeFunc(FE(i,2))+midEdgeFunc(FE(i,1))-midEdgeFunc(FE(i,0)),
      midEdgeFunc(FE(i,2));
    }
    
    igl::parula(funcScalars, funcScalars.minCoeff(), funcScalars.maxCoeff(), CNCMesh);
  }
  
  
}




#endif


