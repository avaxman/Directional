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
  //  VMidEdge:   #E x 3 the vertices of the midedge mesh.
  //  FMidEdge:   #F x 3 the faces of the midedge mesh.
  
  IGL_INLINE void midedge_mesh(const Eigen::MatrixXd& VMesh,
                               const Eigen::MatrixXi& FMesh,
                               const Eigen::MatrixXi& EV,
                               const Eigen::MatrixXi& EF,
                               const Eigen::MatrixXi& FE,
                               Eigen::MatrixXd& VMidEdge,
                               Eigen::MatrixXi& FMidEdge)
  {
    
    using namespace Eigen;
    using namespace std;
    
    FMidEdge = FE;
    VMidEdge.conservativeResize(EV.rows(),3);
    for (int i=0;i<EV.rows();i++)
      VMidEdge.row(i) = (VMesh.row(EV(i,0))+VMesh.row(EV(i,1)))/2.0;
  }
}




#endif


