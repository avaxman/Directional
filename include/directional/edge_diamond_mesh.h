// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2020 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_EDGE_DIAMOND_MESH_H
#define DIRECTIONAL_EDGE_DIAMOND_MESH_H
#include <igl/igl_inline.h>
#include <igl/barycenter.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace directional
  {
  // returns a triangle mesh s.t. every face is tesselated by a central vertex, and then every edge is supported by two triangles. The purpose is to visualize edge-based quantities)
  // Inputs:
  //  V  eigen double matrix     #V by 3 - vertex coordinates
  //  F  eigen int matrix        #F by 3 - vertex indices in face
  //  EV eigen int matrix     #E by 2 - map from edges to end vertices
  //  EF eigen int matrix     #E by 2 - map from edges to adjacent faces
  //
  // Outputs:
  //  edgeV  eigen double matrix  #F+#V by 3 - new vertices
  //  edgeT  eigen int matrix    2*#E-#Boundary - new edge-based triangles
  //  edgeTE eigen int vector     #edgeT edgeT -> original edge in EV.
  IGL_INLINE bool edge_diamond_mesh(const Eigen::MatrixXd& V,
                                    const Eigen::MatrixXi& F,
                                    const Eigen::MatrixXi& EV,
                                    const Eigen::MatrixXi& EF,
                                    Eigen::MatrixXd& edgeV,
                                    Eigen::MatrixXi& edgeT,
                                    Eigen::VectorXi& edgeTE)
  {
    using namespace Eigen;
    MatrixXd faceCenters;
    igl::barycenter(V,F,faceCenters);
    
    edgeV.resize(V.rows()+F.rows(),3);
    edgeV<<V, faceCenters;
    
    std::vector<RowVector3i> edgeTList;
    std::vector<int> edgeTEList;
    for (int i=0;i<EV.rows();i++){
      int f=EF(i,0);
      int g=EF(i,1);
      
      if (f!=-1){
        RowVector3i Tf; Tf<<EV(i,0), EV(i,1), V.rows()+f;
        edgeTList.push_back(Tf);
        edgeTEList.push_back(i);
      }
      
      if (g!=-1){
        RowVector3i Tg; Tg<<EV(i,1), EV(i,0), V.rows()+g;
        edgeTList.push_back(Tg);
        edgeTEList.push_back(i);
      }
    }
    
    edgeT.resize(edgeTList.size(),3);
    edgeTE.resize(edgeTEList.size());
    for (int i=0;i<edgeT.rows();i++){
      edgeT.row(i)<<edgeTList[i];
      edgeTE(i)=edgeTEList[i];
    }
    
    return true;
  }
  
  }


#endif



