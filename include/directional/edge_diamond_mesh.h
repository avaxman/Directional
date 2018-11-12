// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_EDGE_DIAMOND_MESH_H
#define DIRECTIONAL_EDGE_DIAMOND_MESH_H
#include <igl/igl_inline.h>
#include <igl/barycenter.h>
#include <igl/parula.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace directional
{
  // returns a triangle mesh s.t. every face is tesselated by a central vertex, and then every edge is supported by two triangles. The purpose is to visualize edge-based quantities)
  // Inputs:
  //  V  eigen double matrix     #V by 3 - vertex coordinates
  //  D  eigen int vector        #F by 1 - face degrees
  //  F  eigen int matrix        #F by max(D) - vertex indices in face
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
                                    const Eigen::MatrixXi& FE,
                                    const Eigen::MatrixXi& EF,
                                    const Eigen::VectorXd& edgeDiamondFunc,
                                    Eigen::MatrixXd& VEDMesh,
                                    Eigen::MatrixXi& FEDMesh,
                                    Eigen::MatrixXd& CEDMesh)
  {
    using namespace Eigen;
    MatrixXd faceCenters;
    igl::barycenter(V,F,faceCenters);
    
    VEDMesh.conservativeResize(V.rows()+F.rows(),3);
    VEDMesh<<V, faceCenters;
    
    VectorXd funcScalars;
    std::vector<RowVector3i> FEDList;
    std::vector<double> funcScalarList;
    for (int i=0;i<EV.rows();i++){
      int f=EF(i,0);
      int g=EF(i,1);
      
      if (f!=-1){
        RowVector3i Tf; Tf<<EV(i,0), EV(i,1), V.rows()+f;
        FEDList.push_back(Tf);
        funcScalarList.push_back(edgeDiamondFunc(i));
      }
      
      if (g!=-1){
        RowVector3i Tg; Tg<<EV(i,1), EV(i,0), V.rows()+g;
        FEDList.push_back(Tg);
        funcScalarList.push_back(edgeDiamondFunc(i));
      }
    }
    
    FEDMesh.conservativeResize(FEDList.size(),3);
    funcScalars.resize(funcScalarList.size());
    for (int i=0;i<FEDMesh.rows();i++){
      FEDMesh.row(i)<<FEDList[i];
      funcScalars(i)=funcScalarList[i];
    }
    
    igl::parula(funcScalars, funcScalars.minCoeff(), funcScalars.maxCoeff(), CEDMesh);
    
    return true;
  }
}


#endif


