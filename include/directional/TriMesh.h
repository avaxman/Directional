// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_TRIMESH_H
#define DIRECTIONAL_TRIMESH_H

#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <igl/local_basis.h>
#include <directional/polyvector_field.h>
#include <igl/per_face_normals.h>
#include <igl/edge_topology.h>

namespace directional{

class TriMesh{
public:
  
  //Basic quantities
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  
  //Topological quantities
  Eigen::MatrixXi EF, FE, EV;
  
  //Geometric quantities
  Eigen::MatrixXd faceNormals;
  Eigen::MatrixXd Bx,By;  //local basis vectors per face
  Eigen::MatrixXd barycenters;
  
  TriMesh(){}
  ~TriMesh(){}
  
  void IGL_INLINE set_mesh(const Eigen::MatrixXd& _V,
                           const Eigen::MatrixXi& _F,
                           const Eigen::MatrixXi& _EV=Eigen::MatrixXi(),
                           const Eigen::MatrixXi& _FE=Eigen::MatrixXi(),
                           const Eigen::MatrixXi& _EF=Eigen::MatrixXi()){
  
    V=_V;
    F=_F;
    if (_EV.rows()==0){
      igl::edge_topology(V,F,EV,FE,EF);
    } else{
      EV=_EV; FE=_FE; EF=_EF;
    }
    igl::barycenter(V, F, barycenters);
    igl::local_basis(V, F, Bx, By, faceNormals);
  }
};

}



#endif /* DIRECTIONAL_TRIMESH_H */
