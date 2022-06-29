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
#include <igl/barycenter.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/local_basis.h>
#include <igl/per_face_normals.h>
#include <igl/edge_topology.h>
#include <igl/boundary_loop.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/gaussian_curvature.h>

namespace directional{

class TriMesh{
public:
  
  //Basic quantities
  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  
  //Topological quantities
  Eigen::MatrixXi EF, FE, EV,TT, EFi;
  Eigen::MatrixXd FEs;
  Eigen::VectorXi innerEdges, boundEdges;
  
  //Geometric quantities
  Eigen::MatrixXd faceNormals;
  Eigen::MatrixXd Bx,By;  //local basis vectors per face
  Eigen::MatrixXd barycenters;
  Eigen::MatrixXd GaussianCurvature;
  int eulerChar;
  int numGenerators;
  
  std::vector<std::vector<int>> boundaryLoops;
  
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
    std::vector<int> innerEdgesList, boundEdgesList;
    for (int i=0;i<EF.rows();i++)
      if ((EF(i,1)==-1)||(EF(i,0)==-1))
        boundEdgesList.push_back(i);
      else
        innerEdgesList.push_back(i);
    
    
    //computing extra topological information
    //Relative location of edges within faces
    EFi = Eigen::MatrixXi::Constant(EF.rows(), 2, -1); // number of an edge inside the face
    for(int i = 0; i < EF.rows(); i++)
    {
      for (int k = 0; k < 2; k++)
      {
        if (EF(i, k) == -1)
          continue;
        for (int j = 0; j < 3; j++)
          if (FE(EF(i, k), j) == i)
            EFi(i, k) = j;
      }
    }
    
    //sign of edge within face
    FEs = Eigen::MatrixXd::Zero(FE.rows(), FE.cols());
    
    for(int i = 0; i < EF.rows(); i++)
    {
      if(EFi(i, 0) != -1)
        FEs(EF(i, 0), EFi(i, 0)) = 1.0;
      if(EFi(i,1) != -1)
        FEs(EF(i, 1), EFi(i, 1)) = -1.0;
    }
    
    innerEdges = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(innerEdgesList.data(), innerEdgesList.size());
    boundEdges = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(boundEdgesList.data(), boundEdgesList.size());
    igl::barycenter(V, F, barycenters);
    igl::local_basis(V, F, Bx, By, faceNormals);
    igl::triangle_triangle_adjacency(F, TT);
    igl::boundary_loop(F, boundaryLoops);
    igl::gaussian_curvature(V,F,GaussianCurvature);
    eulerChar = V.rows() - EV.rows() + F.rows();
    numGenerators = (2 - eulerChar)/2 - boundaryLoops.size();
  }

};

}



#endif /* DIRECTIONAL_TRIMESH_H */
