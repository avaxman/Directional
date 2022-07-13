// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_INTRINSIC_FACE_TANGENT_BUNDLE_H
#define DIRECTIONAL_INTRINSIC_FACE_TANGENT_BUNDLE_H

#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <igl/boundary_loop.h>
#include <igl/doublearea.h>
#include <directional/dual_cycles.h>
#include <directional/TriMesh.h>


namespace directional{

//This is the interface class for any directional fields represented in cartesian coordinates, of any order N.
class IntrinsicFaceTangentBundle : public TangentBundle{
public:
  
  const TriMesh* mesh;
  Eigen::MatrixXd tangentStartAngles;  //where each edge begins on the intrinsic space
  
  virtual discTangTypeEnum discTangType() const {return discTangTypeEnum::VERTEX_SPACES;}
  
  IntrinsicFaceTangentBundle(){}
  ~IntrinsicFaceTangentBundle(){}
  
  void IGL_INLINE init(const TriMesh& _mesh){
    
    typedef std::complex<double> Complex;
    mesh = &_mesh;
    
    //adjacency relation is by dual edges.
    adjSpaces = mesh->EF;
    oneRing = mesh->FE;
    sources = mesh->barycenters;
    normals = mesh->faceNormals;
    cycleSources = mesh->V;
    cycleNormals = mesh->vertexNormals;
    
    //connection is the ratio of the complex representation of edges
    connection.resize(mesh->EF.rows(),1);  //the difference in the angle representation of edge i from EF(i,0) to EF(i,1)
    Eigen::MatrixXd edgeVectors(mesh->EF.rows(), 3);
    for (int i = 0; i < mesh->EF.rows(); i++) {
      if (mesh->EF(i, 0) == -1 || mesh->EF(i, 1) == -1)
        continue;
      edgeVectors.row(i) = (mesh->V.row(mesh->EV(i, 1)) - mesh->V.row(mesh->EV(i, 0))).normalized();
      Complex ef(edgeVectors.row(i).dot(mesh->FBx.row(mesh->EF(i, 0))), edgeVectors.row(i).dot(mesh->FBy.row(mesh->EF(i, 0))));
      Complex eg(edgeVectors.row(i).dot(mesh->FBx.row(mesh->EF(i, 1))), edgeVectors.row(i).dot(mesh->FBy.row(mesh->EF(i, 1))));
      connection(i) = eg / ef;
    }
    
    //TODO: cycles, cycleCurvature
    directional::dual_cycles(mesh->V, mesh->F, mesh->EV, mesh->EF, cycles, cycleCurvatures, local2Cycle, innerAdjacencies);
    
    //drawing from mesh geometry
    
    /************Smoothness matrices****************/
    stiffnessWeights=Eigen::VectorXd::Zero(mesh->EF.rows());
   
    //mass are face areas
    igl::doublearea(mesh->V,mesh->F,massWeights);
    massWeights.array()/=2.0;

    //The "harmonic" weights from [Brandt et al. 2020].
    for (int i=0;i<mesh->EF.rows();i++){
      if ((mesh->EF(i,0)==-1)||(mesh->EF(i,1)==-1))
        continue;  //boundary edge
      
      double primalLengthSquared = (mesh->V.row(mesh->EV(i,0))-mesh->V.row(mesh->EV(i,1))).squaredNorm();
      stiffnessWeights(i)=3*primalLengthSquared/(massWeights(mesh->EF(i,0))+massWeights(mesh->EF(i,0)));
    }
  
  }
  
  
  //projecting an arbitrary set of extrinsic vectors (e.g. coming from user-prescribed constraints) into intrinsic vectors.
  Eigen::MatrixXd  virtual IGL_INLINE project_to_intrinsic(const Eigen::VectorXi& tangentSpaces, const Eigen::MatrixXd& extDirectionals) const{
    assert(tangentSpaces.rows()==extDirectionals.rows() || tangentSpaces.rows()==0);
    
    Eigen::VectorXi actualTangentSpaces;
    if (tangentSpaces.rows()==0)
      actualTangentSpaces = Eigen::VectorXi::LinSpaced(sources.rows(), 0, sources.rows()-1);
    else
      actualTangentSpaces = tangentSpaces;
    
    int N = extDirectionals.cols()/3;
    Eigen::MatrixXd intDirectionals(actualTangentSpaces.rows(),2*N);
    
    for (int i=0;i<actualTangentSpaces.rows();i++)
      for (int j=0;j<N;j++)
      intDirectionals.block(i,2*j,1,2)<<(extDirectionals.block(i,3*j,1,3).array()*mesh->FBx.row(actualTangentSpaces(i)).array()).sum(),(extDirectionals.block(i,3*j,1,3).array()*mesh->FBy.row(actualTangentSpaces(i)).array()).sum();
    
    return intDirectionals;
  }

  
  //projecting intrinsic to extrinsicc
  Eigen::MatrixXd virtual IGL_INLINE project_to_extrinsic(const Eigen::VectorXi& tangentSpaces, const Eigen::MatrixXd& intDirectionals) const {
    
    assert(tangentSpaces.rows()==intDirectionals.rows() || tangentSpaces.rows()==0);
    Eigen::VectorXi actualTangentSpaces;
    if (tangentSpaces.rows()==0)
      actualTangentSpaces = Eigen::VectorXi::LinSpaced(sources.rows(), 0, sources.rows()-1);
    else
      actualTangentSpaces = tangentSpaces;
    
    int N = intDirectionals.cols()/2;
    Eigen::MatrixXd extDirectionals(actualTangentSpaces.rows(),3);
    
    extDirectionals.conservativeResize(intDirectionals.rows(),intDirectionals.cols()*3/2);
    for (int i=0;i<intDirectionals.rows();i++)
      for (int j=0;j<intDirectionals.cols();j+=2)
        extDirectionals.block(i,3*j/2,1,3)=mesh->FBx.row(actualTangentSpaces(i))*intDirectionals(i,j)+mesh->FBy.row(actualTangentSpaces(i))*intDirectionals(i,j+1);
    
    return extDirectionals;
  }
  
};

}



#endif
  

