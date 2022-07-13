// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_TANGENT_BUNDLE_H
#define DIRECTIONAL_TANGENT_BUNDLE_H

#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <directional/TriMesh.h>

namespace directional{

#define RAW_FIELD 0
#define POWER_FIELD 1
#define POLYVECTOR_FIELD 2

enum class discTangTypeEnum {BASE_CLASS, FACE_SPACES, VERTEX_SPACES};

//This is the interface class for any directional fields represented in cartesian coordinates, of any order N.
class TangentBundle{
public:
    
  //In case some methods are only defined for several classes
  virtual discTangTypeEnum discTangType() const {return discTangTypeEnum::BASE_CLASS;}
  
  Eigen::MatrixXd sources;  //the source point of the extrinsic vectors
  Eigen::MatrixXd normals;  //the normals to the tangent spaces
  Eigen::MatrixXd cycleSources;  //source point of cycles
  Eigen::MatrixXd cycleNormals;  //normals to cycles
  
  Eigen::MatrixXi adjSpaces;  //the adjacent tangent spaces (tangent bundle edges)
  Eigen::MatrixXi oneRing;    //the tangent spaces adjacent around each tangent space
  
  //tangent space geometry
  //the connection between adjacent tangent space. That is, a field is parallel between adjaspaces(i,0) and adjSpaces(i,1) when complex(intField.row(adjSpaceS(i,0))*connection(i))=complex(intField.row(adjSpaceS(i,1))
  Eigen::VectorXcd connection;
  Eigen::VectorXd stiffnessWeights;   //the mass matrix of the inner norm between directional fields
  Eigen::VectorXd massWeights;        //the mass matrix (like "area") of each tangent space
  
  Eigen::SparseMatrix<double> cycles;  //the adjaceny matrix of cycles
  Eigen::VectorXd cycleCurvatures;  //The Gaussian curvature of dual cycles.
  
  Eigen::VectorXi innerAdjacencies; //the indices into adjSpaces that are not boundary
  Eigen::VectorXi local2Cycle;    //a map between local cycles and general cycles
  
  TangentBundle(){}
  ~TangentBundle(){}
  
  //projecting an arbitrary set of extrinsic vectors (e.g. coming from user-prescribed constraints) into intrinsic vectors.
  Eigen::MatrixXd virtual IGL_INLINE project_to_intrinsic(const Eigen::VectorXi& tangentSpaces, const Eigen::MatrixXd& extDirectionals) const {
    return extDirectionals;
  }
  
  //projecting extrinsic to intrinsic
  Eigen::MatrixXd virtual IGL_INLINE project_to_extrinsic(const Eigen::VectorXi& tangentSpaces, const Eigen::MatrixXd& extDirectionals) const {
    return extDirectionals;
  }
  
};

}



#endif
