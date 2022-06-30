// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_CARTESIAN_FIELD_H
#define DIRECTIONAL_CARTESIAN_FIELD_H

#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <directional/TriMesh.h>

namespace directional{

#define RAW_FIELD 0
#define POWER_FIELD 1
#define POLYVECTOR_FIELD 2

enum class discTangTypeEnum {BASE_CLASS, FACE_SPACES, VERTEX_SPACES};

class CartesianField{
public:
  
  const TriMesh* mesh;
  
  int N;  //degree of field (how many vectors are in each point);
  int fieldType;       //What the field actually represents  (for instance, either a raw field or a power/polyvector field)
  
  //In case some methods are only defined for several classes
  virtual static discTangTypeEnum discTangType(){return BASE_CLASS;}
  
  Eigen::MatrixXd sources;  //the source point of the extrinsic vectors
  Eigen::MatrixXd normals;  //the normals to the tangent spaces
  Eigen::VectorXi face;
  
  //TODO intField being a single vector in case of power field
  Eigen::MatrixXd intField;  //the field in intrinsic representation (depending on the local basis of the face).
  Eigen::MatrixXd extField;  //the field in ambient coordinates.
  
  
  Eigen::MatrixXi adjSpaces;  //the adjacent tangent spaces (tangent bundle edges)
  Eigen::MatrixXi oneRing;    //the tangent spaces adjacent around each tangent space
  
  //tangent space geometry
  //the connection between adjacent tangent space. That is, a field is parallel between adjaspaces(i,0) and adjSpaces(i,1) when complex(intField.row(adjSpaceS(i,0))*connection(i))=complex(intField.row(adjSpaceS(i,1))
  Eigen::VectorXcd connection;
  Eigen::VectorXd stiffnessWeights;
  Eigen::VectorXd massWeights;
  
  Eigen::SparseMatrix<double> dualCycles;  //the adjaceny matrix of dual cycles
  Eigen::VectorXd cycleCurvatures;  //The Gaussian curvature of dual cycles.
  Eigen::VectorXi matching;
  Eigen::VectorXd effort;
  Eigen::VectorXi singElements;
  Eigen::VectorXi singIndices;
  Eigen::VectorXi innerAdjacencies;
  Eigen::VectorXi element2Cycle;
  
  CartesianField(){}
  CartesianField(const TriMesh& _mesh):mesh(&_mesh){}
  ~CartesianField(){}
  
  void virtual IGL_INLINE init_field(const TriMesh& _mesh, const int _fieldType, const int _N){
    mesh = &_mesh;
    fieldType = _fieldType;
    N=_N;
  };
  
  void virtual IGL_INLINE set_extrinsic_field(const Eigen::MatrixXd& _extField){
    extField=_extField;
    intField = extField;  //the basic version
    N = extField.cols()/3;
    

    if (face.rows()==0){
      intField=extField;
    } else {
      intField.resize(face.rows(),3*N);
      for (int i=0;i<N;i++)
        for (int j=0;j<face.size();j++)
          intField.block(0,3*i,1,3)=extField.block(0,3*face(i),1,3);
    }
  }
  void virtual IGL_INLINE set_intrinsic_field(const Eigen::MatrixXd& _intField){
    intField = _intField;
    extField = intField;
  }
  
  void virtual IGL_INLINE set_intrinsic_field(const Eigen::MatrixXcd& _intField){
    intField.resize(_intField.rows(),_intField.cols()*2);
    for (int i=0;i<N;i++){
      intField.col(2*i)=_intField.col(i).real();
      intField.col(2*i+1)=_intField.col(i).imag();
    }
    extField = intField;
  }
  
  Eigen::MatrixXd virtual IGL_INLINE project_to_intrinsic(const Eigen::VectorXi& tangentSpaces, const Eigen::MatrixXd& extDirectionals) const {
    return extDirectionals;
  }
  
  void virtual IGL_INLINE set_singularities(const Eigen::VectorXi& _singElements,
                                            const Eigen::VectorXi& _singIndices){}
  

};

}



#endif 
