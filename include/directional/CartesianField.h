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

//This is the interface class for any directional fields represented in cartesian coordinates, of any order N.
class CartesianField{
public:
  
  const TriMesh* mesh;
  
  int N;  //degree of field (how many vectors are in each point);
  int fieldType;       //What the field actually represents  (for instance, either a raw field or a power/polyvector field)
  
  //In case some methods are only defined for several classes
  virtual discTangTypeEnum discTangType() const {return discTangTypeEnum::BASE_CLASS;}
  
  Eigen::MatrixXd sources;  //the source point of the extrinsic vectors
  Eigen::MatrixXd normals;  //the normals to the tangent spaces
  Eigen::MatrixXd dualSources;  //source point of dual cycles
  Eigen::MatrixXd dualNormals;  //normals to dual cycles
  
  //TODO intField being a single vector in case of power field
  Eigen::MatrixXd intField;  //the field in intrinsic representation (depending on the local basis of the face). Size #T x 2N
  Eigen::MatrixXd extField;  //the field in ambient coordinates. Size Size #T x 3N
  
  
  Eigen::MatrixXi adjSpaces;  //the adjacent tangent spaces (tangent bundle edges)
  Eigen::MatrixXi oneRing;    //the tangent spaces adjacent around each tangent space
  
  //tangent space geometry
  //the connection between adjacent tangent space. That is, a field is parallel between adjaspaces(i,0) and adjSpaces(i,1) when complex(intField.row(adjSpaceS(i,0))*connection(i))=complex(intField.row(adjSpaceS(i,1))
  Eigen::VectorXcd connection;
  Eigen::VectorXd stiffnessWeights;   //the mass matrix of the inner norm between directional fields
  Eigen::VectorXd massWeights;        //the mass matrix (like "area") of each tangent space itself
  
  Eigen::SparseMatrix<double> dualCycles;  //the adjaceny matrix of dual cycles
  Eigen::VectorXd cycleCurvatures;  //The Gaussian curvature of dual cycles.
  Eigen::VectorXi matching;         //matching(i)=j when vector k in adjSpaces(i,0) matches to vector (k+j)%N in adjSpaces(i,1)
  Eigen::VectorXd effort;           //The effort of the entire matching (sum of deviations from parallel transport)
  Eigen::VectorXi singElements;     //The singular (dual elements). Only the local cycles! not the generators or boundary cycles
  Eigen::VectorXi singIndices;      //Corresponding indices (this is the numerator where the true fractional index is singIndices/N);
  Eigen::VectorXi innerAdjacencies; //the indices into adjSpaces that are not boundary
  Eigen::VectorXi element2Cycle;    //a map between dual elements and dual cycles
  
  CartesianField(){}
  CartesianField(const TriMesh& _mesh):mesh(&_mesh){}
  ~CartesianField(){}
  
  //Initializing the field with the proper tangent spaces
  void virtual IGL_INLINE init_field(const TriMesh& _mesh, const int _fieldType, const int _N){
    mesh = &_mesh;
    fieldType = _fieldType;
    N=_N;
  };
  
  //Setting the field by the extrinsic ambient field, which will get projected to the intrinsic tangent spaces.
  void virtual IGL_INLINE set_extrinsic_field(const Eigen::MatrixXd& _extField){
    extField=_extField;
    intField = extField;  //the basic version
    N = extField.cols()/3;
  }
  
  //Directly setting the intrinsic tangent spaces
  void virtual IGL_INLINE set_intrinsic_field(const Eigen::MatrixXd& _intField){
    intField = _intField;
    extField = intField;
  }
  
  //The same, just with complex coordinates
  void virtual IGL_INLINE set_intrinsic_field(const Eigen::MatrixXcd& _intField){
    intField.resize(_intField.rows(),_intField.cols()*2);
    for (int i=0;i<N;i++){
      intField.col(2*i)=_intField.col(i).real();
      intField.col(2*i+1)=_intField.col(i).imag();
    }
    extField = intField;
  }
  
  //projecting an arbitrary set of extrinsic vectors (e.g. coming from user-prescribed constraints) into intrinsic vectors.
  Eigen::MatrixXd virtual IGL_INLINE project_to_intrinsic(const Eigen::VectorXi& tangentSpaces, const Eigen::MatrixXd& extDirectionals) const {
    return extDirectionals;
  }
  
  //Directly setting the singularities of the the field (only at the local dual elements; not at generator or boundary cycles).
  void virtual IGL_INLINE set_singularities(const Eigen::VectorXi& _singElements,
                                            const Eigen::VectorXi& _singIndices){}
  

};

}



#endif 
