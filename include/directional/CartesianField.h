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
#define POLYVECTOR_FIELD 1

class CartesianField{
public:
  
  const TriMesh* mesh;
  
  int N;  //degree of field (how many vectors are in each point);
  int fieldType;       //What the field actually represents  (for instance, either a raw field or a power/polyvector field)
  
  Eigen::MatrixXd source;  //as barycentric coordinates in each face.
  Eigen::VectorXi face;
  
  Eigen::MatrixXd intField;  //the field in intrinsic representation (depending on the local basis of the face).
  Eigen::MatrixXd extField;  //the field in ambient coordinates.
  
  Eigen::MatrixXi adjSpaces;  //the adjacent tangent spaces
  //the connection between adjacent tangent space. That is, a field is parallel between adjaspaces(i,0) and adjSpaces(i,1) when complex(intField.row(adjSpaceS(i,0))*connection(i))=complex(intField.row(adjSpaceS(i,1))
  Eigen::MatrixXcd connection;
  
  Eigen::SparseMatrix<double> dualCycles;  //the adjaceny matrix of dual cycles
  Eigen::VectorXd cycleCurvatures;  //The Gaussian curvature of dual cycles.
  Eigen::VectorXi matching;
  Eigen::VectorXd effort;
  Eigen::VectorXi singCycles;
  Eigen::VectorXi singIndices;
  Eigen::VectorXi innerAdjacencies;
  Eigen::VectorXi element2Cycle;
  
  CartesianField(){}
  CartesianField(const TriMesh& _mesh):mesh(&_mesh){}
  ~CartesianField(){}
  
  void virtual IGL_INLINE set_mesh(const TriMesh& _mesh){
    mesh = &_mesh;
  };
  
  void virtual IGL_INLINE set_extrinsic_field(const Eigen::MatrixXd& _extField,
                                              const Eigen::MatrixXd& _source=Eigen::MatrixXd(),
                                              const Eigen::VectorXi& _face=Eigen::VectorXi()){
    extField=_extField;
    source = _source;
    face = _face;
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
  
  void virtual IGL_INLINE set_singularities(const Eigen::VectorXi& _singCycles,
                                            const Eigen::VectorXi& _singIndices){}

};

}



#endif 
