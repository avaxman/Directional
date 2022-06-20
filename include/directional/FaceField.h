// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_FACEFIELD_H
#define DIRECTIONAL_FACEFIELD_H

#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>

namespace directional{

class FaceField{
public:
  
  const TriMesh* mesh;
  
  int N;  //degree of field (how many vectors are in each point);
  
  Eigen::MatrixXd source;  //as barycentric coordinates in each face.
  Eigen::VectorXi face;
  
  Eigen::MatrixXd intField;  //the field in intrinsic representation (depending on the local basis of the face).
  Eigen::MatrixXd extField;  //the field in ambient coordinates.
  
  Eigen::MatrixXi adjSpaces;  //the adjacent tangent spaces
  //the connection between adjacent tangent space. That is, a field is parallel between adjaspaces(i,0) and adjSpaces(i,1) when complex(intField.row(adjSpaceS(i,0))*connection(i))=complex(intField.row(adjSpaceS(i,1))
  Eigen::MatrixXcd connection;
  
  Eigen::MatrixXd cycles;  //the adjaceny matrix of dual cycles
  Eigen::VectorXd cycleCurvatures;  //The Gaussian curvature of dual cycles.
  Eigen::VectorXi matching;
  Eigen::VectorXd effort;
  Eigen::VectorXi singVertices;
  Eigen::VectorXi singIndices;
  
  FaceField(){}
  ~FaceField(){}
  
  void IGL_INLINE set_field(const Eigen::MatrixXd& _extField,
                            const TriMesh& _mesh,
                            const Eigen::MatrixXd& _source=Eigen::MatrixXd(),
                            const Eigen::VectorXi& _face=Eigen::VectorXi()){
  
    extField=_extField;
    mesh = &_mesh;
    source = _source;
    face = _face;
    N = extField.cols()/3;
    intField.resize(extField.rows(),2*N);

    if (face.rows()==0){
      for (int i=0;i<N;i++)
        intField.block(0,2*i,intField.rows(),2)<<(extField.block(0,3*i,extField.rows(),3).array()*mesh->Bx.array()).rowwise().sum(),(extField.block(0,3*i,extField.rows(),3).array()*mesh->By.array()).rowwise().sum();
    } else {
      for (int i=0;i<N;i++)
        for (int j=0;j<face.size();j++)
          intField.block(0,2*i,1,2)<<(extField.block(0,3*i,1,3).array()*mesh->Bx.row(face(j)).array()).sum(),(extField.block(0,3*i,1,3).array()*mesh->By.row(face(j)).array()).sum();
    }
    
    adjSpaces = mesh->EF;
    //TODO: connection, cycles, cycleCurvature
  }
  
  void IGL_INLINE set_singularities(const Eigen::VectorXi& _singVertices,
                                    const Eigen::VectorXi& _singIndices){
    singVertices=_singVertices;
    singIndices=_singIndices;
  }
      
};

}



#endif /* DIRECTIONAL_TRIMESH_H */
