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

namespace Directional{

class FaceField{
public:
  
  TriMesh* tangentSpace;
  
  int N;  //degree of field (how many vectors are in each point);
  
  Eigen::MatrixXd source;  //as barycentric coordinates in each face.
  Eigen::VectorXi face;
  
  Eigen::MatrixXd intField;  //the field in intrinsic representation (depending on the local basis of the face).
  Eigen::MatrixXd extField;  //the field in ambient coordinates.
  
  FaceField(){}
  ~FaceField(){}
  
  void IGL_INLINE set_field(const Eigen::MatrixXd& _extField,
                            const Trimesh& mesh,
                            const Eigen::MatrixXd& source=Eigen::MatrixXd(),
                            const Eigen::MatrixXd& face=Eigen::MatrixXi()){
  
    extField=_extField;
    tangentSpace = &mesh;
    source = _source;
    face = _face;
    N = extField.cols()/3;
    intField.resize(extField.rows(),2*N);

    if (face.rows()==0){
      for (int i=0;i<N;i++)
        intField.block(0,2*i,intField.rows(),2)<<extField.block(0,3*i,extField.rows(),3).rowwise().dot(mesh.Bx.rowwise()),extField.block(0,3*i,extField.rows(),3).rowwise().dot(mesh.By.rowwise());
    } else {
      for (int i=0;i<N;i++)
        for (int j=0;j<face.size();j++)
          intField.block(0,2*i,1,2)<<extField.block(0,3*i,1,3).rowwise().dot(mesh.Bx.row(face(j))),extField.block(0,3*i,1,3).dot(mesh.By.row(face(j)));
    }
  }
      
};

}



#endif /* DIRECTIONAL_TRIMESH_H */
