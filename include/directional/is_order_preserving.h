// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


#ifndef is_order_preserving_h
#define is_order_preserving_h

namespace directional{

IGL_INLINE void is_order_preserving(const Eigen::MatrixXd& V,
                             const Eigen::MatrixXi& F,
                             const Eigen::MatrixXd& rawField,
                             Eigen::VectorXi& isOrderPreserving)
{
  Eigen::MatrixXd normals;
  igl::per_face_normals(V,F, normals);
  
  int N = rawField.cols()/3;
  isOrderPreserving=Eigen::VectorXi::Constant(rawField.rows(),1);
  for (int i=0;i<rawField.rows();i++){
    for (int j=0;j<N;j++){
      Eigen::RowVector3d v1=rawField.block(i,3*j,1,3);
      Eigen::RowVector3d v2=rawField.block(i,3*((j+1)%N),1,3);
      if (normals.row(i).dot(v1.cross(v2))<=0.0){
        isOrderPreserving(i)=0;
        continue;
      }
    }
  }
}

}


#endif /* is_order_preserving_h */
