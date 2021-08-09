// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


#ifndef gradient_h
#define gradient_h

namespace directional{

IGL_INLINE void gradient(const Eigen::MatrixXd& V,
                         const Eigen::MatrixXi& F,
                         const int N,
                         const Eigen::MatrixXd& cornerFunctions,
                         Eigen::MatrixXd& rawField)



{
  
  using namespace Eigen;
  using namespace std;
  
  VectorXd dblA;
  igl::doublearea(V,F,dblA);
  Eigen::MatrixXd normals;
  igl::per_face_normals(V, F, normals);
  rawField=MatrixXd::Zero(F.rows(), 3*N);
  for (int i=0;i<F.rows();i++){
    RowVector3d currNormal=normals.row(i);
    for (int k=0;k<N;k++){
      RowVector3d localGradient(0.0,0.0,0.0);
      for (int j=0;j<3;j++){
        RowVector3d eVec = V.row(F(i,(j+1)%3))-V.row(F(i,j));
        RowVector3d eVecRot = currNormal.cross(eVec);
        double oppCornerFunction = cornerFunctions(i, k+N*((j+2)%3));
        localGradient=localGradient+(oppCornerFunction/dblA(i))*eVecRot;
      }
      rawField.block(i, k*3, 1, 3)=localGradient;
    }
  }
}

}


#endif /* gradient_h */
