// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2020 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_ANGLED_ARROWS_H
#define DIRECTIONAL_ANGLED_ARROWS_H

#include <Eigen/Core>
#include <string>
#include <vector>
#include <cmath> 
#include <complex>
#include <igl/igl_inline.h>

namespace directional
  {
  // creates a mesh of asymmetric rhombus prisms to visualize vector fields
  // Input:
  //  P1,P2:          Each #P by 3 coordinates of the box endpoints
  //  normals:        Normals to the arrows (w.r.t. height).
  //  length, width, height, angle:  angle dimensions
  //  angleColors:      #P by 3 RBG colors per box
  // Output:
  //  V:              #V by 3 arrow mesh coordinates
  //  T:              #T by 3 mesh triangles
  //  C:              #T by 3 colors
  
  IGL_INLINE bool angled_arrows(const Eigen::MatrixXd& P1,
                               const Eigen::MatrixXd& P2,
                               const Eigen::MatrixXd& normals,
                               const double& width,
                               const double& height,
                               const double& angle,
                               const Eigen::MatrixXd& arrowColors,
                               Eigen::MatrixXd& V,
                               Eigen::MatrixXi& T,
                               Eigen::MatrixXd& C)
  {
    using namespace Eigen;
    
    //template mesh
    MatrixXd VArrow(4,3);
    MatrixXi TArrow(2,3);
    
    VArrow<<0.0,0.0,0.0,
    (width/2.0)*(cos(angle/2.0)/sin(angle/2.0)), -width/2.0, 0.0,
    1.0,0.0,0.0,
    (width/2.0)*(cos(angle/2.0)/sin(angle/2.0)), width/2.0, 0.0;
    
    TArrow<<0,1,2,
    2,3,0;
    
    V.resize(VArrow.rows()*P1.rows(),3);
    T.resize(TArrow.rows()*P1.rows(),3);
    int NewColorSize=T.rows();
    C.resize(NewColorSize,3);
    
    for (int i=0;i<P1.rows();i++){
      RowVector3d YAxis=(P2.row(i)-P1.row(i));
      RowVector3d ZAxis=normals.row(i);
      RowVector3d XAxis =YAxis.cross(ZAxis);
      XAxis.rowwise().normalize();
      XAxis*=YAxis.norm();

      Matrix3d R; R<<XAxis, YAxis, ZAxis;
      RowVector3d translation = P1.row(i)+height*normals.row(i);
      
      V.block(VArrow.rows()*i,0,VArrow.rows(),3)=VArrow*R+translation.replicate(VArrow.rows(),1);
      T.block(TArrow.rows()*i,0,TArrow.rows(),3)=TArrow.array()+VArrow.rows()*i;
      C.block(TArrow.rows()*i,0,TArrow.rows(),3)=arrowColors.row(i).replicate(TArrow.rows(),1);
    }
    
    return true;
  }
  
  }




#endif


