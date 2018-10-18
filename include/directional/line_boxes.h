// This file is part of libdirectional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_LINE_BOXES_H
#define DIRECTIONAL_LINE_BOXES_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cmath> 
#include <complex>


namespace directional
{
  // creates a mesh of thin boxes to visualize lines on the overlay of the mesh
  // Inputs:
  //  P1,P2         each #P by 3 coordinates of the box endpoints
  //  normals       normals to the boxes (w.r.t. height).
  //  width, height box dimensions
  //  cyndColors    #P by 3 RBG colors per cylinder
  //  res           the resolution of the cylinder (size of base polygon)
  // colorPerVertex in the output mesh
  // extendMesh     if to extend the V,T,TC, or to overwrite them
  // Outputs:
  //  V             #V by 3 cylinder mesh coordinates
  //  T             #T by 3 mesh triangles
  //  C             #T/#V by 3 colors
  
  IGL_INLINE bool line_boxes(const Eigen::MatrixXd& P1,
                             const Eigen::MatrixXd& P2,
                             const Eigen::MatrixXd& normals,
                             const double& width,
                             const double& height,
                             const Eigen::MatrixXd& boxColors,
                             const bool colorPerVertex,
                             const bool extendMesh,
                             Eigen::MatrixXd& V,
                             Eigen::MatrixXi& T,
                             Eigen::MatrixXd& C)
  {
    using namespace Eigen;
    
    //template mesh
    MatrixXd VBox(8,3);
    MatrixXi TBox(12,3);
    
    VBox<<0.0, 0.0, 0.0,
    1.0, 0.0, 0.0,
    1.0, 0.0, 1.0,
    0.0, 0.0, 1.0,
    0.0, 1.0, 0.0,
    1.0, 1.0, 0.0,
    1.0, 1.0, 1.0,
    0.0, 1.0, 1.0,
    
    TBox<<0,1,2,
    2,3,1,
    2,1,5,
    5,6,2,
    7,6,5,
    5,4,7,
    4,5,1,
    1,0,4,
    3,7,4,
    4,0,3,
    3,2,6,
    6,7,3;
    
    int VOffset, TOffset, COffset;
    if (!extendMesh){
      V.resize(VBox.rows()*P1.rows(),3);
      T.resize(TBox.rows()*P1.rows(),3);
      int NewColorSize=(colorPerVertex ? V.rows() : T.rows());
      C.resize(NewColorSize,3);
      VOffset=TOffset=COffset=0;
    } else {
      VOffset=V.rows();
      TOffset=T.rows();
      COffset=C.rows();
      
      V.conservativeResize(VOffset+VBox.rows()*P1.rows(),3);
      T.conservativeResize(TOffset+TBox.rows()*P1.rows(),3);
      int NewColorSize=(colorPerVertex ? VBox.rows()*P1.rows() : TBox.rows()*P1.rows());
      C.conservativeResize(COffset+NewColorSize,3);
      
    }
   
    for (int i=0;i<P1.rows();i++){
      RowVector3d YAxis=(P2.row(i)-P1.row(i));
      RowVector3d ZAxis=normals.row(i);
      RowVector3d XAxis =YAxis.cross(ZAxis);
      XAxis.rowwise().normalize();
      XAxis*=width;
      ZAxis*=height;
      
      Matrix3d R; R<<XAxis, YAxis, ZAxis;
      
      RowVector3d P1onBox; P1onBox<<0.5, 0, 0.5;
      RowVector3d translation = P1.row(i) - P1onBox*R;
     
      V.block(VOffset+VBox.rows()*i,0,VBox.rows(),3)=VBox*R+translation.replicate(VBox.rows(),1);
      T.block(TOffset+TBox.rows()*i,0,TBox.rows(),3)=TBox.array()+VOffset+VBox.rows()*i;
      if (colorPerVertex)
        C.block(COffset+VBox.rows()*i,0,VBox.rows(),3)=boxColors.row(i).replicate(VBox.rows(),1);
      else
        C.block(TOffset+TBox.rows()*i,0,TBox.rows(),3)=boxColors.row(i).replicate(TBox.rows(),1);
    }
    
    
      
      
      /*for (int j=0;j<res;j++){
        int v1=2*res*i+2*j;
        int v2=2*res*i+2*j+1;
        int v3=2*res*i+2*((j+1)%res);
        int v4=2*res*i+2*((j+1)%res)+1;
        V.row(VOffset+v1)<<P1.row(i)+(PlaneAxis1*PlanePattern(j,0)+PlaneAxis2*PlanePattern(j,1))*radius;
        V.row(VOffset+v2)<<P2.row(i)+(PlaneAxis1*PlanePattern(j,0)+PlaneAxis2*PlanePattern(j,1))*radius;
        
        if (colorPerVertex){
          C.row(COffset+v1)<<cyndColors.row(i);
          C.row(COffset+v2)<<cyndColors.row(i);
        }
        
        
        T.row(TOffset+2*res*i+2*j)<<VOffset+v3,VOffset+v2,VOffset+v1;
        T.row(TOffset+2*res*i+2*j+1)<<VOffset+v4,VOffset+v2,VOffset+v3;
        
        if (!colorPerVertex){
          C.row(COffset+2*res*i+2*j)<<cyndColors.row(i);
          C.row(COffset+2*res*i+2*j+1)<<cyndColors.row(i);
        }
      }
    }*/
    return true;
  }
  
}




#endif


