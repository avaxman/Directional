//This fiel is part of libdirectional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_LINE_CYLINDERS_H
#define DIRECTIONAL_LINE_CYLINDERS_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cmath> 
#include <complex>


namespace directional
{
  // creates small cylinders to visualize lines on the overlay of the mesh
  // Inputs:
  //  P1,P2  eigen double matrix  each #P by 3 - coordinates of the endpoints of the cylinders
  //  R  double               radii of the spheres
  //  C  eigen double matrix  #P by 3 - RBG colors per sphere
  //  res integer             the resolution of the sphere
  // colorPerVertex           speaks for its own
  // extend - if to extend the V,T,TC, or make them new
  // Outputs:
  //  V   eigen double matrix     spheres' meshes coordinates
  //  T   eigen int matrix        meshes triangles
  //  TC  eigen double matrix     #T by 3 colors (allocated from C)
  
  IGL_INLINE bool line_cylinders(const Eigen::MatrixXd& P1,
                                 const Eigen::MatrixXd& P2,
                                 const double& radius,
                                 const Eigen::MatrixXd& C,
                                 const int res,
                                 const bool colorPerVertex,
                                 const bool extend,
                                 Eigen::MatrixXd& V,
                                 Eigen::MatrixXi& T,
                                 Eigen::MatrixXd& TC)
  {
    using namespace Eigen;
    int VOffset, TOffset, TCOffset;
    if (!extend){
      V.resize(2*res*P1.rows(),3);
      T.resize(2*res*P1.rows(),3);
      int NewColorSize=(colorPerVertex ? V.rows() : T.rows());
      TC.resize(NewColorSize,3);
      VOffset=TOffset=TCOffset=0;
    } else {
      VOffset=V.rows();
      TOffset=T.rows();
      TCOffset=TC.rows();
      
      V.conservativeResize(VOffset+2*res*P1.rows(),3);
      T.conservativeResize(TOffset+2*res*P1.rows(),3);
      int NewColorSize=(colorPerVertex ? 2*res*P1.rows() : 2*res*P1.rows());
      TC.conservativeResize(TCOffset+NewColorSize,3);
      
    }
    RowVector3d ZAxis; ZAxis<<0.0,0.0,1.0;
    RowVector3d YAxis; YAxis<<0.0,1.0,0.0;
    
    MatrixXd PlanePattern(res,2);
    for (int i=0;i<res;i++){
      std::complex<double> CurrRoot=exp(2*M_PI*std::complex<double>(0,1)*(double)i/(double)res);
      PlanePattern.row(i)<<CurrRoot.real(), CurrRoot.imag();
    }
    
    
    for (int i=0;i<P1.rows();i++){
      RowVector3d NormAxis=(P2.row(i)-P1.row(i)).normalized();
      RowVector3d PlaneAxis1=NormAxis.cross(ZAxis);
      if (PlaneAxis1.norm()<10e-2)
        PlaneAxis1=NormAxis.cross(YAxis).normalized();
      else
        PlaneAxis1=PlaneAxis1.normalized();
      RowVector3d PlaneAxis2=NormAxis.cross(PlaneAxis1).normalized();
      for (int j=0;j<res;j++){
        int v1=2*res*i+2*j;
        int v2=2*res*i+2*j+1;
        int v3=2*res*i+2*((j+1)%res);
        int v4=2*res*i+2*((j+1)%res)+1;
        V.row(VOffset+v1)<<P1.row(i)+(PlaneAxis1*PlanePattern(j,0)+PlaneAxis2*PlanePattern(j,1))*radius;
        V.row(VOffset+v2)<<P2.row(i)+(PlaneAxis1*PlanePattern(j,0)+PlaneAxis2*PlanePattern(j,1))*radius;
        
        if (colorPerVertex){
          TC.row(TCOffset+v1)<<C.row(i);
          TC.row(TCOffset+v2)<<C.row(i);
        }
        
        
        T.row(TOffset+2*res*i+2*j)<<VOffset+v3,VOffset+v2,VOffset+v1;
        T.row(TOffset+2*res*i+2*j+1)<<VOffset+v4,VOffset+v2,VOffset+v3;
        
        if (!colorPerVertex){
          TC.row(TCOffset+2*res*i+2*j)<<C.row(i);
          TC.row(TCOffset+2*res*i+2*j+1)<<C.row(i);
        }
      }
    }
    return true;
  }
  
}




#endif


