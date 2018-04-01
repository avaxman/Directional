// This file is part of libdirectional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POINT_SPHERES_H
#define DIRECTIONAL_POINT_SPHERES_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cmath> 


namespace directional
{
  // creates small spheres to visualize points on the overlay of the mesh
  // Inputs:
  //  P       #P by 3 coordinates of the centers of spheres
  //  radius  double radii of the spheres
  //  C       #P by 3 - RBG colors per sphere
  //  res integer   the resolution of the sphere
  // colorPerVertex of the output mesh
  // extendMesh if to extend the V,T,TC, or to overwrite them
  // Outputs:
  //  V   #V by 3 cylinder mesh coordinates
  //  T   #T by 3 mesh triangles
  //  TC  #T by 3 colors per triangle
  IGL_INLINE bool point_spheres(const Eigen::MatrixXd& points,
                                const double& radius,
                                const Eigen::MatrixXd& colors,
                                const int res,
                                const bool colorPerVertex,
                                const bool extendMesh,
                                Eigen::MatrixXd& V,
                                Eigen::MatrixXi& T,
                                Eigen::MatrixXd& TC)
  {
    using namespace Eigen;
    int VOffset, TOffset, TCOffset;
    if (!extendMesh){
      V.resize(res*res*points.rows(),3);
      T.resize(2*(res-1)*res*points.rows(),3);
      TC.resize((colorPerVertex ? V.rows() : T.rows()),3);
      VOffset=TOffset=TCOffset=0;
    } else{
      VOffset=V.rows();
      TOffset=T.rows();
      TCOffset=TC.rows();
      V.conservativeResize(V.rows()+res*res*points.rows(),3);
      T.conservativeResize(T.rows()+2*(res-1)*res*points.rows(),3);
      TC.conservativeResize(TC.rows()+(colorPerVertex ? res*res*points.rows() : 2*(res-1)*res*points.rows()),3);
    }
    
    for (int i=0;i<points.rows();i++){
      RowVector3d center=points.row(i);
      
      //creating vertices
      for (int j=0;j<res;j++){
        double z=center(2)+radius*cos(M_PI*(double)j/(double(res-1)));
        for (int k=0;k<res;k++){
          double x=center(0)+radius*sin(M_PI*(double)j/(double(res-1)))*cos(2*M_PI*(double)k/(double(res-1)));
          double y=center(1)+radius*sin(M_PI*(double)j/(double(res-1)))*sin(2*M_PI*(double)k/(double(res-1)));
          V.row(VOffset+(res*res)*i+j*res+k)<<x,y,z;
          if (colorPerVertex)
            TC.row(TCOffset+(res*res)*i+j*res+k)<<colors.row(i);
        }
      }
      
      
      //creating faces
      for (int j=0;j<res-1;j++){
        for (int k=0;k<res;k++){
          int v1=VOffset+(res*res)*i+j*res+k;
          int v2=VOffset+(res*res)*i+(j+1)*res+k;
          int v3=VOffset+(res*res)*i+(j+1)*res+(k+1)%res;
          int v4=VOffset+(res*res)*i+j*res+(k+1)%res;
          T.row(TOffset+2*(((res-1)*res)*i+res*j+k))<<v1,v2,v3;
          T.row(TOffset+2*(((res-1)*res)*i+res*j+k)+1)<<v4,v1,v3;
          if (!colorPerVertex){
            TC.row(TCOffset+2*(((res-1)*res)*i+res*j+k))<<colors.row(i);
            TC.row(TCOffset+2*(((res-1)*res)*i+res*j+k)+1)<<colors.row(i);
          }
        }
      }
    }
    
    return true;
  }
}




#endif


