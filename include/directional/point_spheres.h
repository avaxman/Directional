// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POINT_SPHERES_H
#define DIRECTIONAL_POINT_SPHERES_H

#include <string>
#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <igl/igl_inline.h>


namespace directional
{
  // creates small spheres to visualize P on the overlay of the mesh
  // Input:
  //  P:      #P by 3 coordinates of the centers of spheres
  //  N:      #P by 3 normals (the south-north pole direction of the spheres).
  //  r: radii of the spheres
  //  sphereColors:      #P by 3 - RBG colors per sphere
  //  res:    the resolution of the sphere discretization
  // extendMesh if to extend the V,T,TC, or to overwrite them
  // Output:
  //  V:    #V by 3 sphere mesh coordinates
  //  T     #T by 3 sphere mesh triangles
  //  C:    #T by 3 vertex-based colors
  IGL_INLINE bool point_spheres(const Eigen::MatrixXd& P,
                                const Eigen::MatrixXd& normals,
                                const double& r,
                                const Eigen::MatrixXd& sphereColors,
                                const int res,
                                Eigen::MatrixXd& V,
                                Eigen::MatrixXi& T,
                                Eigen::MatrixXd& C)
  {
    using namespace Eigen;
    /*V.resize(res*res*P.rows(),3);
    T.resize(2*(res-1)*res*P.rows(),3);
    C.resize(V.rows(),3);*/
    
    MatrixXd VSphere(res*res,3);
    MatrixXi TSphere(2*(res-1)*res,3);
    
    //creating template sphere vertices
    for (int j=0;j<res;j++){
      double z=r*cos(igl::PI*(double)j/(double(res-1)));
      for (int k=0;k<res;k++){
        double x=r*sin(igl::PI*(double)j/(double(res-1)))*cos(2*igl::PI*(double)k/(double(res)));
        double y=r*sin(igl::PI*(double)j/(double(res-1)))*sin(2*igl::PI*(double)k/(double(res)));
        VSphere.row(j*res+k)<<x,y,z;
      }
    }
    
  
    for (int j=0;j<res-1;j++){
      for (int k=0;k<res;k++){
        int v1=j*res+k;
        int v2=(j+1)*res+k;
        int v3=(j+1)*res+(k+1)%res;
        int v4=j*res+(k+1)%res;
        TSphere.row(2*(res*j+k))<<v1,v2,v3;
        TSphere.row(2*(res*j+k)+1)<<v4,v1,v3;
      }
    }
    
    //std::cout<<"TSphere: "<<TSphere<<std::endl;
    V.resize(VSphere.rows()*P.rows(),3);
    T.resize(TSphere.rows()*P.rows(),3);
    C.resize(V.rows(),3);
    
    for (int i=0;i<P.rows();i++){
      RowVector3d ZAxis=normals.row(i);
      ZAxis.normalize();
      RowVector3d XAxis; XAxis<<0.0, -normals(i,2), normals(i,1);
      if (XAxis.squaredNorm()<1e-4)
        XAxis<<-normals(i,2),0.0,normals(i,0);
      XAxis.normalize();
      
      RowVector3d YAxis =ZAxis.cross(XAxis);
      YAxis.rowwise().normalize();
    
      Matrix3d R; R<<XAxis, YAxis, ZAxis;
      RowVector3d translation = P.row(i);
      
      V.block(VSphere.rows()*i,0,VSphere.rows(),3)=VSphere*R+translation.replicate(VSphere.rows(),1);
      T.block(TSphere.rows()*i,0,TSphere.rows(),3)=TSphere.array()+VSphere.rows()*i;
      C.block(VSphere.rows()*i,0,VSphere.rows(),3)=sphereColors.row(i).replicate(VSphere.rows(),1);
    }

    return true;
  }

}




#endif


