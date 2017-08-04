// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef POINT_SPHERES_H
#define POINT_SPHERES_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cmath> 


namespace directional
{
    // creates small spheres to visualize points on the overlay of the mesh
    // Inputs:
    //  P  eigen double matrix  #P by 3 - coordinates of the centers of spheres
    //  R  double               radii of the spheres
    //  C  eigen double matrix  #P by 3 - RBG colors per sphere
    //  res integer             the resolution of the sphere
    // colorPerVertex           speaks for its own
    // Outputs:
    //  V   eigen double matrix     spheres' meshes coordinates
    //  T   eigen int matrix        meshes triangles
    //  TC  eigen double matrix     #T by 3 colors (allocated from C)l
    IGL_INLINE bool point_spheres(const Eigen::MatrixXd& points,
                                  const double& radius,
                                  const Eigen::MatrixXd& colors,
                                  const int res,
                                  const bool colorPerVertex,
                                  Eigen::MatrixXd& V,
                                  Eigen::MatrixXi& T,
                                  Eigen::MatrixXd& TC)
    {
        using namespace Eigen;
        V.resize(res*res*points.rows(),3);
        T.resize(2*(res-1)*res*points.rows(),3);
        TC.resize((colorPerVertex ? V.rows() : T.rows()),3);
        
        for (int i=0;i<points.rows();i++){
            RowVector3d center=points.row(i);
            
            //creating vertices
            for (int j=0;j<res;j++){
                double z=center(2)+radius*cos(M_PI*(double)j/(double(res-1)));
                for (int k=0;k<res;k++){
                    double x=center(0)+radius*sin(M_PI*(double)j/(double(res-1)))*cos(2*M_PI*(double)k/(double(res-1)));
                    double y=center(1)+radius*sin(M_PI*(double)j/(double(res-1)))*sin(2*M_PI*(double)k/(double(res-1)));
                    V.row((res*res)*i+j*res+k)<<x,y,z;
                    if (colorPerVertex)
                        TC.row((res*res)*i+j*res+k)<<colors.row(i);
                }
            }
            
 
            //creating faces
            for (int j=0;j<res-1;j++){
                for (int k=0;k<res;k++){
                    int v1=(res*res)*i+j*res+k;
                    int v2=(res*res)*i+(j+1)*res+k;
                    int v3=(res*res)*i+(j+1)*res+(k+1)%res;
                    int v4=(res*res)*i+j*res+(k+1)%res;
                    T.row(2*(((res-1)*res)*i+res*j+k))<<v1,v2,v3;
                    T.row(2*(((res-1)*res)*i+res*j+k)+1)<<v4,v1,v3;
                    if (!colorPerVertex){
                        TC.row(2*(((res-1)*res)*i+res*j+k))<<colors.row(i);
                        TC.row(2*(((res-1)*res)*i+res*j+k)+1)<<colors.row(i);
                    }
                }
            }
        }
        
        return true;
    }
}

    


#endif


