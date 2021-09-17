// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_BAR_CHAIN_H
#define DIRECTIONAL_BAR_CHAIN_H

#include <Eigen/Core>
#include <string>
#include <vector>
#include <cmath>
#include <complex>
#include <igl/igl_inline.h>

namespace directional
  {
  // creates a mesh of conforming bar chains to visualize polylines (like function isolines) on meshes)
  // Input:
  //  isoV:          Each #P by 3 coordinates of the bar chains
  //  isoE:          #E by 2 edges of the bar
  //  isoN:          #E by 3 normals to the bars (corresponding to face normals at these points)
  //  width           bars widths
  //  heights:        bar heights (individual in #isoE)
  //  barColors:      #E by 3 RGB colors per bar
  // Output:
  //  V:              #V by 3 bar mesh coordinates
  //  T:              #T by 3 mesh triangles
  //  C:              #T by 3 colors
  
  bool IGL_INLINE bar_chains(const Eigen::MatrixXd& isoV,
                             const Eigen::MatrixXi& isoE,
                             const Eigen::MatrixXd& isoN,
                             const double& width,
                             const Eigen::MatrixXd& heights,
                             const Eigen::MatrixXd& barColors,
                             Eigen::MatrixXd& V,
                             Eigen::MatrixXi& T,
                             Eigen::MatrixXd& C,
                             const double margin = 0.1)
  {
    using namespace Eigen;
    
    //template mesh
    MatrixXd VBar(4,3);
    MatrixXi TBar(2,3);
    
    //Using some margins to offset the intersection problem with face sides - that's a hack until figuring out the edge trapezoid thing
    VBar<<0.0,0.0,0.0,
    0.0,-1.0,0.0,
    1.0,-1.0,0.0,
    1.0,0.0,0.0;
    
    TBar<<0,1,2,
    2,3,0;
    
    V.resize(VBar.rows()*isoE.rows(),3);
    T.resize(TBar.rows()*isoE.rows(),3);
    int NewColorSize=T.rows();
    C.resize(NewColorSize,3);
    
    for (int i=0;i<isoE.rows();i++){
      RowVector3d XAxis=(isoV.row(isoE(i,1))-isoV.row(isoE(i,0)));
      RowVector3d ZAxis=isoN.row(i);
      RowVector3d YAxis =ZAxis.cross(XAxis);
      YAxis.normalize();
      YAxis*=width;
      
      Matrix3d R; R<<(1.0+margin*2.0)*XAxis, YAxis, ZAxis;
      RowVector3d midway=YAxis/2.0-XAxis*margin;
      RowVector3d translation = isoV.row(isoE(i,0))+midway+heights(i)*isoN.row(i);
      
      V.block(VBar.rows()*i,0,VBar.rows(),3)=VBar*R+translation.replicate(VBar.rows(),1);
      T.block(TBar.rows()*i,0,TBar.rows(),3)=TBar.array()+VBar.rows()*i;
      C.block(TBar.rows()*i,0,TBar.rows(),3)=barColors.row(i).replicate(TBar.rows(),1);
    }
    
    return true;
  }
  

  }




#endif


