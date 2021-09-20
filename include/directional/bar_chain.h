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
#include <igl/doublearea.h>

namespace directional
  {
  // creates a mesh of conforming bar chains to visualize polylines (like function isolines) on meshes)
  // Input:
  //  V, F:           The original mesh
  //  isoV:          Each #P by 3 coordinates of the bar chains
  //  isoE:          #E by 2 edges of the bar
  //  isoOrigE:       #E by 3 edges of each bar in (f,if1,if2) format  (if - index of halfedge opposite to vertex if).
  //  isoN:          #E by 3 normals to the bars (corresponding to face normals at these points)
  //  width           bars widths
  //  heights:        bar heights (individual in #isoE)
  //  barColors:      #E by 3 RGB colors per bar
  // Output:
  //  V:              #V by 3 bar mesh coordinates
  //  T:              #T by 3 mesh triangles
  //  C:              #T by 3 colors
  
  bool IGL_INLINE bar_chains(const Eigen::MatrixXd& origV,
                             const Eigen::MatrixXi& origF,
                             const Eigen::MatrixXd& isoV,
                             const Eigen::MatrixXi& isoE,
                             const Eigen::MatrixXi& isoOrigE,
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
    
    VectorXd origA;
    igl::doublearea(origV,origF, origA);
    

    TBar<<0,1,2,
    2,3,0;
    
    V.resize(VBar.rows()*isoE.rows(),3);
    T.resize(TBar.rows()*isoE.rows(),3);
    int NewColorSize=T.rows();
    C.resize(NewColorSize,3);
    
    for (int i=0;i<isoE.rows();i++){
      
      RowVector3d translation = heights(i)*isoN.row(i);
      
      //std::cout<<"isoOrigE.row(i): "<<isoOrigE.row(i)<<std::endl;
      int v0e0=origF(isoOrigE(i,0),(isoOrigE(i,1)+1)%3);
      int v1e0=origF(isoOrigE(i,0),(isoOrigE(i,1)+2)%3);
      int v0e1=origF(isoOrigE(i,0),(isoOrigE(i,2)+1)%3);
      int v1e1=origF(isoOrigE(i,0),(isoOrigE(i,2)+2)%3);
      
      RowVector3d e0Vec = origV.row(v1e0)-origV.row(v0e0);
      RowVector3d e1Vec = origV.row(v1e1)-origV.row(v0e1);
      double e0Length = e0Vec.norm();
      double e1Length = e1Vec.norm();
      RowVector3d e0NormVec = e0Vec/e0Length;
      RowVector3d e1NormVec = e1Vec/e1Length;
      
      RowVector3d p12 = isoV.row(isoE(i,1))-isoV.row(isoE(i,0));
      double p12Length = p12.norm();
      
      double baryCoords0 = (isoV.row(isoE(i,0))-origV.row(v0e0)).norm()/e0Length;
      double baryCoords1 = (isoV.row(isoE(i,1))-origV.row(v0e1)).norm()/e1Length;
      
      
      double cosWith0 = p12.dot(e0NormVec)/p12Length;
      double cosWith1 = p12.dot(e1NormVec)/p12Length;
      
      
      double lengthDeviation0=width/(sqrt(1.0-cosWith0*cosWith0))/2.0;
      double lengthDeviation1=width/(sqrt(1.0-cosWith1*cosWith1))/2.0;
      
      double baryDeviation0 =lengthDeviation0/e0Length;
      double baryDeviation1 =lengthDeviation1/e1Length;
      
      double baryTop0 = (baryCoords0+baryDeviation0 > 1.0 ? 1.0 : baryCoords0+baryDeviation0);
      double baryBottom0 = (baryCoords0-baryDeviation0 < 0.0 ? 0.0 : baryCoords0-baryDeviation0);
      double baryTop1 = (baryCoords1+baryDeviation1 > 1.0 ? 1.0 : baryCoords1+baryDeviation1);
      double baryBottom1 = (baryCoords1-baryDeviation1 < 0.0 ? 0.0 : baryCoords1-baryDeviation1);
      
      
      VBar<<origV.row(v0e0).array() + (origV.row(v1e0)-origV.row(v0e0)).array()*baryBottom0,
      origV.row(v0e0).array() + (origV.row(v1e0)-origV.row(v0e0)).array()*baryTop0,
      origV.row(v0e1).array() + (origV.row(v1e1)-origV.row(v0e1)).array()*baryBottom1,
      origV.row(v0e1).array() + (origV.row(v1e1)-origV.row(v0e1)).array()*baryTop1;
      
      V.block(VBar.rows()*i,0,VBar.rows(),3)=VBar+translation.replicate(VBar.rows(),1);
      T.block(TBar.rows()*i,0,TBar.rows(),3)=TBar.array()+VBar.rows()*i;
      C.block(TBar.rows()*i,0,TBar.rows(),3)=barColors.row(i).replicate(TBar.rows(),1);
    }
    
    return true;
  }
  

  }




#endif


