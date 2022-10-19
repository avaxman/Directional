// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_HALFEDGE_HIGHLIGHTS_H
#define DIRECTIONAL_HALFEDGE_HIGHLIGHTS_H

#include <Eigen/Core>
#include <string>
#include <vector>
#include <cmath>
#include <complex>
#include <igl/igl_inline.h>

namespace directional
  {
  // creates a mesh of conforming halfedge trapezoidal bars that highlight a halfedge with some color (for instance for showing seams)
  // Input:
  //  V, F:          original mesh
  //  hlHalfedges:    #F by 3 color indicators of the highlighted halfedges (corresponding to face corners where the edge i os corner [i,(i+1)%3]). -1 means not highlighted
  //  hlColors:      #scolors for highlighted, into which hlHalfedges indexes
  //  widthRatio:     width of highlight bar w.r.t. opposite edge height in face
  // Output:
  //  hlV:              #V by 3 highlight mesh coordinates
  //  hlT:              #T by 3 mesh triangles
  //  hlC:              #T by 3 colors
  
  bool IGL_INLINE halfedge_highlights(const Eigen::MatrixXd& V,
                                      const Eigen::MatrixXi& F,
                                      const Eigen::MatrixXi& hlHalfedges,
                                      const Eigen::MatrixXd& hlColors,
                                      Eigen::MatrixXd& hlV,
                                      Eigen::MatrixXi& hlT,
                                      Eigen::MatrixXd& hlC,
                                      const double widthRatio,
                                      const double height)
  {
    using namespace Eigen;

    MatrixXd normals;
     igl::per_face_normals(V, F, normals);
    
    //counting highlights
    int numHighlights=0;
    for (int i=0;i<hlHalfedges.rows();i++)
      for (int j=0;j<3;j++)
        if (hlHalfedges(i,j)!=-1)
          numHighlights++;
    
    hlV.resize(numHighlights*4,3);
    hlT.resize(numHighlights*2,3);
    hlC.resize(numHighlights*2,3);
    
    numHighlights = 0;
    for (int i=0;i<hlHalfedges.rows();i++){
      for (int j=0;j<3;j++)
      {
        if (hlHalfedges(i,j)==-1)
          continue;
        
        MatrixXd VTrapeze(4,3);
        MatrixXi TTrapeze(2,3);
        
        //edge coordinates
        VTrapeze.row(0)=V.row(F(i,j));
        VTrapeze.row(1)=V.row(F(i,(j+1)%3));
        //coordinates into the other edges
        VTrapeze.row(2) = V.row(F(i,(j+1)%3))*(1.0-widthRatio)+V.row(F(i,(j+2)%3))*(widthRatio);
        VTrapeze.row(3) = V.row(F(i,(j+2)%3))*(widthRatio)+V.row(F(i,j))*(1.0-widthRatio);
        
        VTrapeze=VTrapeze.array()+normals.row(i).replicate(VTrapeze.rows(),1).array()*height;
        
        TTrapeze<<0,1,2,
        2,3,0;
        
        hlV.block(numHighlights*VTrapeze.rows(),0,VTrapeze.rows(),3)=VTrapeze;
        hlT.block(numHighlights*TTrapeze.rows(),0,TTrapeze.rows(),3)=TTrapeze.array()+numHighlights*4;
        hlC.block(numHighlights*TTrapeze.rows(),0,TTrapeze.rows(),3)=hlColors.row(hlHalfedges(i,j)).replicate(TTrapeze.rows(),1);
        
        numHighlights++;
      }
      
    }
    
    return true;
  }
  

  }




#endif


