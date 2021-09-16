// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_VERTEX_HIGHLIGHTS_H
#define DIRECTIONAL_VERTEX_HIGHLIGHTS_H

#include <Eigen/Core>
#include <string>
#include <vector>
#include <cmath>
#include <complex>
#include <igl/igl_inline.h>
#include <igl/vertex_triangle_adjacency.h>

namespace directional
  {
  // creates a mesh of vertex small corner triangles that highlight a vertex, and that can complement a halfedge highlights for a more complete look
  // Input:
  //  V, F:          original mesh
  //  hlVertices:    indices of the highlighted vertices into V
  //  hlColors:      #hlVertices colors
  //  widthRatio:     width of highlight bar w.r.t. width of highlight bar w.r.t. opposite edge height in face
  // Output:
  //  hlV:              #V by 3 highlight mesh coordinates
  //  hlT:              #T by 3 mesh triangles
  //  hlC:              #T by 3 colors
  
  bool IGL_INLINE vertex_highlights(const Eigen::MatrixXd& V,
                                    const Eigen::MatrixXi& F,
                                    const Eigen::MatrixXi& hlVertices,
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
    
    std::vector<std::vector<int>> VF, VFi;
    igl::vertex_triangle_adjacency(V, F,VF,VFi);
    
    //counting highlights
    int numHighlights=0;
    for (int i=0;i<hlVertices.size();i++)
      numHighlights+=VF[hlVertices(i)].size();
    
    hlV.resize(numHighlights*3,3);
    hlT.resize(numHighlights,3);
    hlC.resize(numHighlights,3);
    
    numHighlights = 0;
    for (int i=0;i<hlVertices.size();i++){
      for (int j=0;j<VF[hlVertices(i)].size();j++)
      {
       
        MatrixXd VCornerTri(3,3);
        
        int face=VF[hlVertices(i)][j], inFace=VFi[hlVertices(i)][j];
        
        //std::cout<<"F(face,inFace):"<<F(face,inFace)<<std::endl;
        //std::cout<<"hlVertices(i): "<<hlVertices(i)<<std::endl;

        //edge coordinates
        VCornerTri.row(0)=V.row(F(face,inFace));
        //coordinates into the other edges
        VCornerTri.row(1) = V.row(F(face,(inFace+1)%3))*widthRatio+V.row(F(face,inFace))*(1.0-widthRatio);
        VCornerTri.row(2) = V.row(F(face,(inFace+2)%3))*widthRatio+V.row(F(face,inFace))*(1.0-widthRatio);
        
        VCornerTri=VCornerTri.array()+normals.row(i).replicate(VCornerTri.rows(),1).array()*height;

        hlV.block(numHighlights*VCornerTri.rows(),0,VCornerTri.rows(),3)=VCornerTri;
        hlT.row(numHighlights)<<3*numHighlights, 3*numHighlights+1, 3*numHighlights+2;
        hlC.row(numHighlights)=hlColors.row(i);
        
        numHighlights++;
      }
      
    }
    
    return true;
  }
  

  }




#endif


