// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.



#ifndef DIRECTIONAL_BRANCHED_ISOLINES_H
#define DIRECTIONAL_BRANCHED_ISOLINES_H
#include <igl/igl_inline.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <vector>
#include <igl/remove_duplicate_vertices.h>
#include <directional/isolines.h>
#include <directional/visualization_schemes.h>


namespace directional{

//Traces isolines for a branched function defined on the vertices of a (cut) mesh.
//Input:
//V:          |V| x 3 vertex coordinates for the mesh
//F:          |F| x 3 face vertex indices (into V)
//NFunction:  |V| x N branched function values
//Output:
//VIsolines:  vertices of isoline visualization mesh
//FIsolines:  faces -"-
//CIsolines:  colors -"-

  void branched_isolines(const Eigen::MatrixXd& V,
                         const Eigen::MatrixXi& F,
                         const Eigen::MatrixXd& NFunction,
                         Eigen::MatrixXd& VIsolines,
                         Eigen::MatrixXi& FIsolines,
                         Eigen::MatrixXd& CIsolines)
  {
    
    int N = NFunction.cols();
    Eigen::MatrixXd funcColors=directional::default_glyph_colors(8);
    double isolineRadius=0.02;
    int jumps = (N%2 == 0 ? 2 : 1);
    Eigen::MatrixXd isoV;
    Eigen::MatrixXi isoE;
    VIsolines.resize(0,3); FIsolines.resize(0,3); CIsolines.resize(0,3);
    double l = 1.25*igl::avg_edge_length(V, F);
    
    Eigen::MatrixXd P1(0,3), P2(0,3), C(0,3);
    for (int i=0;i<NFunction.cols()/jumps;i++){
      Eigen::VectorXd currFunc = NFunction.col(i);
      igl::isolines(V,F, currFunc, 100, isoV, isoE);
      
      int oldPSize = P1.rows();
      P1.conservativeResize(P1.rows()+isoE.rows(),3);
      P2.conservativeResize(P2.rows()+isoE.rows(),3);
      C.conservativeResize(C.rows()+isoE.rows(),3);
      for (int j=0;j<isoE.rows();j++){
        P1.row(oldPSize+j)=isoV.row(isoE(j,0));
        P2.row(oldPSize+j)=isoV.row(isoE(j,1));
        C.row(oldPSize+j)=funcColors.row(i);
      }

    }
    
    directional::line_cylinders(P1, P2, l*isolineRadius,C,4, VIsolines, FIsolines, CIsolines);
          
  }


}

#endif
