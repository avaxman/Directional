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

namespace directional{

//Traces isolines for a branched function defined on the vertices of a (cut) mesh.
//Input:
//V:          |V| x 3 vertex coordinates for the mesh
//F:          |F| x 3 face vertex indices (into V)
//NFunction:  |V| x N branched function values
//Output:
//P1,P2:  start and end points of isoline
//funcNum: identity of function (as NFunction #col) of the corresponding P1 (or P2) entry.

  void branched_isolines(const Eigen::MatrixXd& V,
                         const Eigen::MatrixXi& F,
                         const Eigen::MatrixXd& NFunction,
                         Eigen::MatrixXd& P1,
                         Eigen::MatrixXd& P2,
                         Eigen::VectorXi& funcNum)
  {
    
    int N = NFunction.cols();
    int jumps = (N%2 == 0 ? 2 : 1);
    Eigen::MatrixXd isoV;
    Eigen::MatrixXi isoE;
    
    P1.resize(0,3), P2.resize(0,3), funcNum.resize(0);
    for (int i=0;i<NFunction.cols()/jumps;i++){
      Eigen::VectorXd currFunc = NFunction.col(i);
      igl::isolines(V,F, currFunc, 100, isoV, isoE);
      
      int oldPSize = P1.rows();
      P1.conservativeResize(P1.rows()+isoE.rows(),3);
      P2.conservativeResize(P2.rows()+isoE.rows(),3);
      funcNum.conservativeResize(funcNum.rows()+isoE.rows());
      for (int j=0;j<isoE.rows();j++){
        P1.row(oldPSize+j)=isoV.row(isoE(j,0));
        P2.row(oldPSize+j)=isoV.row(isoE(j,1));
        funcNum(oldPSize+j)=i;
      }
    }
  }
}

#endif
