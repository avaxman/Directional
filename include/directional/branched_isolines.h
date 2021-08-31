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
//isoV:       coordinates of isolines
//isoE:       #E by 2 connectivity edges of isolines
//isoN:       #E by 3 normals to isolines (normals to mesh at these lines).
//funcNum: identity of function (as NFunction #col) of the corresponding P1 (or P2) entry.

  void branched_isolines(const Eigen::MatrixXd& V,
                         const Eigen::MatrixXi& F,
                         const Eigen::MatrixXd& NFunction,
                         Eigen::MatrixXd& isoV,
                         Eigen::MatrixXi& isoE,
                         Eigen::MatrixXd& isoN,
                         Eigen::VectorXi& funcNum)
  {
    
    int N = NFunction.cols();
    int jumps = (N%2 == 0 ? 2 : 1);
    Eigen::MatrixXd isoVPart, isoNPart;
    Eigen::MatrixXi isoEPart;
    
    for (int i=0;i<NFunction.cols()/jumps;i++){
      Eigen::VectorXd currFunc = NFunction.col(i);
      igl::isolines(V,F, currFunc, 100, isoVPart, isoEPart, isoNPart);
      
      int oldVSize = isoV.rows();
      int oldESize = isoE.rows();
      isoV.conservativeResize(oldVSize+isoVPart.rows(),3);
      isoE.conservativeResize(oldESize+isoEPart.rows(),2);
      isoN.conservativeResize(oldESize+isoNPart.rows(),3);
      funcNum.conservativeResize(oldESize+isoEPart.rows());
      
      isoV.block(oldVSize,0,isoVPart.rows(),3) = isoVPart;
      isoE.block(oldESize,0,isoEPart.rows(),2) = isoEPart.array()+oldVSize;
      isoN.block(oldESize,0,isoNPart.rows(),3) = isoNPart;
      funcNum.tail(isoNPart.rows()).setConstant(i);
    }
  }
}

#endif
