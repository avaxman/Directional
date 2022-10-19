// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef HEDRA_POLYGONAL_WRITE_OFF_H
#define HEDRA_POLYGONAL_WRITE_OFF_H
#include <igl/igl_inline.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace hedra
{
  // writes a polygonal mesh as an ascii OFF file
  // Inputs:
  //   str  path to .off file
  //  V  eigen double matrix  #V by 3 - vertex coordinates
  //  D  eigen int vector     #F by 1 - face degrees
  //  F  eigen int matrix     #F by max(D) - vertex indices in face
  IGL_INLINE bool polygonal_write_OFF(const std::string& str,
                                      const Eigen::MatrixXd& V,
                                      const Eigen::VectorXi& D,
                                      const Eigen::MatrixXi& F)
  {
    
    using namespace std;
    using namespace Eigen;
    ofstream FileHandle;
    FileHandle.open(str);
    if (!FileHandle.is_open())
      return false;
    FileHandle<<"OFF"<<endl<<V.rows()<<" "<<F.rows()<<" 0"<<endl;
    FileHandle<<V<<endl;
    MatrixXi FD(D.rows(), D.cols()+F.cols());
    FD<<D, F;
    for (int i=0;i<F.rows();i++)
      FileHandle<<FD.block(i,0,1,D(i)+1)<<endl;
    FileHandle.close();
    return true;
  }
}


#endif


