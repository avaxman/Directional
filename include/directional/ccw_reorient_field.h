// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_CCW_REORIENT_FIELD_H
#define DIRECTIONAL_CCW_REORIENT_FIELD_H

#include <iostream>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/local_basis.h>
#include <unsupported/Eigen/Polynomials>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>


namespace directional
{
  // Converts a field in PolyVector representation to raw represenation.
  // Inputs:
  //  B1, B2:           #F by 3 matrices representing the local base of each face.
  //  polyVectorField:  #F by N complex PolyVectors
  //  N:                The degree of the field.
  // Outputs:
  //  raw:              #F by 3*N matrix with all N explicit vectors of each directional in raw format xyzxyz
  IGL_INLINE void ccw_reorient_field(const Eigen::MatrixXd& B1,
                                    const Eigen::MatrixXd& B2,
                                    const Eigen::MatrixXd& rawField,
                                    Eigen::MatrixXd& reOrientRawField)
  {
    reOrientRawField.resize(rawField.rows(), rawField.cols());
    int N=rawField.cols()/3;
    for (int f = 0; f < B1.rows(); f++)
    {
      Eigen::VectorXcd complexField(N);
      for (int i=0;i<N;i++){
        Eigen::RowVector3d currVector=rawField.block(f,3*i,1,3);
        complexField(i)=std::complex<double>(currVector.dot(B1.row(f)), currVector.dot(B2.row(f)));
      }
      std::sort(complexField.data(), complexField.data() + complexField.size(), [](std::complex<double> a, std::complex<double> b){return arg(a) < arg(b);});
      
      for (int i=0;i<N;i++)
        reOrientRawField.block(f,3*i,1,3)=B1.row(f)*complexField(i).real()+B2.row(f)*complexField(i).imag();
    }
  }
}
  
#endif
