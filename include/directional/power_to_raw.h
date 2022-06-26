// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POWER_TO_RAW_H
#define DIRECTIONAL_POWER_TO_RAW_H

#include <igl/igl_inline.h>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>

namespace directional
{
  // Converts the power complex representation to raw representation.
  // Input:
  //  B1, B2, B3: bases for each face from igl::local_base().
  //  powerField: #F x 1 Representation of the field as complex numbers
  //  N: the degree of the field.
  // normalize: whether to produce a normalized result (length = 1)
  // Output:
  //  rawField: #F by 3*N matrix with all N explicit vectors of each directional in the order X,Y,Z,X,Y,Z, ...
  IGL_INLINE void power_to_raw(const directional::CartesianField& powerField,
                               int N,
                               directional::CartesianField& rawField,
                               bool normalize=false)
  {
    assert(powerField.fieldType==POWER_FIELD && "The input field should be a power/PolyVector field");
    rawField.init_field(*(powerField.mesh), RAW_FIELD,N);
    //std::cout<<"powerField.intField: "<<powerField.intField<<std::endl;
    Eigen::MatrixXcd intFieldComplex(powerField.intField.rows(),N);
    Eigen::VectorXcd complexPowerField(powerField.intField.rows());
    //power fields are represented as -u^N since they are a special case of PVs.
    complexPowerField.array().real()=-powerField.intField.col(0);
    complexPowerField.array().imag()=-powerField.intField.col(1);
    intFieldComplex.col(0)=pow(complexPowerField.array(),1.0/(double)N);
    for (int i=1;i<N;i++)
      intFieldComplex.col(i)=intFieldComplex.col(0)*exp(std::complex<double>(0,2*igl::PI*(double)i/(double)N));
    
    //std::cout<<"intFieldComplex: "<<intFieldComplex<<std::endl;
   
    if (normalize)
      intFieldComplex.array()/=intFieldComplex.array().abs();
    
    //std::cout<<"intFieldComplex: "<<intFieldComplex<<std::endl;
    
    Eigen::MatrixXd intField(intFieldComplex.rows(),2*N);
    for (int i=0;i<N;i++){
      intField.col(2*i)=intFieldComplex.col(i).real();
      intField.col(2*i+1)=intFieldComplex.col(i).imag();
    }
    //std::cout<<"intField: "<<intField<<std::endl;
    
    rawField.set_intrinsic_field(intField);
  }
  
}

#endif
