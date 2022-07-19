// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_PRINCIPAL_MATCHING_H
#define DIRECTIONAL_PRINCIPAL_MATCHING_H

#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <directional/representative_to_raw.h>
#include <directional/effort_to_indices.h>
#include <directional/TangentBundle.h>
#include <directional/definitions.h>

namespace directional
{
  // Takes a field in raw form and computes both the principal effort and the consequent principal matching on every edge.
  // Important: if the Raw field in not CCW ordered, the result is meaningless.
  // The input and output are both a RAW_FIELD type cartesian field, in which the matching, effort, and singularities are set.
  IGL_INLINE void principal_matching(directional::CartesianField& field)
  {
    
    typedef std::complex<double> Complex;
    using namespace Eigen;
    using namespace std;
    
    field.matching.conservativeResize(field.tb->adjSpaces.rows());
    field.matching.setConstant(-1);
    
    field.effort = VectorXd::Zero(field.tb->adjSpaces.rows());
    for (int i = 0; i < field.tb->adjSpaces.rows(); i++) {
      if (field.tb->adjSpaces(i, 0) == -1 || field.tb->adjSpaces(i, 1) == -1)
        continue;
    
      double minRotAngle=10000.0;
      int indexMinFromZero=0;
      
      //computing some effort and the extracting principal one
      Complex freeCoeff(1.0,0.0);
      //finding where the 0 vector in EF(i,0) goes to with smallest rotation angle in EF(i,1), computing the effort, and then adjusting the matching to have principal effort.
      RowVector2d vec0f = field.intField.block(field.tb->adjSpaces(i, 0), 0, 1, 2);
      Complex vec0fc = Complex(vec0f(0), vec0f(1));
      Complex transvec0fc = vec0fc*field.tb->connection(i);
      for (int j = 0; j < field.N; j++) {
        RowVector2d vecjf = field.intField.block(field.tb->adjSpaces(i, 0), 2 * j, 1, 2);
        Complex vecjfc = Complex(vecjf(0),vecjf(1));
        RowVector2d vecjg = field.intField.block(field.tb->adjSpaces(i, 1), 2 * j, 1, 2);
        Complex vecjgc = Complex(vecjg(0),vecjg(1));
        Complex transvecjfc = vecjfc*field.tb->connection(i);
        freeCoeff *= (vecjgc / transvecjfc);
        double currRotAngle =arg(vecjgc / transvec0fc);
        if (abs(currRotAngle)<abs(minRotAngle)){
          indexMinFromZero=j;
          minRotAngle=currRotAngle;
        }
        
        //taking principal effort
        
      }
      field.effort(i) = arg(freeCoeff);
      
      //finding the matching that implements effort(i)
      //This is still not perfect
      double currEffort=0;
      for (int j = 0; j < field.N; j++) {
        RowVector2d vecjf = field.intField.block(field.tb->adjSpaces(i, 0), 2*j, 1, 2);
        Complex vecjfc = Complex(vecjf(0), vecjf(1));
        RowVector2d vecjg = field.intField.block(field.tb->adjSpaces(i, 1), 2 *((j+indexMinFromZero+field.N)%field.N), 1, 2);
        Complex vecjgc = Complex(vecjg(0), vecjg(1));
        Complex transvecjfc = vecjfc*field.tb->connection(i);
        currEffort+= arg(vecjgc / transvecjfc);
      }
   
      field.matching(i)=indexMinFromZero-round((currEffort-field.effort(i))/(2.0*igl::PI));
    }
          
    //Getting final singularities and their indices
    effort_to_indices(field);
    
  }
}




#endif


