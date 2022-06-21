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
#include <directional/TriMesh.h>
#include <directional/FaceField.h>

namespace directional
{
  // Takes a field in raw form and computes both the principal effort and the consequent principal matching on every edge.
  // Important: if the Raw field in not CCW ordered, the result is meaningless.
  // Input:
  //  V:      #V x 3 vertex coordinates
  //  F:      #F x 3 face vertex indices
  //  EV:     #E x 2 edges to vertices indices
  //  EF:     #E x 2 edges to faces indices
  //  raw:    The directional field, assumed to be ordered CCW, and in xyzxyzxyz...xyz (3*N cols) form. The degree is inferred by the size.
  // Output:
  //  matching: #E matching function, where vector k in EF(i,0) matches to vector (k+matching(k))%N in EF(i,1). In case of boundary, there is a -1.
  //  effort: #E principal matching efforts.
  //  singVertices: indices (into V) of which vertices are singular; including boundary vertices which carry the singularity of their loop
  //  singIndices: the index of the singular vertices (corresponding with singIndices), relative to N (the true index is then i/N). This discludes boundary vertices (boundary cycles have their own index along generator cycles)
  IGL_INLINE void principal_matching(directional::CartesianField& field)
  {
    
    typedef std::complex<double> Complex;
    using namespace Eigen;
    using namespace std;
    
    field.matching.conservativeResize(field.adjSpaces.rows());
    field.matching.setConstant(-1);
    
    field.effort = VectorXd::Zero(field.adjSpaces.rows());
    for (int i = 0; i < field.adjSpaces.rows(); i++) {
      if (field.adjSpaces(i, 0) == -1 || field.adjSpaces(i, 1) == -1)
        continue;
    
      double minRotAngle=10000.0;
      int indexMinFromZero=0;
      
      //computing some effort and the extracting principal one
      Complex freeCoeff(1.0,0.0);
      //finding where the 0 vector in EF(i,0) goes to with smallest rotation angle in EF(i,1), computing the effort, and then adjusting the matching to have principal effort.
      RowVector2d vec0f = field.intField.block(field.adjSpaces(i, 0), 0, 1, 2);
      Complex vec0fc = Complex(vec0f(0), vec0f(1));
      Complex transvec0fc = vec0fc*field.connection(i);
      for (int j = 0; j < field.N; j++) {
        RowVector2d vecjf = field.intField.block(field.adjSpaces(i, 0), 2 * j, 1, 2);
        Complex vecjfc = Complex(vecjf(0),vecjf(1));
        RowVector2d vecjg = field.intField.block(field.adjSpaces(i, 1), 2 * j, 1, 2);
        Complex vecjgc = Complex(vecjg(0),vecjg(1));
        Complex transvecjfc = vecjfc*field.connection(i);
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
        RowVector2d vecjf = field.intField.block(field.adjSpaces(i, 0), 2*j, 1, 2);
        Complex vecjfc = Complex(vecjf(0), vecjf(1));
        RowVector2d vecjg = field.intField.block(field.adjSpaces(i, 1), 2 *((j+indexMinFromZero+field.N)%field.N), 1, 2);
        Complex vecjgc = Complex(vecjg(0), vecjg(1));
        Complex transvecjfc = vecjfc*field.connection(i);
        currEffort+= arg(vecjgc / transvecjfc);
      }
   
      field.matching(i)=indexMinFromZero-round((currEffort-field.effort(i))/(2.0*igl::PI));
    }
          
    //Getting final singularities and their indices
    effort_to_indices(field);
    
  }
  
  //Version with representative vector (for N-RoSy alone) as input.
  /*IGL_INLINE void principal_matching(directional::FaceField& field,
                                     const Eigen::MatrixXd& representativeField,
                                     const int N)
  {
    Eigen::MatrixXd rawField;
    representative_to_raw(field.mesh->V, field.mesh->F, representativeField, N, rawField);
    field.set_field(rawField,field.mesh);
    principal_matching(field);
  }*/
}




#endif


