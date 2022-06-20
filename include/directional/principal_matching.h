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
  IGL_INLINE void principal_matching(directional::FaceField& field)
  {
    
    typedef std::complex<double> Complex;
    using namespace Eigen;
    using namespace std;
    
    field.matching.conservativeResize(field.mesh->EF.rows());
    field.matching.setConstant(-1);
    
    VectorXcd edgeTransport(field.mesh->EF.rows());  //the difference in the angle representation of edge i from EF(i,0) to EF(i,1)
    MatrixXd edgeVectors(field.mesh->EF.rows(), 3);
    for (int i = 0; i < field.mesh->EF.rows(); i++) {
      if (field.mesh->EF(i, 0) == -1 || field.mesh->EF(i, 1) == -1)
        continue;
      edgeVectors.row(i) = (field.mesh->V.row(field.mesh->EV(i, 1)) - field.mesh->V.row(field.mesh->EV(i, 0))).normalized();
      Complex ef(edgeVectors.row(i).dot(field.mesh->Bx.row(field.mesh->EF(i, 0))), edgeVectors.row(i).dot(field.mesh->By.row(field.mesh->EF(i, 0))));
      Complex eg(edgeVectors.row(i).dot(field.mesh->Bx.row(field.mesh->EF(i, 1))), edgeVectors.row(i).dot(field.mesh->By.row(field.mesh->EF(i, 1))));
      edgeTransport(i) = eg / ef;
    }
    
    field.effort = VectorXd::Zero(field.mesh->EF.rows());
    for (int i = 0; i < field.mesh->EF.rows(); i++) {
      if (field.mesh->EF(i, 0) == -1 || field.mesh->EF(i, 1) == -1)
        continue;
      //computing free coefficient effort (a.k.a. [Diamanti et al. 2014])
      //Complex freeCoeffEffort(1.0, 0.0);
      double minRotAngle=10000.0;
      int indexMinFromZero=0;
      
      //computing some effort and the extracting principal one
      Complex freeCoeff(1.0,0.0);
      //finding where the 0 vector in EF(i,0) goes to with smallest rotation angle in EF(i,1), computing the effort, and then adjusting the matching to have principal effort.
      
      RowVector3d vec0f = field.extField.block(field.mesh->EF(i, 0), 0, 1, 3);
      Complex vec0fc = Complex(vec0f.dot(field.mesh->Bx.row(field.mesh->EF(i, 0))), vec0f.dot(field.mesh->By.row(field.mesh->EF(i, 0))));
      Complex transvec0fc = vec0fc*edgeTransport(i);
      for (int j = 0; j < field.N; j++) {
        RowVector3d vecjf = field.extField.block(field.mesh->EF(i, 0), 3 * j, 1, 3);
        Complex vecjfc = Complex(vecjf.dot(field.mesh->Bx.row(field.mesh->EF(i, 0))), vecjf.dot(field.mesh->By.row(field.mesh->EF(i, 0))));
        RowVector3d vecjg = field.extField.block(field.mesh->EF(i, 1), 3 * j, 1, 3);
        Complex vecjgc = Complex(vecjg.dot(field.mesh->Bx.row(field.mesh->EF(i, 1))), vecjg.dot(field.mesh->By.row(field.mesh->EF(i, 1))));
        Complex transvecjfc = vecjfc*edgeTransport(i);
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
        RowVector3d vecjf = field.extField.block(field.mesh->EF(i, 0), 3*j, 1, 3);
        Complex vecjfc = Complex(vecjf.dot(field.mesh->Bx.row(field.mesh->EF(i, 0))), vecjf.dot(field.mesh->By.row(field.mesh->EF(i, 0))));
        RowVector3d vecjg = field.extField.block(field.mesh->EF(i, 1), 3 *((j+indexMinFromZero+field.N)%field.N), 1, 3);
        Complex vecjgc = Complex(vecjg.dot(field.mesh->Bx.row(field.mesh->EF(i, 1))), vecjg.dot(field.mesh->By.row(field.mesh->EF(i, 1))));
        Complex transvecjfc = vecjfc*edgeTransport(i);
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


