// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can

// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_CURL_MATCHING_H
#define DIRECTIONAL_CURL_MATCHING_H

#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <igl/igl_inline.h>
#include <directional/TriMesh.h>
#include <directional/FaceField.h>
#include <directional/effort_to_indices.h>

namespace directional
  {
  // Takes a field in raw form and computes both the curl-matching effort and the consequent curl matching on every tangent-space adjacency.
  // Important: if the Raw field in not CCW ordered, the result is meaningless.
  // Note: curl is only (future work...) defined for face-based fields.
  // Input:
  //  rawField:  a RAW_FIELD type *face-based* field.
  // Output:
  // curlNorm: the L2-norm of the curl vector
  // rawField: the input field matching, effort, and singularities are altered.

  IGL_INLINE void curl_matching(directional::FaceField& rawField,
                                Eigen::VectorXd& curlNorm)
  {
    
    typedef std::complex<double> Complex;
    using namespace Eigen;
    using namespace std;
    
    rawField.matching.conservativeResize(rawField.mesh->EF.rows());
    rawField.matching.setConstant(-1);
    curlNorm.conservativeResize(rawField.mesh->EF.rows());
    
    MatrixXd edgeVectors(rawField.mesh->EF.rows(), 3);
    for (int i = 0; i < rawField.mesh->EF.rows(); i++) {
      if (rawField.mesh->EF(i, 0) == -1 || rawField.mesh->EF(i, 1) == -1)
        continue;
      edgeVectors.row(i) = (rawField.mesh->V.row(rawField.mesh->EV(i, 1)) - rawField.mesh->V.row(rawField.mesh->EV(i, 0))).normalized();

    }
    
    //effort = VectorXd::Zero(rawField.mesh->EF.rows());
    for (int i = 0; i < rawField.mesh->EF.rows(); i++) {
      if (rawField.mesh->EF(i, 0) == -1 || rawField.mesh->EF(i, 1) == -1)
        continue;
      //computing free coefficient effort (a.k.a. [Diamanti et al. 2014])
      //Complex freeCoeffEffort(1.0, 0.0);
      int indexMinFromZero=0;
      //finding where the 0 vector in EF(i,0) goes to with smallest rotation angle in EF(i,1), computing the effort, and then adjusting the matching to have principal effort.
      double minCurl = 32767000.0;
      for (int j = 0; j < rawField.N; j++) {
        double currCurl = 0;
        for (int k=0;k<rawField.N;k++){
          RowVector3d vecDiff =rawField.extField.block(rawField.mesh->EF(i, 1), 3 * ((j+k)%rawField.N), 1, 3)-rawField.extField.block(rawField.mesh->EF(i, 0), 3*k, 1, 3);
          currCurl +=pow(edgeVectors.row(i).dot(vecDiff),2.0);
        }
        
        if (currCurl < minCurl){
          indexMinFromZero=j;
          minCurl=currCurl;
        }
      }
      
      rawField.matching(i) =indexMinFromZero;
      curlNorm(i)= sqrt(minCurl);
      
      //computing the full effort for 0->indexMinFromZero, and readjusting the matching to fit principal effort
      Complex freeCoeff(1,0);
      rawField.effort.resize(rawField.matching.size());
      for (int j = 0; j < rawField.N; j++) {
        //RowVector3d vecjf = rawField.extField.block(rawField.adjSpaces(i, 0), 3*j, 1, 3);
        RowVector2d vecjf = rawField.intField.block(rawField.adjSpaces(i, 0), 2 * j, 1, 2);
        Complex vecjfc = Complex(vecjf(0),vecjf(1));
        RowVector2d vecjg = rawField.intField.block(rawField.adjSpaces(i, 1), 2 * (rawField.matching(i)+j+rawField.N)%rawField.N, 1, 2);
        Complex vecjgc = Complex(vecjg(0),vecjg(1));
        Complex transvecjfc = vecjfc*rawField.connection(i);
        freeCoeff *= (vecjgc / transvecjfc);
      }
      
      rawField.effort(i) = arg(freeCoeff);
      
    }
    
    //Getting final singularities and their indices
    effort_to_indices(rawField);
    
  }
  
}




#endif


