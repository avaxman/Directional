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
  // Takes a field in raw form and computes both the curl-matching effort and the consequent curl matching on every edge.
  // Important: if the Raw field in not CCW ordered, the result is meaningless.
  // Input:
  //  V:      #V x 3 vertex coordinates
  //  F:      #F x 3 face vertex indices
  //  EV:     #E x 2 edges to vertices indices
  //  EF:     #E x 2 edges to faces indices
  //  raw:    The directional field, assumed to be ordered CCW, and in xyzxyzxyz...xyz (3*N cols) form. The degree is inferred by the size.
  // Output:
  // matching: #E matching function, where vector k in EF(i,0) matches to vector (k+matching(k))%N in EF(i,1). In case of boundary, there is a -1.
  //  effort: #E principal matching efforts.
  // curlNorm: the L2-norm of the curl vector
  //  singVertices: indices (into V) of which vertices are singular; including boundary vertices which carry the singularity of their loop
  //  singIndices: the index of the singular vertices (corresponding with singIndices), relative to N (the true index is then i/N).
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
      /*double currEffort=0;
      for (int j = 0; j < rawField.N; j++) {
        RowVector3d vecjf = rawField.block(EF(i, 0), 3*j, 1, 3);
        Complex vecjfc = Complex(vecjf.dot(B1.row(EF(i, 0))), vecjf.dot(B2.row(EF(i, 0))));
        RowVector3d vecjg = rawField.block(EF(i, 1), 3 * ((matching(i)+j+N)%N), 1, 3);
        Complex vecjgc = Complex(vecjg.dot(B1.row(EF(i, 1))), vecjg.dot(B2.row(EF(i, 1))));
        Complex transvecjfc = vecjfc*edgeTransport(i);
        //cout<<"transvecjfc, vecjgc: "<<transvecjfc<<","<<vecjgc<<endl;
        currEffort+= arg(vecjgc / transvecjfc);
        //cout<<"arg(vecjgc / transvecjfc): "<<arg(vecjgc / transvecjfc)<<endl;
      }
      
      effort(i) = currEffort;
      */
    }
    
    //Getting final singularities and their indices
    effort_to_indices(rawField);
  }
  
  //version with representative vector (for N-RoSy) as input.
  /*IGL_INLINE void curl_matching(const Eigen::MatrixXd& V,
                                const Eigen::MatrixXi& F,
                                const Eigen::MatrixXi& EV,
                                const Eigen::MatrixXi& EF,
                                const Eigen::MatrixXi& FE,
                                const Eigen::MatrixXd& representativeField,
                                const int N,
                                Eigen::VectorXi& matching,
                                Eigen::VectorXd& effort,
                                Eigen::VectorXd& curlNorm,
                                Eigen::VectorXi& singVertices,
                                Eigen::VectorXi& singIndices)
  {
    Eigen::MatrixXd rawField;
    representative_to_raw(V, F, representativeField, N, rawField);
    curl_matching(V, F, EV, EF, FE, rawField, matching, effort, curlNorm,singVertices, singIndices);
  }*/
  }




#endif


