// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef PRINCIPAL_COMBING_H
#define PRINCIPAL_COMBING_H
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_topology.h>
#include <directional/tree.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>

#include <Eigen/Core>
#include <queue>
#include <vector>
#include <cmath>


namespace directional
{
  // Takes a field in raw form and computes both the principal effort and the consequent principal matching on every edge.
  // Note: equals to rotation_angle*N in N-RoSy fields
  // Important: if the Raw field in not CCW ordered (e.., resulting from all functions X_t_raw() in libdirectional), the result is unpredictable.
  // Input:
  //  V:      #V x 3 vertex coordinates
  //  F:      #F x 3 face vertex indices
  //  EV:     #E x 2 edges to vertices indices
  //  EF:     #E x 2 edges to faces indices
  //  raw:    The directional field, assumed to be ordered CCW, and in xyzxyzxyz...xyz (3*N cols) form. The degree is inferred by the size.
  // Output:
  // matching: #E matching function, where vector k in EF(i,0) matches to vector (i+matching(i))%N in EF(i,1). In case of boundary, there is a -1.
  //  effort: #E principal matching efforts.
  //TODO: also return matching
  IGL_INLINE void principal_combing(const Eigen::MatrixXd& V,
                             const Eigen::MatrixXi& F,
                             const Eigen::MatrixXi& EV,
                             const Eigen::MatrixXi& EF,
                             const Eigen::MatrixXi& FE,
                             const Eigen::MatrixXd& rawField,
                             Eigen::MatrixXd& combedField,
                             Eigen::VectorXi& combedMatching,
                             Eigen::VectorXd& combedEffort)
  {
    using namespace Eigen;
    VectorXi matching;
    VectorXd effort;
    principal_matching(V,F,EV, EF,FE,rawField,matching,effort);
    //flood-filling through the matching to comb field
    combedField.conservativeResize(rawField.rows(), rawField.cols());
    int N=rawField.cols()/3;
    //dual tree to find combing routes
    VectorXi visitedFaces=VectorXi::Constant(F.rows(),1,0);
    std::queue<std::pair<int,int> > faceMatchingQueue;
    faceMatchingQueue.push(std::pair<int,int>(0,0));
    do{
      std::pair<int,int> currFaceMatching=faceMatchingQueue.front();
      faceMatchingQueue.pop();
      if (visitedFaces(currFaceMatching.first))
        continue;
      visitedFaces(currFaceMatching.first)=1;
      
      //combing field to start from the matching index
      combedField.block(currFaceMatching.first, 0, 1, 3*(N-currFaceMatching.second))=rawField.block(currFaceMatching.first, 3*currFaceMatching.second, 1, 3*(N-currFaceMatching.second));
      combedField.block(currFaceMatching.first, 3*(N-currFaceMatching.second), 1, 3*currFaceMatching.second)=rawField.block(currFaceMatching.first, 0, 1, 3*currFaceMatching.second);
      
      for (int i=0;i<3;i++){
        int nextMatching=(matching(FE(currFaceMatching.first,i)));
        int nextFace=(EF(FE(currFaceMatching.first,i),0)==currFaceMatching.first ? EF(FE(currFaceMatching.first,i),1) : EF(FE(currFaceMatching.first,i),0));
        nextMatching*=(EF(FE(currFaceMatching.first,i),0)==currFaceMatching.first ? 1.0 : -1.0);
        nextMatching=(nextMatching+currFaceMatching.second+10*N)%N;  //killing negatives
        if ((nextFace!=-1)&&(!visitedFaces(nextFace)))
          faceMatchingQueue.push(std::pair<int,int>(nextFace, nextMatching));
        
      }
      
    }while (!faceMatchingQueue.empty());
    
    //can be produced from the combing, but better kept for sanity check.
    principal_matching(V, F, EV, EF, FE, combedField, combedMatching, combedEffort);
  }
  
  //representative version
  IGL_INLINE void principal_combing(const Eigen::MatrixXd& V,
                             const Eigen::MatrixXi& F,
                             const Eigen::MatrixXi& EV,
                             const Eigen::MatrixXi& EF,
                             const Eigen::MatrixXi& FE,
                             const Eigen::MatrixXd& representativeField,
                             const int N,
                             Eigen::MatrixXd& combedField,
                             Eigen::VectorXi& combedMatching,
                             Eigen::VectorXd& combedEffort)
  {
    Eigen::MatrixXd raw;
    representative_to_raw(V, F, representativeField, N, raw);
    principal_combing(V, F, EV, EF, FE, raw, combedField, combedMatching, combedEffort);
  }
  
  
  
  // (V, F) only raw version
  IGL_INLINE void principal_combing(const Eigen::MatrixXd& V,
                             const Eigen::MatrixXi& F,
                             const Eigen::MatrixXd& raw,
                             Eigen::MatrixXd& combedField,
                             Eigen::VectorXi& combedMatching,
                             Eigen::VectorXd& combedEffort)
  {
    Eigen::MatrixXi EV, FE, EF;
    igl::edge_topology(V, F, EV, FE, EF);
    principal_combing(V, F, EV, EF, FE, raw, combedField, combedMatching, combedEffort);
  }
  
  // (V,F) only representative version
  IGL_INLINE void principal_combing(const Eigen::MatrixXd& V,
                             const Eigen::MatrixXi& F,
                             const Eigen::MatrixXd& representativeField,
                             const int N,
                             Eigen::MatrixXd& combedField,
                             Eigen::VectorXi& combedMatching,
                             Eigen::VectorXd& combedEffort)
  {
    Eigen::MatrixXi EV, FE, EF;
    igl::edge_topology(V, F, EV, FE, EF);
    principal_combing(V, F, EV, EF, FE, representativeField, N, combedField, combedMatching, combedEffort);
  }
}




#endif


