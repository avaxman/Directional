// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_COMBING_H
#define DIRECTIONAL_COMBING_H


#include <Eigen/Core>
#include <queue>
#include <vector>
#include <cmath>
#include <igl/igl_inline.h>
#include <directional/CartesianField.h>
#include <directional/tree.h>
#include <directional/principal_matching.h>

namespace directional
{
  // Reorders the vectors in a face (preserving CCW) so that the prescribed matching across most edges, except a small set (called a cut), is an identity, making it ready for cutting and parameterization.
  // Important: if the Raw field in not CCW ordered, the result is unpredictable.
  // Input:
  //  V:        #V x 3 vertex coordinates
  //  F:        #F x 3 face vertex indices
  //  EV:       #E x 2 edges to vertices indices
  //  EF:       #E x 2 edges to faces indices
  //  rawField: #F by 3*N  The directional field, assumed to be ordered CCW, and in xyzxyz raw format. The degree is inferred by the size.
  //  matching: #E matching function, where vector k in EF(i,0) matches to vector (k+matching(k))%N in EF(i,1). In case of boundary, there is a -1.
  // Output:
  //  combedField: #F by 3*N reindexed field
  IGL_INLINE void combing(const directional::CartesianField& rawField,
                          directional::CartesianField& combedField,
                          const Eigen::MatrixXi& _faceIsCut=Eigen::MatrixXi())
  {
    using namespace Eigen;
    combedField.init_field(*(rawField.mesh), RAW_FIELD, rawField.N);
    Eigen::MatrixXi faceIsCut(rawField.intField.rows(),3);
    if (_faceisCut.rows()==0)
      faceIsCut.setZero();
    else
      faceIsCut=_faceIsCut;
    //flood-filling through the matching to comb field
    //combedField.extField.conservativeResize(rawField.rows(), rawField.cols());
    //int N=rawField.cols()/3;
    //dual tree to find combing routes
    VectorXi visitedSpaces=VectorXi::Constant(rawField.intField.rows(),1,0);
    std::queue<std::pair<int,int> > spaceMatchingQueue;
    spaceMatchingQueue.push(std::pair<int,int>(0,0));
    MatrixXd combedIntField(combedField.intField.rows(), combedField.intField.cols());
    do{
      std::pair<int,int> currSpaceMatching=spaceMatchingQueue.front();
      spaceMatchingQueue.pop();
      if (visitedSpaces(currSpaceMatching.first))
        continue;
      visitedSpaces(currSpaceMatching.first)=1;
      
      //combing field to start from the matching index
      combedIntField.block(currSpaceMatching.first, 0, 1, 2*(rawField.N-currSpaceMatching.second))=rawField.intField.block(currSpaceMatching.first, 2*currSpaceMatching.second, 1, 2*(rawField.N-currSpaceMatching.second));
      combedIntField.block(currSpaceMatching.first, 2*(rawField.N-currSpaceMatching.second), 1, 2*currSpaceMatching.second)=rawField.intField.block(currSpaceMatching.first, 0, 1, 2*currSpaceMatching.second);
      
      for (int i=0;i<3;i++){
        int nextMatching=(rawField.matching(rawField.oneRing(currSpaceMatching.first,i)));
        int nextFace=(rawField.adjSpaces(rawField.oneRing(currSpaceMatching.first,i),0)==currSpaceMatching.first ? rawField.adjSpaces(rawField.oneRing(currSpaceMatching.first,i),1) : rawField.adjSpaces(rawField.oneRing(currSpaceMatching.first,i),0));
        nextMatching*=(rawField.adjSpaces(rawField.oneRing(currSpaceMatching.first,i),0)==currSpaceMatching.first ? 1.0 : -1.0);
        nextMatching=(nextMatching+currSpaceMatching.second+10*rawField.N)%rawField.N;  //killing negatives
        if ((nextFace!=-1)&&(!visitedSpaces(nextFace))&&(!faceIsCut(currSpaceMatching.first,i)))
          spaceMatchingQueue.push(std::pair<int,int>(nextFace, nextMatching));
        
      }
      
    }while (!spaceMatchingQueue.empty());
    
    combedField.set_intrinsic_field(combedIntField);
  }


  
  //version for input in representative format (for N-RoSy directionals).
  /*IGL_INLINE void combing(const Eigen::MatrixXd& V,
                          const Eigen::MatrixXi& F,
                          const Eigen::MatrixXi& EV,
                          const Eigen::MatrixXi& EF,
                          const Eigen::MatrixXi& FE,
                          const Eigen::MatrixXd& representativeField,
                          const int N,
                          const Eigen::VectorXi& matching,
                          Eigen::MatrixXd& combedField)
  {
    Eigen::MatrixXd rawField;
    representative_to_raw(V, F, representativeField, N, rawField);
    combing(V, F, EV, EF, FE, rawField, matching, combedField);
  }*/
  
  //version with prescribed cuts from faces
  /*IGL_INLINE void combing(const directional::CartesianField& rawField,
                          const Eigen::MatrixXi& faceIsCut,
                          directional::CartesianField& combedField)
  {
    using namespace Eigen;
    //flood-filling through the matching to comb field
    combedField.conservativeResize(rawField.rows(), rawField.cols());
    combedMatching.conservativeResize(EF.rows());
    int N=rawField.cols()/3;
    //dual tree to find combing routes
    VectorXi visitedSpaces=VectorXi::Constant(F.rows(),1,0);
    std::queue<std::pair<int,int> > spaceMatchingQueue;
    spaceMatchingQueue.push(std::pair<int,int>(0,0));
    VectorXi faceTurns(rawField.rows());
    do{
      std::pair<int,int> currSpaceMatching=spaceMatchingQueue.front();
      spaceMatchingQueue.pop();
      if (visitedSpaces(currSpaceMatching.first))
        continue;
      visitedSpaces(currSpaceMatching.first)=1;
      
      //combing field to start from the matching index
      combedField.block(currSpaceMatching.first, 0, 1, 3*(N-currSpaceMatching.second))=rawField.block(currSpaceMatching.first, 3*currSpaceMatching.second, 1, 3*(N-currSpaceMatching.second));
      combedField.block(currSpaceMatching.first, 3*(N-currSpaceMatching.second), 1, 3*currSpaceMatching.second)=rawField.block(currSpaceMatching.first, 0, 1, 3*currSpaceMatching.second);
      
      faceTurns(currSpaceMatching.first)=currSpaceMatching.second;
      
      for (int i=0;i<3;i++){
        int nextMatching=(matching(FE(currSpaceMatching.first,i)));
        int nextFace=(EF(FE(currSpaceMatching.first,i),0)==currSpaceMatching.first ? EF(FE(currSpaceMatching.first,i),1) : EF(FE(currSpaceMatching.first,i),0));
        nextMatching*=(EF(FE(currSpaceMatching.first,i),0)==currSpaceMatching.first ? 1.0 : -1.0);
        nextMatching=(nextMatching+currSpaceMatching.second+10*N)%N;  //killing negatives
        if ((nextFace!=-1)&&(!visitedSpaces(nextFace))&&(!faceIsCut(currSpaceMatching.first,i)))
          spaceMatchingQueue.push(std::pair<int,int>(nextFace, nextMatching));
        
      }
      
    }while (!spaceMatchingQueue.empty());
    
    //giving combed matching
    for (int i=0;i<EF.rows();i++){
      if ((EF(i,0)==-1)||(EF(i,1)==-1))
        combedMatching(i)=-1;
      else
        combedMatching(i)=(faceTurns(EF(i,0))-faceTurns(EF(i,1))+matching(i)+1000000*N)%N;
    }
  }*/
  
}




#endif


