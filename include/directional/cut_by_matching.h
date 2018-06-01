// This file is part of libdirectional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_CUT_BY_MATCHING_H
#define DIRECTIONAL_CUT_BY_MATCHING_H
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_topology.h>
#include <directional/tree.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/dcel.h>
#include <directional/polygonal_edge_topology.h>

#include <Eigen/Core>
#include <queue>
#include <vector>
#include <cmath>


namespace directional
{
  //This doesn't work with boundaries!
  // Reorders the vectors in a face (preserving CCW) so that the principal matching across most edges, except a small set (called a cut), is an identity, making it ready for cutting and parameterization.
  // Important: if the Raw field in not CCW ordered, the result is unpredictable.
  // Input:
  //  wholeV:      #V x 3 vertex coordinates
  //  wholeF:      #F x 3 face vertex indices
  //  EV:     #E x 2 edges to vertices indices
  //  EF:     #E x 2 edges to faces indices
  // matching: #E matching function, where vector k in EF(i,0) matches to vector (k+matching(k))%N in EF(i,1). In case of boundary, there is a -1. Expect most matching =0 due to the combing.
  // Output:

  //
  IGL_INLINE void cut_by_matching(const int N,
                                  const Eigen::MatrixXd& wholeV,
                                  const Eigen::MatrixXi& wholeF,
                                  const Eigen::VectorXi& matching,
                                  const Eigen::VectorXi& singIndices,
                                  Eigen::MatrixXd& cutV,
                                  Eigen::MatrixXi& cutF,
                                  Eigen::VectorXi& cut2wholeIndices,
                                  Eigen::VectorXi& edge2TransitionIndices,
                                  Eigen::SparseMatrix<double>& vt2cMat,
                                  Eigen::SparseMatrix<double>& constraintMat)
  {
    
    using namespace Eigen;
    using namespace std;

    MatrixXi EV, FE, EF, EFi,EH, FH;
    MatrixXd FEs;
    VectorXi VH, HV, HE, HF, nextH,prevH,twinH,innerEdges;
    
    VectorXi D=VectorXi::Constant(wholeV.rows(),3);
    
    hedra::polygonal_edge_topology(D,wholeF,EV,FE,EF,EFi, FEs,innerEdges);
    hedra::dcel(D,wholeF,EV,EF,EFi,innerEdges,VH,EH,FH,HV,HE,HF,nextH,prevH, twinH);
    
    vector<MatrixXi> constParmMatrices(N);
    MatrixXi unitPermMatrix=MatrixXi::Zero(N,N);
    for (int i=0;i<N;i++)
      unitPermMatrix((i+1)%N,i)=1.0;
    
    constParmMatrices[0]=MatrixXi::Identity(N,N);
    for (int i=1;i<N;i++)
      constParmMatrices[i]=unitPermMatrix*constParmMatrices[i-1];
    
    //establishing transition variables by tracing cut cut curves
    VectorXi Halfedge2TransitionIndices=VectorXi::Constant(HE.rows(),-1);
    VectorXi Halfedge2Matching=VectorXi::Constant(HE.rows(),-1);
    int currTransition=0;
    
    VectorXi cutValence=VectorXi::Zero(EV.rows());
    for (int i=0;i<EV.rows();i++){
      if (matching(i)!=0){
        cutValence(EV(i,0))++;
        cutValence(EV(i,1))++;
      }
    }
    
    //matching: EF(i,0)->EF(i,1), halfedge on i get the incoming matching, and so opposite.
    for (int i=0;i<HE.rows();i++)
      Halfedge2Matching(i)=(EH(HE(i),0)==i ? -matching(HE(i)) : matching(HE(i)));
    
    //starting from each node, we trace curves
    for (int i=0;i<wholeV.rows();i++){
      if ((cutValence(i)==2)||(cutValence(i)==0))
        continue;  //either mid-cut curve or non at all
      
      //tracing curves until next node, if not already filled
      int beginH=VH(i);
      int currH=beginH;
      int nextHalfedgeInCut=-1;
      do{
        if ((matching(HE(currH))!=0)&&(Halfedge2TransitionIndices(currH)==-1)) { //unclaimed halfedge
          nextHalfedgeInCut=currH;
          break;
        }
        currH=twinH(prevH(currH));
      }while (beginH!=currH);
      
      if (nextHalfedgeInCut==-1)
        continue;  //all are assigned already from other directions
      
      Halfedge2TransitionIndices(nextHalfedgeInCut)=currTransition;
      Halfedge2TransitionIndices(twinH(nextHalfedgeInCut))=-currTransition;
      int nextCutVertex=HV(nextH(nextHalfedgeInCut));
      //advancing on the cut until next node
      while (cutValence(nextCutVertex)==2){
        int beginH=VH(i);
        int currH=beginH;
        int nextHalfedgeInCut=-1;
        do{
          if ((matching(HE(currH))!=0)&&(Halfedge2TransitionIndices(currH)==-1)) { //unclaimed cut halfedge
            nextHalfedgeInCut=currH;
            break;
          }
          currH=twinH(prevH(currH));
        }while (beginH!=currH);
        Halfedge2TransitionIndices(nextHalfedgeInCut)=currTransition;
        Halfedge2TransitionIndices(twinH(nextHalfedgeInCut))=-currTransition;
        nextCutVertex=HV(nextH(nextHalfedgeInCut));
      }
      
      currTransition++;
    }
    
    //m2v is minimal amount of variables after constraint factorization
    vector<Triplet<double> > v2cTriplets, m2vTriplets, constTriplets;
    int currConst=0;
    for (int i=0;i<VH.rows();i++){
      std::vector<MatrixXi> permMatrices;
      std::vector<int> permIndices;  //in the space #V + #transitions
      //The initial corner gets the identity without any transition
      permMatrices.push_back(MatrixXi::Identity(N,N));
      permIndices.push_back(i);
      int beginH=VH(i);
      int currH=beginH;
      do{
        int currFace=HF(currH);
        int currCorner=-1;
        for (int j=0;j<3;j++)
          if (wholeF(currFace,j)==i)
            currCorner=3*currFace+j;
        
        //currCorner gets the permutations so far
        for (int i=0;i<permIndices.size();i++)
          for (int j=0;j<N;j++)
            for (int k=0;k<N;k++)
              v2cTriplets.push_back(Triplet<double>(N*permIndices[i]+j, N*permIndices[i]+k, (double)permMatrices[i](j,k)));
        
        //updating the matrices for the next corner
        int nextHalfedge=twinH(prevH(currH));
        MatrixXi nextPermMatrix = constParmMatrices[(Halfedge2Matching(nextHalfedge)+N)%N];
        int nextTransition = Halfedge2TransitionIndices(nextHalfedge);
        if (nextTransition >=0){  //Pe*f + Je
          for (int j=0;j<permMatrices.size();j++)
            permMatrices[j]=nextPermMatrix*permMatrices[j];
          
          //and identity on the fresh transition
          permMatrices.push_back(MatrixXi::Identity(N,N));
          permIndices.push_back(nextTransition);
        } else { // (Pe*(f-Je)) (Pe is alreay inverse since matching is signed on halfedges
          //reverse order
          permMatrices.push_back(-MatrixXi::Identity(N,N));
          permIndices.push_back(-nextTransition);
          
          for (int j=0;j<permMatrices.size();j++)
            permMatrices[j]=nextPermMatrix*permMatrices[j];
          
        }

        currH=nextHalfedge;
      }while(currH!=beginH);
      
      //cleaning parmMatrices and permIndices to see if there is a constraint
      
      std::set<int> cleanPermIndicesSet(permIndices.begin(), permIndices.end());
      std::vector<int> cleanPermIndices(cleanPermIndicesSet.begin(), cleanPermIndicesSet.end());
      std::vector<MatrixXi> cleanPermMatrices(cleanPermIndices.size());
      
      for (int j=0;j<cleanPermIndices.size();j++){
        cleanPermMatrices[j]=MatrixXi::Zero(N,N);
        for (int k=0;k<permIndices.size();k++)
          if (cleanPermIndices[j]==permIndices[k])
            cleanPermMatrices[j]+=permMatrices[k];
        if (cleanPermIndices[j]==i)
          cleanPermMatrices[j]-=MatrixXi::Identity(N,N);
      }
      
      //if not all matrices are zero, there is a constraint
      bool isConstraint=false;
      for (int j=0;j<cleanPermMatrices.size();j++)
        if (cleanPermMatrices[j].cwiseAbs().maxCoeff()!=0)
          isConstraint=true;
      
      if (isConstraint){
        for (int j=0;j<cleanPermMatrices.size();j++)
          for (int k=0;k<N;k++)
            for (int l=0;l<N;l++)
              constTriplets.push_back(Triplet<double>(N*currConst+k, N*cleanPermIndices[j]+l, (double)cleanPermMatrices[j](k,l)));
        currConst++;
            
      }
    }
    
    vt2cMat.conservativeResize(3*N*wholeF.rows(), wholeV.rows()+currTransition);
    vt2cMat.setFromTriplets(v2cTriplets.begin(), v2cTriplets.end());
    
    constraintMat.conservativeResize(N*currConst, wholeV.rows()+currTransition);
    constraintMat.setFromTriplets(constTriplets.begin(), constTriplets.end());
    
    //cutting the mesh
    vector<int> cut2whole;
    vector<RowVector3d> cutVlist;
    cutF.conservativeResize(wholeF.rows(),3);
    for (int i=0;i<VH.rows();i++){
      //creating corners whereever we have non-trivial matching
      int beginH=VH(i);
      int currH=beginH;
      
      do{
        if (matching(HE(currH))!=0){
          cut2whole.push_back(i);
          cutVlist.push_back(wholeV.row(i));
        }
        
        for (int j=0;j<3;j++)
          if (wholeF(HF(currH),j)==i)
            cutF(HF(currH),j)=cut2whole.size()-1;
        
        currH=twinH(prevH(currH));
      }while (beginH!=currH);
    }
    
    cut2wholeIndices.conservativeResize(cut2whole.size());
    for (int i=0;i<cut2wholeIndices.size();i++)
      cut2wholeIndices(i)=cut2whole[i];
    
    cutV.conservativeResize(cutVlist.size(),3);
    for (int i=0;i<cutVlist.size();i++)
      cutV.row(i)=cutVlist[i];
    
  }
}




#endif


