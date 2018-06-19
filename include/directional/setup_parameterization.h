// This file is part of libdirectional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_SETUP_PARAMETERIZATION_H
#define DIRECTIONAL_SETUP_PARAMETERIZATION_H
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
  //Eigen::MatrixXd& cutV,
  // Eigen::MatrixXi& cutF,
  //Eigen::VectorXi& cut2wholeIndices,
  //Eigen::VectorXi& edge2TransitionIndices,
  IGL_INLINE void setup_parameterization(const int N,
                                         const Eigen::MatrixXd& wholeV,
                                         const Eigen::MatrixXi& wholeF,
                                         const Eigen::VectorXi& matching,
                                         const Eigen::VectorXi& singPositions,
                                         const Eigen::MatrixXi& face2cut,
                                         Eigen::SparseMatrix<double>& ind2vertexTransMat,
                                         Eigen::SparseMatrix<double>& vertexTrans2CutMat,
                                         Eigen::SparseMatrix<double>& constraintMat,
                                         Eigen::VectorXi& constrainedVertices,
                                         Eigen::MatrixXd& cutV,
                                         Eigen::MatrixXi& cutF)
  {
    
    using namespace Eigen;
    using namespace std;
    
    MatrixXi EV, FE, EF, EFi,EH, FH;
    MatrixXd FEs;
    VectorXi VH, HV, HE, HF, nextH,prevH,twinH,innerEdges;
    
    VectorXi D=VectorXi::Constant(wholeF.rows(),3);
    VectorXi isSingular=VectorXi::Zero(wholeV.rows());
    for (int i=0;i<singPositions.size();i++)
      isSingular(singPositions(i))=1;
    
    constrainedVertices=VectorXi::Zero(wholeV.rows());
    hedra::polygonal_edge_topology(D,wholeF,EV,FE,EF,EFi, FEs,innerEdges);
    hedra::dcel(D,wholeF,EV,EF,EFi,innerEdges,VH,EH,FH,HV,HE,HF,nextH,prevH, twinH);
    
    vector<MatrixXi> constParmMatrices(N);
    MatrixXi unitPermMatrix=MatrixXi::Zero(N,N);
    for (int i=0;i<N;i++)
      unitPermMatrix((i+1)%N,i)=1;
    
    constParmMatrices[0]=MatrixXi::Identity(N,N);
    for (int i=1;i<N;i++)
      constParmMatrices[i]=unitPermMatrix*constParmMatrices[i-1];
    
    VectorXi isHEcut=VectorXi::Zero(HE.rows());
    VectorXi cutValence=VectorXi::Zero(wholeV.rows());
    for (int i=0;i<wholeF.rows();i++)
      for (int j=0;j<3;j++)
        if (face2cut(i,j)){
          cutValence(wholeF(i,j))++;
          isHEcut(FH(i,j))=1;
        }
    
    
    //establishing transition variables by tracing cut curves
    VectorXi Halfedge2TransitionIndices=VectorXi::Constant(HE.rows(),32767);
    VectorXi Halfedge2Matching(HE.rows());
    VectorXi isHEClaimed=VectorXi::Zero(HE.rows());
    //matching: EF(i,0)->EF(i,1), halfedge on i get the incoming matching, and so opposite.
    
    for (int i=0;i<HE.rows();i++){
      Halfedge2Matching(i)=(EH(HE(i),0)==i ? -matching(HE(i)) : matching(HE(i)));
      while (Halfedge2Matching(i)<0) Halfedge2Matching(i)+=N;
      //while (Halfedge2Matching(i)>=N/2) Halfedge2Matching(i)-=N;
    }
    
    int currTransition=1;
    
    //cutting mesh and creating map between wholeF and cutF
    //cutting the mesh
    vector<int> cut2whole;
    vector<RowVector3d> cutVlist;
    cutF.conservativeResize(wholeF.rows(),3);
    for (int i=0;i<VH.rows();i++){
      //creating corners whereever we have non-trivial matching
      int beginH=VH(i);
      int currH=beginH;
      
      //reseting to first cut, if exists
      do{
        if (isHEcut(currH)!=0)
          break;
        currH=twinH(prevH(currH));
      }while (beginH!=currH);
      
      beginH=currH;

      do{
        if ((isHEcut(currH)!=0)||(beginH==currH)){
          cut2whole.push_back(i);
          cutVlist.push_back(wholeV.row(i));
        }
        
        for (int j=0;j<3;j++)
          if (wholeF(HF(currH),j)==i)
            cutF(HF(currH),j)=cut2whole.size()-1;
        
        currH=twinH(prevH(currH));
      }while (beginH!=currH);
    }
    
    cout<<"cutF: "<<cutF<<endl;
    /*cut2wholeIndices.conservativeResize(cut2whole.size());
    for (int i=0;i<cut2wholeIndices.size();i++)
      cut2wholeIndices(i)=cut2whole[i];*/
    
    cutV.conservativeResize(cutVlist.size(),3);
    for (int i=0;i<cutVlist.size();i++)
      cutV.row(i)=cutVlist[i];
    
    //starting from each cut-graph node, we trace cut curves
    cout<<"wholeV.rows(): "<<wholeV.rows()<<endl;
    cout<<"wholeF.rows(): "<<wholeF.rows()<<endl;
    for (int i=0;i<wholeV.rows();i++){
      if (((cutValence(i)==2)&&(!isSingular(i)))||(cutValence(i)==0))
        continue;  //either mid-cut curve or non at all
      
      //cout<<"Starting to trace from vertex "<<i<<" cut valence "<<cutValence(i)<<endl;
      //cout<<"currTransition: "<<currTransition<<endl;
      //tracing curves until next node, if not already filled
      int beginH=VH(i);
      int currH=beginH;
      int nextHalfedgeInCut=-1;
      do{
        if ((isHEcut(currH)!=0)&&(isHEClaimed(currH)==0)) { //unclaimed halfedge
          nextHalfedgeInCut=currH;
          Halfedge2TransitionIndices(nextHalfedgeInCut)=currTransition;
          Halfedge2TransitionIndices(twinH(nextHalfedgeInCut))=-currTransition;
          isHEClaimed(nextHalfedgeInCut)=1;
          isHEClaimed(twinH(nextHalfedgeInCut))=1;
          int nextCutVertex=HV(nextH(nextHalfedgeInCut));
          //cout<<"nextCutVertex: "<<nextCutVertex<<endl;
          //cout<<"cutValence(nextCutVertex): "<<cutValence(nextCutVertex)<<endl;
          //cout<<"isSingular(nextCutVertex)"<<isSingular(nextCutVertex)<<endl;
          //advancing on the cut until next node
          while ((cutValence(nextCutVertex)==2)&&(!isSingular(nextCutVertex))){
            int beginH=VH(nextCutVertex);
            int currH=beginH;
            int nextHalfedgeInCut=-1;
            do{
              if ((isHEcut(currH)!=0)&&(isHEClaimed(currH)==0)) { //unclaimed cut halfedge
                nextHalfedgeInCut=currH;
                break;
              }
              currH=twinH(prevH(currH));
            }while (beginH!=currH);
            //cout<<"nextHalfedgeInCut: "<<nextHalfedgeInCut<<endl;
            Halfedge2TransitionIndices(nextHalfedgeInCut)=currTransition;
            Halfedge2TransitionIndices(twinH(nextHalfedgeInCut))=-currTransition;
            isHEClaimed(nextHalfedgeInCut)=1;
            isHEClaimed(twinH(nextHalfedgeInCut))=1;
            nextCutVertex=HV(nextH(nextHalfedgeInCut));
            //cout<<"nextCutVertex: "<<nextCutVertex<<endl;
            //cout<<"cutValence(nextCutVertex): "<<cutValence(nextCutVertex)<<endl;
            //cout<<"isSingular(nextCutVertex)"<<isSingular(nextCutVertex)<<endl;
          }
          
          currTransition++;
        }
        currH=twinH(prevH(currH));
      }while (beginH!=currH);
    }
    
    int numTransitions=currTransition-1;
    
    //checking
    /*for (int i=0;i<HE.rows();i++)
      if ((matching(HE(i))!=0)&&(Halfedge2TransitionIndices(i)==32767)){
        cout<<"HV(i), cutValence(HV(i)): "<<HV(i)<<","<<cutValence(HV(i))<<endl;
        cout<<"HV(nextH(i)), cutValence(HV(nextH(i))): "<<HV(nextH(i))<<","<<cutValence(HV(nextH(i)))<<endl;
      }*/
    
    //establishing independent variables
    //all singularity positions are dependent of transition variables, and the number of independent variables is therefore regular vertices + transitions - 1
    VectorXi independentVars(wholeV.rows()+numTransitions-singPositions.rows());
    int currIndVar=0;
    for (int i=0;i<wholeV.rows();i++){
      if (!isSingular(i))
          independentVars(currIndVar++)=i;
    }
    
    for (int i=0;i<numTransitions;i++)
      independentVars(currIndVar++)=wholeV.rows()+i;
    
    vector<Triplet<double> > vertexTrans2CutTriplets, constTriplets,  ind2vertexTransTriplets;
    
    //the identity parts of ind2vertexTranMat
    for (int i=0;i<wholeV.rows()-singPositions.rows()+numTransitions;i++)
      for (int k=0;k<N;k++)
        ind2vertexTransTriplets.push_back(Triplet<double>(N*independentVars(i)+k, N*i+k, 1.0));
    
    //forming the constraints and the singularity positions
    int currConst=0;
    for (int i=0;i<VH.rows();i++){
      std::vector<MatrixXi> permMatrices;
      std::vector<int> permIndices;  //in the space #V + #transitions
      //The initial corner gets the identity without any transition
      permMatrices.push_back(MatrixXi::Identity(N,N));
      permIndices.push_back(i);
      int beginH=VH(i);
      int currH=beginH;
       //if (cutValence(i)==2)
        // cout<<"tracking a 2-valence vertex"<<endl;
      int currCutVertex=-1;
      do{
        int currFace=HF(currH);
        int newCutVertex=-1;
        for (int j=0;j<3;j++)
          if (wholeF(currFace,j)==i)
            newCutVertex=cutF(currFace,j);
        
        //currCorner gets the permutations so far
        if (newCutVertex!=currCutVertex){
          currCutVertex=newCutVertex;
          
          for (int i=0;i<permIndices.size();i++)
            for (int j=0;j<N;j++)
              for (int k=0;k<N;k++)
                vertexTrans2CutTriplets.push_back(Triplet<double>(N*currCutVertex+j, N*permIndices[i]+k, (double)permMatrices[i](j,k)));
        }
        //updating the matrices for the next corner
        int nextHalfedge=twinH(prevH(currH));
        MatrixXi nextPermMatrix = constParmMatrices[Halfedge2Matching(nextHalfedge)%N];
        //cout<<"Halfedge2Matching(nextHalfedge): "<<Halfedge2Matching(nextHalfedge)<<endl;
        if (isHEcut(nextHalfedge)==0) { //no update needed
          currH=nextHalfedge;
          continue;
        }
        
        //otherwise, updating matrices with transition
        int nextTransition = Halfedge2TransitionIndices(nextHalfedge);
        // cout<<"nextTransition: "<<nextTransition<<endl;
        if (nextTransition>0){  //Pe*f + Je
          for (int j=0;j<permMatrices.size();j++)
            permMatrices[j]=nextPermMatrix*permMatrices[j];
          
          //and identity on the fresh transition
          permMatrices.push_back(MatrixXi::Identity(N,N));
          permIndices.push_back(wholeV.rows()+nextTransition-1);
        } else { // (Pe*(f-Je))  matrix is already inverse since halfedge matching is minused
          //reverse order
          permMatrices.push_back(-MatrixXi::Identity(N,N));
          permIndices.push_back(wholeV.rows()-nextTransition-1);
          
          for (int j=0;j<permMatrices.size();j++)
            permMatrices[j]=nextPermMatrix*permMatrices[j];
          
        }
        currH=nextHalfedge;
      }while(currH!=beginH);
      
      //cleaning parmMatrices and permIndices to see if there is a constraint or reveal singularity-from-ransition
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
        if (isSingular(i)){  //the position of the singularity depends on a combination of transition variables
          MatrixXi singMatrix=cleanPermMatrices[0];  //assuming that's the first vertex
          for (int j=1;j<cleanPermMatrices.size();j++){
            MatrixXi currMatrix=-singMatrix.transpose()*cleanPermMatrices[j];
            for (int k=0;k<N;k++)
              for (int l=0;l<N;l++)
                ind2vertexTransTriplets.push_back(Triplet<double>(N*i+k, N*(cleanPermIndices[j]-singPositions.rows())+l, (double)currMatrix(k,l)));
          }
        }else{  //non-singular node ties between transition variables
          for (int j=0;j<cleanPermMatrices.size();j++)
            for (int k=0;k<N;k++)
              for (int l=0;l<N;l++)
                constTriplets.push_back(Triplet<double>(N*currConst+k, N*cleanPermIndices[j]+l, (double)cleanPermMatrices[j](k,l)));
          currConst++;
        }
        cout<<"found constraint with: "<<endl;
        cout<<"cutValence(i): "<<cutValence(i)<<endl;
        cout<<"isSingular(i): "<<isSingular(i)<<endl;
        for (int j=0;j<permIndices.size();j++){
          cout<<"permIndices[j]: "<<permIndices[j]<<endl;
          cout<<"permMatrices[j]: "<<permMatrices[j]<<endl;
        }
        for (int j=0;j<cleanPermIndices.size();j++){
          cout<<"cleanPermIndices[j]: "<<cleanPermIndices[j]<<endl;
          cout<<"cleanPermMatrices[j]: "<<cleanPermMatrices[j]<<endl;
        }
        constrainedVertices(i)=1;
      }
    }
    
    cout<<"currConst: "<<currConst<<endl;
    
    ind2vertexTransMat.conservativeResize(N*(wholeV.rows()+numTransitions), N*(wholeV.rows()-singPositions.rows()+numTransitions));
    vector<Triplet<double>> cleanTriplets;
    for (int i=0;i<ind2vertexTransTriplets.size();i++)
      if (ind2vertexTransTriplets[i].value()!=0.0)
        cleanTriplets.push_back(ind2vertexTransTriplets[i]);
    ind2vertexTransMat.setFromTriplets(cleanTriplets.begin(), cleanTriplets.end());
    
    vertexTrans2CutMat.conservativeResize(N*cutV.rows(), N*(wholeV.rows()+numTransitions));
    cleanTriplets.clear();
    for (int i=0;i<vertexTrans2CutTriplets.size();i++)
      if (vertexTrans2CutTriplets[i].value()!=0.0)
        cleanTriplets.push_back(vertexTrans2CutTriplets[i]);
    vertexTrans2CutMat.setFromTriplets(cleanTriplets.begin(), cleanTriplets.end());
    
    constraintMat.conservativeResize(N*currConst, N*(wholeV.rows()+numTransitions));
    cleanTriplets.clear();
    for (int i=0;i<constTriplets.size();i++)
      if (constTriplets[i].value()!=0.0)
        cleanTriplets.push_back(constTriplets[i]);
    constraintMat.setFromTriplets(cleanTriplets.begin(), cleanTriplets.end());
    
  }
}




#endif


