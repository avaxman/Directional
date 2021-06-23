// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_SETUP_INTEGRATION_H
#define DIRECTIONAL_SETUP_INTEGRATION_H

#include <queue>
#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_topology.h>
#include <directional/tree.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/dcel.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/combing.h>

namespace directional
{
    
  struct IntegrationData
  {
    int N;  //# uncompressed parametric functions
    int n;  //# independent parameteric functions
    Eigen::MatrixXi linRed; //Linear Reduction tying the n dofs to the full N
    Eigen::MatrixXi periodMat;  //function spanning integers
    Eigen::SparseMatrix<double> vertexTrans2CutMat;   //a map between the whole mesh (vertex + translational jump) representation to the vertex-based representation on the cut mesh
    Eigen::SparseMatrix<double> constraintMat;  //linear constraints (resulting from non-singular nodes)
    Eigen::SparseMatrix<double> linRedMat;  //the global uncompression of n->N
    Eigen::SparseMatrix<double> intSpanMat;  //Spanning the translational jump lattice
    Eigen::SparseMatrix<double> singIntSpanMat;  //the layer for the singularities
    Eigen::VectorXi constrainedVertices;         //constrained vertices (fixed points in the parameterization)
    Eigen::VectorXi integerVars;                 //variables that are to be rounded.
    Eigen::MatrixXi face2cut;                    //|F|x3 map of which edges of faces are seams
    
    Eigen::VectorXi fixedIndices;  //the translation fixing indices
    Eigen::VectorXd fixedValues;   //translation fixed values
    Eigen::VectorXi singularIndices;   //the singular-vertex indices
    
    //integer versions, for pure seamless parameterizations
    Eigen::SparseMatrix<int> vertexTrans2CutMatInteger;
    Eigen::SparseMatrix<int> constraintMatInteger;
    Eigen::SparseMatrix<int> linRedMatInteger;
    Eigen::SparseMatrix<int> intSpanMatInteger;
    Eigen::SparseMatrix<int> singIntSpanMatInteger;
    
    double lengthRatio;     //global scaling of functions
    
    //Flags
    bool integralSeamless;  //If to the full translational seamless.
    bool roundSeams;        //If to round seams or round singularities
    bool verbose;           //output the integration optimization log.
    bool localInjectivity;  //Enforce local injectivity; might result in failure!
    
    IntegrationData(int _N):lengthRatio(0.02), integralSeamless(false), roundSeams(true), verbose(false), localInjectivity(false){
      N=_N;
      n=(N%2==0 ? N/2 : N);
      if (N%2==0)
        set_sign_symmetry(N);
      else linRed=Eigen::MatrixXi::Identity(N,n);
      set_default_period_matrix(n);
    }
    ~IntegrationData(){}
    
    IGL_INLINE void set_linear_reduction(const Eigen::MatrixXi& _linRed, const Eigen::MatrixXi& _periodMat){linRed =_linRed; N=linRed.rows(); n=linRed.cols(); periodMat=_periodMat;}
    
    //the default symmetry, where for even N there are N/2 lines
    IGL_INLINE void set_sign_symmetry(int N){
      assert(N%2==0);
      linRed.resize(N,N/2);
      linRed<<Eigen::MatrixXi::Identity(N/2,N/2),-Eigen::MatrixXi::Identity(N/2,N/2);
      n=N/2;
      set_default_period_matrix(n);
    }
    
    //the entire first N/3 lines are symmetric w.r.t. to the next two (N/3) packets, and where if N is even we also add sign symmetry.
    IGL_INLINE void set_triangular_symmetry(int N){
      assert(N%3==0);
      if (N%2==0){
        linRed.resize(N,N/3);
        linRed.block(0,0,N/2,N/3)<<Eigen::MatrixXi::Identity(N/3,N/3),-Eigen::MatrixXi::Identity(N/6,N/6),Eigen::MatrixXi::Identity(N/6,N/6);
        linRed.block(N/2,0,N/2,N/3)=-linRed.block(0,0,N/2,N/3);
        n=N/3;
      } else {
        linRed.resize(N,2*N/3);
        linRed<<Eigen::MatrixXi::Identity(2*N/3,2*N/3),-Eigen::MatrixXi::Identity(N/3,N/3),-Eigen::MatrixXi::Identity(N/3,N/3);
        n=2*N/3;
      }
      set_default_period_matrix(n);
    }
    
    IGL_INLINE void set_default_period_matrix(int n){
      periodMat=Eigen::MatrixXi::Identity(n,n);
    }
  };
  
  
  
  
  
  
  // Setting up the seamless integration algorithm
  // Input:
  //  wholeV:       #V x 3 vertex coordinates
  //  wholeF:       #F x 3 face vertex indices
  //  EV:           #E x 2 edges to vertices indices
  //  EF:           #E x 2 edges to faces indices
  //  FE:           #F x 3 faces to edges tindices
  // matching:      #E matching function, where vector k in EF(i,0) matches to vector (k+matching(k))%N in EF(i,1). In case of boundary, there is a -1. Most matching should be zero due to prior combing.
  // singVertices:  list of singular vertices in wholeV.
  // intData:       Integration data structure
  // Output:
  //  intData:      updated integration data.
  //  cutV:         the Vertices of the cut mesh.
  //  cutF:         the Faces of the cut mesh (1-1 correspondence with wholeF, but vertices indexed into cutV).
  //  combedField:  The raw field combed into N different fields on the cut mesh (every column is a single-vf).
  //  combedMatching: the new matching of the combed field when given on the whole mesh (mostly zero except on cuts).
  
  IGL_INLINE void setup_integration(const Eigen::MatrixXd& wholeV,
                                    const Eigen::MatrixXi& wholeF,
                                    const Eigen::MatrixXi& EV,
                                    const Eigen::MatrixXi& EF,
                                    const Eigen::MatrixXi& FE,
                                    const Eigen::MatrixXd& rawField,
                                    const Eigen::VectorXi& matching,
                                    const Eigen::VectorXi& singVertices,
                                    IntegrationData& intData,
                                    Eigen::MatrixXd& cutV,
                                    Eigen::MatrixXi& cutF,
                                    Eigen::MatrixXd& combedField,
                                    Eigen::VectorXi& combedMatching)
  {
    
    using namespace Eigen;
    using namespace std;
    
    //cutting mesh and combing field.
    cut_mesh_with_singularities(wholeV, wholeF, singVertices, intData.face2cut);
    combing(wholeV,wholeF, EV, EF, FE, intData.face2cut, rawField, matching, combedField, combedMatching);
    
    MatrixXi EFi,EH, FH;
    MatrixXd FEs;
    VectorXi VH, HV, HE, HF, nextH, prevH, twinH, innerEdges;
    
    // it stores number of edges per face, for now only tirangular
    VectorXi D = VectorXi::Constant(wholeF.rows(), 3);
    
    // mark vertices as being a singularity vertex of the vector field
    VectorXi isSingular = VectorXi::Zero(wholeV.rows());
    for (int i = 0; i < singVertices.size(); i++)
      isSingular(singVertices(i)) = 1;
    
    //cout<<"singVertices: "<<singVertices<<endl;
    
    intData.constrainedVertices = VectorXi::Zero(wholeV.rows());
    
    //computing extra topological information
    std::vector<int> innerEdgesVec; // collects ids of inner edges
    EFi = Eigen::MatrixXi::Constant(EF.rows(), 2, -1); // number of an edge inside the face
    
    /* used later for internal edges there is 1 or  -1 ie if two faces are adjacent then for a given edge we
     * will have 1 in the frst face and -1 in the second
     */
    FEs = Eigen::MatrixXd::Zero(FE.rows(), FE.cols());
    
    /*
     * here we collect information about position of an edge inside each face containing it. Each triangular face
     * has three edges of ids 0, 1, 2. So EFi(i, k) = j means that the edge i is inside the face k \in [0,1]
     * at the position j.
     */
    for(int i = 0; i < EF.rows(); i++)
    {
      for (int k = 0; k < 2; k++)
      {
        if (EF(i, k) == -1)
          continue;
        for (int j = 0; j < D(EF(i, k)); j++)
          if (FE(EF(i, k), j) == i)
            EFi(i, k) = j;
      }
    }
    
    // collect information about inner edges
    for(int i = 0; i < EF.rows(); i++)
    {
      if(EFi(i, 0) != -1)
        FEs(EF(i, 0), EFi(i, 0)) = 1.0;
      if(EFi(i,1) != -1)
        FEs(EF(i, 1), EFi(i, 1)) = -1.0;
      if ((EF(i, 0) !=-1) && (EF(i,1)!=-1))
        innerEdgesVec.push_back(i);
    }
    
    // copy the information into  Eigen vector
    innerEdges.resize(innerEdgesVec.size());
    for (int i = 0; i < innerEdgesVec.size(); i++)
      innerEdges(i) = innerEdgesVec[i];
    
    // compute the half-edge representation
    hedra::dcel(D, wholeF, EV, EF, EFi, innerEdges, VH, EH, FH, HV, HE, HF, nextH, prevH, twinH);
    
    // find boundary vertices and mark them
    VectorXi isBoundary = VectorXi::Zero(wholeV.rows());
    for (int i = 0; i < HV.rows(); i++)
      if (twinH(i) == -1){
        isBoundary(HV(i)) = 1;
        isSingular(HV(i)) = 0; //boundary vertices cannot be singular
      }
    
    
    
    /*for (int i=0;i<wholeV.rows();i++){
     if (isSingular(i))
     cout<<"vertex "<<i<<" is singular "<<endl;
     }*/
    
    // here we compute a permutation matrix
    vector<MatrixXi> constParmMatrices(intData.N);
    MatrixXi unitPermMatrix = MatrixXi::Zero(intData.N, intData.N);
    for (int i = 0; i < intData.N; i++)
      unitPermMatrix((i + 1) % intData.N, i) = 1;
    
    // generate all the members of the permutation group
    constParmMatrices[0] = MatrixXi::Identity(intData.N, intData.N);
    for (int i = 1; i < intData.N; i++)
      constParmMatrices[i] = unitPermMatrix * constParmMatrices[i - 1];
    
    // each edge which is on the cut seam is marked by 1 and 0 otherwise
    VectorXi isSeam = VectorXi::Zero(EV.rows());
    for(int i = 0; i < FE.rows(); i++)
    {
      for (int j = 0; j < 3; j++)
        if (intData.face2cut(i, j)) // face2cut is initalized by directional::cut_mesh_with_singularities
          isSeam(FE(i, j)) = 1;
    }
    
    // do the same for the half-edges, mark edges which correspond to the cut seam
    VectorXi isHEcut = VectorXi::Zero(HE.rows());
    for(int i = 0; i < wholeF.rows(); i++)
    {
      for (int j = 0; j < 3; j++)
        if (intData.face2cut(i, j)) // face2cut is initalized by directional::cut_mesh_with_singularities
          isHEcut(FH(i, j)) = 1; // FH is face to half-edge mapping
    }
    
    // calculate valency of the vertices which lay on the seam
    VectorXi cutValence = VectorXi::Zero(wholeV.rows());
    for(int i = 0; i < EV.rows(); i++)
    {
      if (isSeam(i))
      {
        cutValence(EV(i, 0))++;
        cutValence(EV(i, 1))++;
      }
    }
    
    
    //establishing transition variables by tracing cut curves
    VectorXi Halfedge2TransitionIndices = VectorXi::Constant(HE.rows(), 32767);
    VectorXi Halfedge2Matching(HE.rows());
    VectorXi isHEClaimed = VectorXi::Zero(HE.rows());
    
    // here we convert the matching that was calculated for the vector field over edges to half-edges
    for (int i = 0; i < HE.rows(); i++)
    {
      // HE is a map between half-edges to edges, but it does not carry the direction
      // EH edge to half-edge mapping
      Halfedge2Matching(i) = (EH(HE(i), 0) == i ? -combedMatching(HE(i)) : combedMatching(HE(i)));
      if(Halfedge2Matching(i) < 0)
        Halfedge2Matching(i) = (intData.N + (Halfedge2Matching(i) % intData.N)) % intData.N;
    }
    
    int currTransition = 1;
    
    /*
     * Next steps: cutting mesh and creating map between wholeF and cutF
     */
    
    //cutting the mesh
    vector<int> cut2whole;
    vector<RowVector3d> cutVlist;
    cutF.resize(wholeF.rows(),3);
    for (int i = 0; i < VH.rows(); i++)
    {
      //creating corners whereever we have non-trivial matching
      int beginH = VH(i);
      int currH = beginH;
      
      //reseting to first cut or first boundary, if exists
      if (!isBoundary(i))
      {
        do
        {
          if (isHEcut(currH)!=0)
            break;
          currH=nextH(twinH(currH));
        } while (beginH!=currH);
      }
      else
      {
        do
        {
          if (twinH(currH)==-1)
            break;
          currH=nextH(twinH(currH));
        } while(twinH(currH)!=-1);
      }
      
      beginH = currH;
      
      do
      {
        if ((isHEcut(currH) != 0) || (beginH == currH))
        {
          cut2whole.push_back(i);
          cutVlist.push_back(wholeV.row(i));
        }
        
        for (int j = 0; j < 3; j++)
          if (wholeF(HF(currH), j) == i)
            cutF(HF(currH), j) = cut2whole.size() - 1;
        currH = twinH(prevH(currH));
      } while((beginH != currH) && (currH != -1));
    }
    
    cutV.resize(cutVlist.size(), 3);
    for(int i = 0; i < cutVlist.size(); i++)
      cutV.row(i) = cutVlist[i];
    
    //starting from each cut-graph node, we trace cut curves
    for(int i = 0;  i < wholeV.rows(); i++)
    {
      if (((cutValence(i) == 2) && (!isSingular(i))) || (cutValence(i) == 0))
        continue;  //either mid-cut curve or non at all
      
      //tracing curves until next node, if not already filled
      int beginH = VH(i);
      
      //reseting to first boundary
      int currH = beginH;
      
      if (isBoundary(i))
      {
        do
        {
          if (twinH(currH) == -1)
            break;
          currH = nextH(twinH(currH));
        } while(twinH(currH) != -1);
      }
      
      beginH = currH;
      
      int nextHalfedgeInCut = -1;
      do
      {
        //unclaimed inner halfedge
        if ((isHEcut(currH) != 0) && (isHEClaimed(currH) == 0) && (twinH(currH) != -1))
        {
          nextHalfedgeInCut = currH;
          Halfedge2TransitionIndices(nextHalfedgeInCut) = currTransition;
          Halfedge2TransitionIndices(twinH(nextHalfedgeInCut)) = -currTransition;
          isHEClaimed(nextHalfedgeInCut) = 1;
          isHEClaimed(twinH(nextHalfedgeInCut)) = 1;
          int nextCutVertex=HV(nextH(nextHalfedgeInCut));
          //advancing on the cut until next node
          while ((cutValence(nextCutVertex) == 2) && (!isSingular(nextCutVertex)) && (!isBoundary(nextCutVertex)))
          {
            int beginH = VH(nextCutVertex);
            int currH = beginH;
            int nextHalfedgeInCut = -1;
            do
            {
              //unclaimed cut halfedge
              if ((isHEcut(currH) != 0) && (isHEClaimed(currH) == 0))
              {
                nextHalfedgeInCut = currH;
                break;
              }
              currH=twinH(prevH(currH));
            } while (beginH != currH);
            Halfedge2TransitionIndices(nextHalfedgeInCut) = currTransition;
            Halfedge2TransitionIndices(twinH(nextHalfedgeInCut)) = -currTransition;
            isHEClaimed(nextHalfedgeInCut) = 1;
            isHEClaimed(twinH(nextHalfedgeInCut)) = 1;
            nextCutVertex = HV(nextH(nextHalfedgeInCut));
          }
          currTransition++;
        }
        currH = twinH(prevH(currH));
      } while((beginH != currH) && (currH != -1));
    }
    // end of cutting
    
    int numTransitions = currTransition - 1;
    //cout<<"numtransitions: "<<numTransitions<<endl;
    vector<Triplet<double> > vertexTrans2CutTriplets, constTriplets;
    vector<Triplet<int> > vertexTrans2CutTripletsInteger, constTripletsInteger;
    //forming the constraints and the singularity positions
    int currConst = 0;
    // this loop set up the transtions (vector field matching) across the cuts
    for (int i = 0; i < VH.rows(); i++)
    {
      std::vector<MatrixXi> permMatrices;
      std::vector<int> permIndices;  //in the space #V + #transitions
      //The initial corner gets the identity without any transition
      permMatrices.push_back(MatrixXi::Identity(intData.N, intData.N));
      permIndices.push_back(i);
      
      int beginH = VH(i);
      int currH = beginH;
      
      //reseting to first cut or boundary, if exists
      if (!isBoundary(i))
      {
        // travel throu the start of the vertex and stop once the edge on the cut is found
        do
        {
          if (isHEcut(currH) != 0)
            break;
          currH = nextH(twinH(currH));
        } while(beginH != currH);
      }
      else
      {
        do
        {
          // travel until an edge without a twin is found, i.e., boundary
          if (twinH(currH) == -1)
            break;
          currH = nextH(twinH(currH));
        } while(twinH(currH) != -1);
      }
      
      // set the beginning to the edge on the cut or on the boundary
      beginH = currH;
      
      int currCutVertex = -1;
      do
      {
        int currFace = HF(currH); // face containing the half-edge
        int newCutVertex = -1;
        //find position of the vertex i in the face of the initial mesh
        for (int j = 0; j < 3; j++)
        {
          if (wholeF(currFace, j) == i)
            newCutVertex = cutF(currFace, j);
        }
        
        //currCorner gets the permutations so far
        if (newCutVertex != currCutVertex)
        {
          currCutVertex = newCutVertex;
          for(int i = 0; i < permIndices.size(); i++)
          {
            // place the perumtation matrix in a bigger matrix, we need to know how things are connected along the cut, no?
            for(int j = 0; j < intData.N; j++)
              for(int k = 0; k < intData.N; k++){
                vertexTrans2CutTriplets.emplace_back(intData.N * currCutVertex + j, intData.N * permIndices[i] + k, (double) permMatrices[i](j, k));
                vertexTrans2CutTripletsInteger.emplace_back(intData.N * currCutVertex + j, intData.N * permIndices[i] + k, permMatrices[i](j, k));
              }
          }
        }
        
        //updating the matrices for the next corner
        int nextHalfedge = twinH(prevH(currH));
        //reached a boundary
        if(nextHalfedge == -1)
        {
          currH = nextHalfedge;
          continue;
        }
        
        // constParmMatrices contains all the members of the permutation group
        MatrixXi nextPermMatrix = constParmMatrices[Halfedge2Matching(nextHalfedge) % intData.N];
        //no update needed
        if(isHEcut(nextHalfedge) == 0)
        {
          currH = nextHalfedge;
          continue;
        }
        
        //otherwise, updating matrices with transition
        int nextTransition = Halfedge2TransitionIndices(nextHalfedge);
        //Pe*f + Je
        if(nextTransition > 0)
        {
          for(int j = 0; j < permMatrices.size(); j++)
            permMatrices[j] = nextPermMatrix * permMatrices[j];
          
          //and identity on the fresh transition
          permMatrices.push_back(MatrixXi::Identity(intData.N, intData.N));
          permIndices.push_back(wholeV.rows() + nextTransition - 1);
        }
        // (Pe*(f-Je))  matrix is already inverse since halfedge matching is minused
        else
        {
          //reverse order
          permMatrices.push_back(-MatrixXi::Identity(intData.N, intData.N));
          permIndices.push_back(wholeV.rows() - nextTransition - 1);
          
          for(int j = 0; j < permMatrices.size(); j++)
            permMatrices[j] = nextPermMatrix * permMatrices[j];
        }
        currH = nextHalfedge;
      } while((currH != beginH) && (currH != -1));
      
      //cleaning parmMatrices and permIndices to see if there is a constraint or reveal singularity-from-transition
      std::set<int> cleanPermIndicesSet(permIndices.begin(), permIndices.end());
      std::vector<int> cleanPermIndices(cleanPermIndicesSet.begin(), cleanPermIndicesSet.end());
      std::vector<MatrixXi> cleanPermMatrices(cleanPermIndices.size());
      
      for (int j = 0; j < cleanPermIndices.size(); j++)
      {
        cleanPermMatrices[j] = MatrixXi::Zero(intData.N, intData.N);
        for(int k = 0;k < permIndices.size(); k++)
          if(cleanPermIndices[j] == permIndices[k])
            cleanPermMatrices[j] += permMatrices[k];
        if(cleanPermIndices[j] == i)
          cleanPermMatrices[j] -= MatrixXi::Identity(intData.N, intData.N);
      }
      
      //if not all matrices are zero, there is a constraint
      bool isConstraint = false;
      for(int j = 0; j < cleanPermMatrices.size(); j++)
        if (cleanPermMatrices[j].cwiseAbs().maxCoeff() != 0)
          isConstraint = true;
      
      if((isConstraint) && (!isBoundary(i)))
      {
        for(int j = 0; j < cleanPermMatrices.size(); j++)
        {
          for(int k = 0; k < intData.N; k++)
            for(int l = 0; l < intData.N; l++){
              constTriplets.emplace_back(intData.N * currConst + k, intData.N * cleanPermIndices[j] + l, (double) cleanPermMatrices[j](k, l));
              constTripletsInteger.emplace_back(intData.N * currConst + k, intData.N * cleanPermIndices[j] + l, cleanPermMatrices[j](k, l));
            }
        }
        currConst++;
        intData.constrainedVertices(i) = 1;
      }
    }
    
    vector< Triplet< double > > cleanTriplets;
    vector< Triplet< int > > cleanTripletsInteger;
    
    intData.vertexTrans2CutMat.resize(intData.N * cutV.rows(), intData.N * (wholeV.rows() + numTransitions));
    intData.vertexTrans2CutMatInteger.resize(intData.N * cutV.rows(), intData.N * (wholeV.rows() + numTransitions));
    cleanTriplets.clear();
    cleanTripletsInteger.clear();
    for(int i = 0; i < vertexTrans2CutTriplets.size(); i++){
      if(vertexTrans2CutTripletsInteger[i].value() != 0){
        cleanTripletsInteger.push_back(vertexTrans2CutTripletsInteger[i]);
        cleanTriplets.push_back(vertexTrans2CutTriplets[i]);
      }
      // if(std::abs((float)vertexTrans2CutTriplets[i].value())>10e-7)
    }
    intData.vertexTrans2CutMat.setFromTriplets(cleanTriplets.begin(), cleanTriplets.end());
    intData.vertexTrans2CutMatInteger.setFromTriplets(cleanTripletsInteger.begin(), cleanTripletsInteger.end());
    
    //
    
    intData.constraintMat.resize(intData.N * currConst, intData.N * (wholeV.rows() + numTransitions));
    intData.constraintMatInteger.resize(intData.N * currConst, intData.N * (wholeV.rows() + numTransitions));
    cleanTriplets.clear();
    cleanTripletsInteger.clear();
    for(int i = 0; i < constTriplets.size(); i++){
      if(constTripletsInteger[i].value() != 0){
        cleanTripletsInteger.push_back(constTripletsInteger[i]);
        cleanTriplets.push_back(constTriplets[i]);
      }
    }
    intData.constraintMat.setFromTriplets(cleanTriplets.begin(), cleanTriplets.end());
    intData.constraintMatInteger.setFromTriplets(cleanTripletsInteger.begin(), cleanTripletsInteger.end());
    
    //doing the integer spanning matrix
    intData.intSpanMat.resize(intData.n * (wholeV.rows() + numTransitions), intData.n * (wholeV.rows() + numTransitions));
    intData.intSpanMatInteger.resize(intData.n * (wholeV.rows() + numTransitions), intData.n * (wholeV.rows() + numTransitions));
    vector<Triplet<double> > intSpanMatTriplets;
    vector<Triplet<int> > intSpanMatTripletsInteger;
    for (int i=0;i<intData.n*numTransitions;i+=intData.n){
      for(int k = 0; k < intData.n; k++)
        for(int l = 0; l < intData.n; l++){
          if (intData.periodMat(k,l)!=0){
            intSpanMatTriplets.emplace_back(intData.n * wholeV.rows()+i+k, intData.n * wholeV.rows()+i+l, (double)intData.periodMat(k,l));
            intSpanMatTripletsInteger.emplace_back(intData.n * wholeV.rows()+i+k, intData.n * wholeV.rows()+i+l, intData.periodMat(k,l));
          }
        }
    }
    for (int i=0;i<intData.n * wholeV.rows();i++){
      intSpanMatTriplets.emplace_back(i,i,1.0);
      intSpanMatTripletsInteger.emplace_back(i,i,1);
    }
    
    intData.intSpanMat.setFromTriplets(intSpanMatTriplets.begin(), intSpanMatTriplets.end());
    intData.intSpanMatInteger.setFromTriplets(intSpanMatTripletsInteger.begin(), intSpanMatTripletsInteger.end());
    
    //filtering out barycentric symmetry, including sign symmetry. The parameterization should always only include n dof for the surface
    //TODO: this assumes n divides N!
    intData.linRedMat.resize(intData.N * (wholeV.rows() + numTransitions), intData.n * (wholeV.rows() + numTransitions));
    intData.linRedMatInteger.resize(intData.N * (wholeV.rows() + numTransitions), intData.n * (wholeV.rows() + numTransitions));
    vector<Triplet<double> > linRedMatTriplets;
    vector<Triplet<int> > linRedMatTripletsInteger;
    for(int i = 0; i < intData.N*(wholeV.rows() + numTransitions); i +=intData.N)
      for(int k = 0; k < intData.N; k++)
        for(int l = 0; l < intData.n; l++){
          if (intData.linRed(k,l)!=0){
            linRedMatTriplets.emplace_back(i + k, i*intData.n/intData.N + l, (double)intData.linRed(k,l));
            linRedMatTripletsInteger.emplace_back(i + k, i*intData.n/intData.N + l, intData.linRed(k,l));
          }
        }
    
    intData.linRedMat.setFromTriplets(linRedMatTriplets.begin(), linRedMatTriplets.end());
    intData.linRedMatInteger.setFromTriplets(linRedMatTripletsInteger.begin(), linRedMatTripletsInteger.end());
    
    //integer variables are per single "d" packet, and the rounding is done for the N functions with projection over linRed
    intData.integerVars.resize(numTransitions);
    intData.integerVars.setZero();
    for(int i = 0; i < numTransitions; i++)
      intData.integerVars(i) = wholeV.rows() + i;
    
    //fixed values
    intData.fixedIndices.resize(intData.n);
    if (isSingular.sum()==0){  //no inner singular vertices; vertex 0 is set to (0....0)
      for (int j=0;j<intData.n;j++)
        intData.fixedIndices(j)=j;
    }else {  //fixing first singularity to (0.5,....0.5)
      int firstSing;
      for (firstSing=0;firstSing<isSingular.size();firstSing++)
        if (isSingular(firstSing))
          break;
      //firstSing=0; //like before
      for (int j=0;j<intData.n;j++)
        intData.fixedIndices(j)=intData.n*firstSing+j;
    }
    
    //creating list of singular corners and singular integer matrix
    VectorXi singularIndices(intData.n * isSingular.sum());
    int counter=0;
    for (int i=0;i<isSingular.size();i++){
      if (isSingular(i))
        for (int j=0;j<intData.n;j++)
          singularIndices(counter++)=intData.n*i+j;
    }
    
    //doing the integer spanning matrix
    intData.singIntSpanMat.resize(intData.n * (wholeV.rows() + numTransitions), intData.n * (wholeV.rows() + numTransitions));
    intData.singIntSpanMatInteger.resize(intData.n * (wholeV.rows() + numTransitions), intData.n * (wholeV.rows() + numTransitions));
    vector<Triplet<double> > singIntSpanMatTriplets;
    vector<Triplet<int> > singIntSpanMatTripletsInteger;
    for (int i=0;i<isSingular.size();i++){
      if (!isSingular(i)){
        for (int j=0;j<intData.n;j++){
          singIntSpanMatTriplets.emplace_back(intData.n * i+j,intData.n * i+j,1.0);
          singIntSpanMatTripletsInteger.emplace_back(intData.n * i+j,intData.n * i+j,1);
        }
      } else {
        for(int k = 0; k < intData.n; k++)
          for(int l = 0; l < intData.n; l++){
            if (intData.periodMat(k,l)!=0){
              singIntSpanMatTriplets.emplace_back(intData.n*i+k, intData.n*i+l, (double)intData.periodMat(k,l));
              singIntSpanMatTripletsInteger.emplace_back(intData.n*i+k, intData.n*i+l, intData.periodMat(k,l));
            }
          }
      }
    }
    
    for (int i=intData.n * wholeV.rows() ; i<intData.n*(wholeV.rows()+numTransitions);i++){
      singIntSpanMatTriplets.emplace_back(i,i,1.0);
      singIntSpanMatTripletsInteger.emplace_back(i,i,1);
    }
    
    intData.singIntSpanMat.setFromTriplets(singIntSpanMatTriplets.begin(), singIntSpanMatTriplets.end());
    intData.singIntSpanMatInteger.setFromTriplets(singIntSpanMatTripletsInteger.begin(), singIntSpanMatTripletsInteger.end());
    
    intData.singularIndices=singularIndices;
    intData.fixedValues.resize(intData.n);
    intData.fixedValues.setConstant(0);
    
  }
}

#endif


