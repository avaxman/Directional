// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_SETUP_PARAMETERIZATION_H
#define DIRECTIONAL_SETUP_PARAMETERIZATION_H

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

namespace directional
{
  struct ParameterizationData // writing this as a class in C++ is a bit strange struct and class are the same only the default access is opposite
  {
    Eigen::SparseMatrix<double> vertexTrans2CutMat;
    Eigen::SparseMatrix<double> constraintMat;
    Eigen::SparseMatrix<double> symmMat;
    Eigen::VectorXi constrainedVertices;
    Eigen::VectorXi integerVars;
    Eigen::MatrixXi face2cut;
    
    ParameterizationData(){}
    ~ParameterizationData(){}
  };
  
  
  // Setting up the seamless parameterization algorithm
  // Input:
  //  N:          The degree of the field.
  //  wholeV:     #V x 3 vertex coordinates
  //  wholeF:     #F x 3 face vertex indices
  //  EV:         #E x 2 edges to vertices indices
  //  EF:         #E x 2 edges to faces indices
  // matching:    #E matching function, where vector k in EF(i,0) matches to vector (k+matching(k))%N in EF(i,1). In case of boundary, there is a -1. Most matching should be zero due to prior combing.
  // singVertices:list of singular vertices in wholeV.
  // Output:
  //  pd:         parameterization data subsequently used in directional::parameterize();
  //  cutV:       the Vertices of the cut mesh.
  //  cutF:       The Vaces of the cut mesh.
  
  IGL_INLINE void setup_parameterization(const int N,
                                         const Eigen::MatrixXd& wholeV,
                                         const Eigen::MatrixXi& wholeF,
                                         const Eigen::MatrixXi& EV,
                                         const Eigen::MatrixXi& EF,
                                         const Eigen::MatrixXi& FE,
                                         const Eigen::VectorXi& matching,
                                         const Eigen::VectorXi& singVertices,
                                         ParameterizationData& pd,
                                         Eigen::MatrixXd& cutV,
                                         Eigen::MatrixXi& cutF)
  {
    
    using namespace Eigen;
    using namespace std;
    
    MatrixXi EFi,EH, FH;
    MatrixXd FEs;
    VectorXi VH, HV, HE, HF, nextH, prevH, twinH, innerEdges;

    // it stores number of edges per face, for now only tirangular
    VectorXi D = VectorXi::Constant(wholeF.rows(), 3);

    // mark vertices as being a singularity vertex of the vector field
    VectorXi isSingular = VectorXi::Zero(wholeV.rows());
    for (int i = 0; i < singVertices.size(); i++)
      isSingular(singVertices(i)) = 1;
    
    pd.constrainedVertices = VectorXi::Zero(wholeV.rows());
    
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
      if (twinH(i) == -1)
        isBoundary(HV(i)) = 1;

    // here we compute a permutation matrix
    vector<MatrixXi> constParmMatrices(N);
    MatrixXi unitPermMatrix = MatrixXi::Zero(N, N);
    for (int i = 0; i < N; i++)
      unitPermMatrix((i + 1) % N, i) = 1;

    // generate all the members of the permutation group
    constParmMatrices[0] = MatrixXi::Identity(N, N);
    for (int i = 1; i < N; i++)
      constParmMatrices[i] = unitPermMatrix * constParmMatrices[i - 1];

    // each edge which is on the cut seam is marked by 1 and 0 otherwise
    VectorXi isSeam = VectorXi::Zero(EV.rows());
    for(int i = 0; i < FE.rows(); i++)
    {
      for (int j = 0; j < 3; j++)
        if (pd.face2cut(i, j)) // face2cut is initalized by directional::cut_mesh_with_singularities
          isSeam(FE(i, j)) = 1;
    }

    // do the same for the half-edges, mark edges which correspond to the cut seam
    VectorXi isHEcut = VectorXi::Zero(HE.rows());
    for(int i = 0; i < wholeF.rows(); i++)
    {
      for (int j = 0; j < 3; j++)
        if (pd.face2cut(i, j)) // face2cut is initalized by directional::cut_mesh_with_singularities
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
      Halfedge2Matching(i) = (EH(HE(i), 0) == i ? -matching(HE(i)) : matching(HE(i)));
      if(Halfedge2Matching(i) < 0)
        Halfedge2Matching(i) = (N + (Halfedge2Matching(i) % N)) % N;
    }
    
    int currTransition = 1;
    
    /*
    * Next steps: cutting mesh and creating map between wholeF and cutF
    */

    //cutting the mesh
    vector<int> cut2whole;
    vector<RowVector3d> cutVlist;
    cutF.conservativeResize(wholeF.rows(),3);
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
    
    cutV.conservativeResize(cutVlist.size(), 3);
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
    vector<Triplet<double> > vertexTrans2CutTriplets, constTriplets;
    //forming the constraints and the singularity positions
    int currConst = 0;
    // this loop set up the transtions (vector field matching) across the cuts
    for (int i = 0; i < VH.rows(); i++)
    {
      std::vector<MatrixXi> permMatrices;
      std::vector<int> permIndices;  //in the space #V + #transitions
      //The initial corner gets the identity without any transition
      permMatrices.push_back(MatrixXi::Identity(N, N));
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
            for(int j = 0; j < N; j++)
              for(int k = 0; k < N; k++)
                vertexTrans2CutTriplets.emplace_back(N * currCutVertex + j, N * permIndices[i] + k, (double) permMatrices[i](j, k));
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
        MatrixXi nextPermMatrix = constParmMatrices[Halfedge2Matching(nextHalfedge) % N];
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
          permMatrices.push_back(MatrixXi::Identity(N, N));
          permIndices.push_back(wholeV.rows() + nextTransition - 1);
        }
        // (Pe*(f-Je))  matrix is already inverse since halfedge matching is minused
        else
        {
          //reverse order
          permMatrices.push_back(-MatrixXi::Identity(N, N));
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
        cleanPermMatrices[j] = MatrixXi::Zero(N, N);
        for(int k = 0;k < permIndices.size(); k++)
          if(cleanPermIndices[j] == permIndices[k])
            cleanPermMatrices[j] += permMatrices[k];
        if(cleanPermIndices[j] == i)
          cleanPermMatrices[j] -= MatrixXi::Identity(N, N);
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
          for(int k = 0; k < N; k++)
            for(int l = 0; l < N; l++)
              constTriplets.emplace_back(N * currConst + k, N * cleanPermIndices[j] + l, (double) cleanPermMatrices[j](k, l));
        }
        currConst++;
        pd.constrainedVertices(i) = 1;
      }
    }

    vector< Triplet< double > > cleanTriplets;
   
    pd.vertexTrans2CutMat.conservativeResize(N * cutV.rows(), N * (wholeV.rows() + numTransitions));
    cleanTriplets.clear();
    for(int i = 0; i < vertexTrans2CutTriplets.size(); i++)
      if((float)vertexTrans2CutTriplets[i].value() != 0.0f)
        cleanTriplets.push_back(vertexTrans2CutTriplets[i]);
    pd.vertexTrans2CutMat.setFromTriplets(cleanTriplets.begin(), cleanTriplets.end());
    
    pd.constraintMat.conservativeResize(N * currConst, N * (wholeV.rows() + numTransitions));
    cleanTriplets.clear();
    for(int i = 0; i < constTriplets.size(); i++)
      if((float)constTriplets[i].value() != 0.0f)
        cleanTriplets.push_back(constTriplets[i]);
    pd.constraintMat.setFromTriplets(cleanTriplets.begin(), cleanTriplets.end());
    
    //filtering out sign symmetry
    // symmetry between (u, v, w) and (-u, -v, -w)
    pd.symmMat.conservativeResize(N * (wholeV.rows() + numTransitions), (int)(N * (wholeV.rows() + numTransitions) / 2.));
    vector<Triplet<double> > symmMatTriplets;
    for(int i = 0; i < N * (wholeV.rows() + numTransitions); i += N)
    {
     for(int j = 0; j < N / 2.; j++)
      {
        symmMatTriplets.emplace_back(i + j, (int)(i / 2. + j), 1.0);
        symmMatTriplets.emplace_back(i + j + (int)(N / 2.), (int)(i / 2. + j), -1.0);
      }
    }
    pd.symmMat.setFromTriplets(symmMatTriplets.begin(), symmMatTriplets.end());
    
    
    //in this case, also doing UV->UVW packing. This only works for N=6.
    if(N == 6)
    {
      /* this is just conversion between axial and cube coordinate systems, with the exception that
       * normally in the axial corrdianetes r is z and q is x then y is -r -q
       * In fact the plane has here the equation x - y + z = 0 and NOT as usually x + y + z = 0,
       * therefore q is x and r is y. More information can be found here: https://www.redblobgames.com/grids/hexagons/
       */
      SparseMatrix<double> baryMat(3 * (wholeV.rows() + numTransitions), 2 * (wholeV.rows() + numTransitions));
      vector<Triplet<double> > baryMatTriplets;
      for(int i = 0; i < 3 * (wholeV.rows() + numTransitions); i += 3)
      {
        baryMatTriplets.emplace_back(i, (int)((i * 2.) / 3.), 1.0);
        baryMatTriplets.emplace_back(i + 1, (int)((i * 2.) / 3. + 1.), 1.0);

        baryMatTriplets.emplace_back(i + 2, (int)((i * 2.) / 3.), -1.0);
        baryMatTriplets.emplace_back(i + 2, (int)((i * 2.) / 3. + 1), 1.0);
      }
      baryMat.setFromTriplets(baryMatTriplets.begin(), baryMatTriplets.end());
      pd.symmMat = pd.symmMat * baryMat;
      
      pd.integerVars.conservativeResize(2 * numTransitions);
      pd.integerVars.setZero();
      for(int i = 0; i < numTransitions; i++)
      {
        for (int j = 0; j < 2; j++)
          pd.integerVars(2 * i + j) = 2 * (wholeV.rows() + i) + j;
      }
    }
    else
    {
      pd.integerVars.conservativeResize((int)(N / 2. * numTransitions));
      pd.integerVars.setZero();
      for(int i = 0; i < numTransitions; i++)
      {
        for (int j = 0; j < N / 2.; j++)
          pd.integerVars((int)(N / 2. * i) + j) = (int)(N / 2.) * (wholeV.rows() + i) + j;

      }
    }
  }
}

#endif


