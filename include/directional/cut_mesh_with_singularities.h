// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2015 Olga Diamanti <olga.diam@gmail.com>, 2018 Amir vaxman
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_CUT_MESH_WITH_SINGULARITIES
#define DIRECTIONAL_CUT_MESH_WITH_SINGULARITIES

#include <Eigen/Core>
#include <vector>
#include <set>
#include <directional/shortest_path.h>

namespace directional {
  // Given a mesh and the singularities of a polyvector field, cut the mesh
  // to disk topology in such a way that the singularities lie at the boundary of
  // the disk, as described in the paper "Mixed Integer Quadrangulation" by
  // Bommes et al. 2009.
  // Inputs:
  //   V                #V by 3 list of the vertex positions
  //   F                #F by 3 list of the faces (must be triangles)
  //   VF               #V list of lists of incident faces (adjacency list), e.g.
  //                    as returned by igl::vertex_triangle_adjacency
  //   VV               #V list of lists of incident vertices (adjacency list), e.g.
  //                    as returned by igl::adjacency_list
  //   TT               #F by 3 triangle to triangle adjacent matrix (e.g. computed
  //                    via igl:triangle_triangle_adjacency)
  //   TTi              #F by 3 adjacent matrix, the element i,j is the id of edge of the
  //                    triangle TT(i,j) that is adjacent with triangle i (e.g. computed
  //                    via igl:triangle_triangle_adjacency)
  //   singularities    #S by 1 list of the indices of the singular vertices
  // Outputs:
  //   cuts             #F by 3 list of boolean flags, indicating the edges that need to be cut
  //                    (has 1 at the face edges that are to be cut, 0 otherwise)
  //
  inline void cut_mesh_with_singularities(const TriMesh& mesh,
                                          const Eigen::VectorXi &singularities,
                                          Eigen::MatrixXi &face2cut)
  {

      //doing a flood-fill to cut mesh into a topological discs
      std::queue<std::pair<int, int>> faceQueue;
      Eigen::VectorXi isHECut = Eigen::VectorXi::Ones(mesh.HV.rows());
      Eigen::VectorXi isFaceVisited = Eigen::VectorXi::Zero(mesh.F.rows());
      std::vector<int> cutVertices;
      int currFace = 0;
      faceQueue.push(std::pair<int,int>(currFace,-1));
      while (!faceQueue.empty()){
          std::pair<int, int> currFaceHE = faceQueue.front();
          int currFace =currFaceHE.first;
          int prevHE = currFaceHE.second;
          faceQueue.pop();
          if (isFaceVisited[currFace]==1)
              continue;

          isFaceVisited[currFace]=1;
          if (prevHE!=-1) {
              isHECut(prevHE) = 0;
              isHECut(mesh.twinH(prevHE)) = 0;
          }
          for (int i=0;i<3;i++){
              int currHE = mesh.FH(currFace,i);
              if (mesh.twinH(currHE)==-1){
                  isHECut(currHE)=0;
                  continue;
              }
              int nextFace = mesh.HF(mesh.twinH(currHE));
              if (!isFaceVisited[nextFace]) //can spill into that face
                  faceQueue.push(std::pair<int, int>(nextFace, currHE));
          }

      }

      //retract valence 1 vertices
      Eigen::VectorXi cutValences = Eigen::VectorXi::Zero(mesh.V.rows());

      //gathering cut vertices
      Eigen::VectorXi isSingularity = Eigen::VectorXi::Zero(mesh.V.rows());
      for (int i=0;i<singularities.size();i++)
          isSingularity(singularities(i))=1;
      for (int i=0;i<mesh.HV.size();i++) {
          if ((isHECut(i))||(mesh.twinH(i)==-1))  //if cut or a boundary
              cutValences(mesh.HV(i))++;  //the twin should already be inside
      }

      std::queue<int> cutQueue;
      for (int i=0;i<mesh.HV.rows();i++)
          if ((isHECut(i))&&(cutValences(mesh.HV(i))==1)&&(!isSingularity(mesh.HV(i))))
              cutQueue.push(i);

      //std::cout<<"isHECut.sum(): "<<isHECut.sum()<<std::endl;
      //int stop = 3000;
      while (!cutQueue.empty()){
          /*stop--;
          if (stop==2000)
              break;*/
          int currHE = cutQueue.front();
          cutQueue.pop();
          if (!isHECut(currHE))
              continue;
          if (cutValences(mesh.HV(currHE))!=1)
              continue;
          if (isSingularity(mesh.HV(currHE)))
              continue;

          isHECut(currHE)=0;
          isHECut(mesh.twinH(currHE))=0;
          //finding the next edge
          int nextVertex = mesh.HV(mesh.nextH(currHE));
          if (mesh.isBoundaryVertex(nextVertex))
            continue;
          cutValences(nextVertex)--;
          cutValences(mesh.HV(currHE))--;
          if (cutValences(nextVertex)==1) {  //finding next edge
              int hebegin = mesh.VH(nextVertex);
              int hecurr=hebegin;
              do{
                  if (isHECut(hecurr))
                      break;
                  hecurr = mesh.twinH(mesh.prevH(hecurr));
              }while (hecurr!=hebegin);
              assert("hecurr is not cut!" && isHECut(hecurr));
              cutQueue.push(hecurr);
          }
      }

      //std::cout<<"isHECut.sum(): "<<isHECut.sum()<<std::endl;

      //Connecting all singularities to the cut graph
      /*for (int i=0;i<singularities.size();i++) {
          std::vector<int> HEPath;
         // shortest_path(mesh, singularities(i), cutVertices, HEPath);
          for (int j=0;j<HEPath.size();j++) {
              //isHECut[HEPath[j]] = 1;
              //if (mesh.twinH(HEPath[j]) != -1)
              //    isHECut[mesh.twinH(HEPath[j])] = 1;
          }
      }*/

      face2cut.resize(mesh.F.rows(),3);
      for (int i=0;i<mesh.F.rows();i++)
          for (int j=0;j<3;j++)
              face2cut(i,j)=isHECut(mesh.FH(i,j));

  }

};


#endif
