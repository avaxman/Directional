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
      std::queue<int> HEQueue;
      Eigen::VectorXi isHECut = Eigen::VectorXi::Zero(mesh.HV.rows());
      Eigen::VectorXi isFaceVisited = Eigen::VectorXi::Zero(mesh.F.rows());
      std::vector<int> cutVertices;
      int currFace = 0;
      for (int i=0;i<3;i++)
        HEQueue.push(mesh.FH(currFace, i));
      isFaceVisited(0)=1;
      while (!HEQueue.empty()){
          int currHE = HEQueue.front();
          HEQueue.pop();
          if ((isHECut(currHE))||(mesh.twinH(currHE)==-1))
              continue;

          if (isFaceVisited(mesh.HF(mesh.twinH(currHE)))){ //two fronts meeting
              isHECut(currHE)=1;
              isHECut(mesh.twinH(currHE))=1;
              cutVertices.push_back(mesh.HV(currHE));
              cutVertices.push_back(mesh.HV(mesh.twinH(currHE)));
              continue;
          }

          //Otherwise flood into the next face
          HEQueue.push(mesh.nextH(mesh.twinH(currHE)));
          HEQueue.push(mesh.nextH(mesh.nextH(mesh.twinH(currHE))));
      }

      //Connecting all singularities to the cut graph
      for (int i=0;i<singularities.size();i++) {
          std::vector<int> HEPath;
          shortest_path(mesh, singularities(i), cutVertices, HEPath);
          for (int j=0;j<HEPath.size();j++) {
              isHECut[HEPath[j]] = 1;
              if (mesh.twinH(HEPath[j]) != -1)
                  isHECut[mesh.twinH(HEPath[j])] = 1;
          }
      }

      face2cut.resize(mesh.F.rows(),3);
      for (int i=0;i<mesh.F.rows();i++)
          for (int j=0;j<3;j++)
              face2cut(i,j)=isHECut(mesh.FH(i,j));

  }

};


#endif
