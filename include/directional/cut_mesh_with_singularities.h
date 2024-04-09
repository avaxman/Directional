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
  IGL_INLINE void directional::cut_mesh_with_singularities(const TriMesh& mesh,
                                                           const Eigen::VectorXi &singularities,
                                                           Eigen::MatrixXi &cuts)
  {

      //doing a flood-fill to cut mesh into a topological discs
      std::queue<int> HEQueue;
      Eigen::VectorXi isHECut = Eigen::VectorXi::Zero(mesh.HV.rows());
      Eigen::VectorXi isFaceVisited = Eigen::VectorXi::Zero(mesh.F.rows());
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
              continue;
          }

          //Otherwise flood into the next face
          HEQueue.push(mesh.nextH(mesh.twinH(currHE)));
          HEQueue.push(mesh.nextH(mesh.nextH(mesh.twinH(currHE))));
      }

      //Connecting all singularities to the cut graph



      //first, get a spanning tree for the mesh (no missmatch needed)
      /*igl::cut_mesh_from_singularities(V, F, Eigen::MatrixXd::Zero(F.rows(), 3).eval(), cuts);

      std::set<int> vertices_in_cut;
      for (int i =0; i< cuts.rows(); ++i)
          for (int j =0;j< cuts.cols(); ++j)
              if (cuts(i,j))
                  vertices_in_cut.insert(F(i,j));*/

      //then, add all singularities one by one by using Dijkstra's algorithm
      for (int i = 0; i<singularities.rows(); ++i)
      {
          std::vector<int> path;
          Eigen::VectorXd min_distance;
          Eigen::VectorXi previous;
          int vertex_found = igl::dijkstra(singularities[i], vertices_in_cut, VV, min_distance, previous);
          if(vertex_found ==-1)
              // this means that there are no cuts
              path.push_back(singularities[i]);
          else
              igl::dijkstra(vertex_found, previous, path);

          vertices_in_cut.insert(path.begin(), path.end());

          //insert to cut
          for (int ii = 0; ii<path.size()-1; ++ii)
          {
              const int &v0 = path[ii];
              const int &v1 = path[ii+1];

              std::vector<int> vf0 = VF[v0];
              std::sort(vf0.begin(), vf0.end());
              std::vector<int> vf1 = VF[v1];
              std::sort(vf1.begin(), vf1.end());
              std::vector<int> common_face_v(std::max(vf0.size(),vf1.size()));
              std::vector<int>::iterator it;
              it=std::set_intersection (vf0.begin(), vf0.end(), vf1.begin(), vf1.end(), common_face_v.begin());
              common_face_v.resize(it-common_face_v.begin());
              assert(common_face_v.size() == 2);

              const int &fi = common_face_v[0];
              int j=-1;
              for (unsigned z=0; z<3; ++z)
                  if (((F(fi,z) == v0) && (F(fi,(z+1)%3) == v1)) ||((F(fi,z) == v1) && (F(fi,(z+1)%3) == v0)))
                      j=z;
              assert(j!=-1);
              cuts(fi,j) = 1;
              cuts(TT(fi,j), TTi(fi,j)) = 1;

          }
      }

  }

  
};


#endif
