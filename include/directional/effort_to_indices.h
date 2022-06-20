//This file is part of Directional, a library for directional field processing.
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_EFFORT_TO_INDICES_H
#define DIRECTIONAL_EFFORT_TO_INDICES_H

#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <igl/igl_inline.h>
#include <igl/boundary_loop.h>
#include <directional/dual_cycles.h>
#include <directional/TriMesh.h>
#include <directional/FaceField.h>


namespace directional
{
  // Computes cycle-based indices from dual-edge-based efforts.
  // Note: input is effort (sum of rotation angles), and not individual rotation angles
  // Input:
  //  basisCycles:    #c by #iE (inner edges of the mesh) the oriented basis cycles around which the indices are measured
  //  effort:         #iE the effort (sum of rotation angles) of matched vectors across the dual edge. Equal to N*rotation angles for N-RoSy fields.
  //  cycleCurvature: #c the cycle curvature (for instance, from directional::dual_cycles)
  //  N:              The degree of the field
  // Output:
  //  indices:     #c the index of the cycle x N (always an integer).
  IGL_INLINE void effort_to_indices(const Eigen::SparseMatrix<double>& basisCycles,
                                    const Eigen::VectorXd& effort,
                                    const Eigen::VectorXi& matching,
                                    const Eigen::VectorXd& cycleCurvature,
                                    const int N,
                                    Eigen::VectorXi& indices)
  {
    using namespace std;
    Eigen::VectorXd dIndices = ((basisCycles * effort + N*cycleCurvature).array() / (2.0*igl::PI));  //this should already be an integer up to numerical precision
        
    indices.conservativeResize(dIndices.size());
    for (int i=0;i<indices.size();i++)
      indices(i)=std::round(dIndices(i));

  }
  
  
  // minimal version without precomputed cycles or inner edges, returning only inner-vertex singularities
  IGL_INLINE void effort_to_indices(directional::FaceField& field)
  {
    Eigen::SparseMatrix<double> basisCycles;
    Eigen::VectorXd cycleCurvature;
    Eigen::VectorXi vertex2cycle;
    Eigen::VectorXi innerEdges;
    directional::dual_cycles(field.mesh->V, field.mesh->F,field.mesh->EV, field.mesh->EF, basisCycles, cycleCurvature, vertex2cycle, innerEdges);
    Eigen::VectorXd effortInner(innerEdges.size());
    for (int i=0;i<innerEdges.size();i++)
      effortInner(i)=field.effort(innerEdges(i));
    Eigen::VectorXi fullIndices;
    directional::effort_to_indices(basisCycles, effortInner, field.matching, cycleCurvature, field.N, fullIndices);
   
    Eigen::VectorXi indices(field.mesh->V.rows());
    for (int i=0;i<field.mesh->V.rows();i++)
      indices(i)=fullIndices(vertex2cycle(i));
    
    //removing boundary indices
   std::vector<std::vector<int> > L;
    igl::boundary_loop(field.mesh->F, L);
    for (int j=0;j<L.size();j++)
      for (int k=0;k<L[j].size();k++)
        indices(L[j][k])=0;
                    
    std::vector<int> singVerticesList;
    std::vector<int> singIndicesList;
    for (int i=0;i<field.mesh->V.rows();i++)
      if (indices(i)!=0){
        singVerticesList.push_back(i);
        singIndicesList.push_back(indices(i));
      }
    
    field.singVertices.resize(singVerticesList.size());
    field.singIndices.resize(singIndicesList.size());
    for (int i=0;i<singVerticesList.size();i++){
      field.singVertices(i)=singVerticesList[i];
      field.singIndices(i)=singIndicesList[i];
    }
  }
}

#endif


