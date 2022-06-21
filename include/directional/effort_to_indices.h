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
#include <directional/CartesianField.h>


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
  IGL_INLINE void effort_to_indices(directional::CartesianField& field)
  {
    Eigen::VectorXd effortInner(field.innerAdjacencies.size());
    for (int i=0;i<field.innerAdjacencies.size();i++)
      effortInner(i)=field.effort(field.innerAdjacencies(i));
    Eigen::VectorXi fullIndices;
    directional::effort_to_indices(field.dualCycles, effortInner, field.matching, field.cycleCurvatures, field.N, fullIndices);
   
    Eigen::VectorXi indices(field.element2Cycle.size());
    for (int i=0;i<field.element2Cycle.size();i++)
      indices(i)=fullIndices(field.element2Cycle(i));
                        
    std::vector<int> singCyclesList;
    std::vector<int> singIndicesList;
    for (int i=0;i<field.element2Cycle.size();i++)
      if (indices(i)!=0){
        singCyclesList.push_back(i);
        singIndicesList.push_back(indices(i));
      }
    
    Eigen::VectorXi singCycles(singCyclesList.size());
    Eigen::VectorXi singIndices(singIndicesList.size());
    for (int i=0;i<singCyclesList.size();i++){
      singCycles(i)=singCyclesList[i];
      singIndices(i)=singIndicesList[i];
    }
    field.set_singularities(singCycles, singIndices);
  }
}

#endif


