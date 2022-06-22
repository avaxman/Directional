// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_INDEX_PRESCRIPTION_H
#define DIRECTIONAL_INDEX_PRESCRIPTION_H

#include <Eigen/Core>
#include <vector>
#include <cmath>
#include <igl/igl_inline.h>
#include <directional/CartesianField.h>
#include <directional/rotation_to_raw.h>


namespace directional
{    
  // Computes the dual-edge-based rotation angles that are required to reproduce a prescribed set of indices on the dual cycles of the mesh.
  // In case the sum of curvature is not consistent with the topology, the system is solved in least squares and unexpected singularities may appear elsewhere. linfError will mostl like be far from zero.
  // Inputs:
  //  V:          #V by 3 vertex coordinates
  //  F:          #F by 3 face vertex indices
  //  EV:         #E by 3 edges
  //  innerEdges: #iE the subset from EV of inner (non-boundary) edges.
  //  basisCycles:#c X #E the basis cycles matrix (obtained from directional::dual_cycles
  //  indices:    #c the prescribed index around each cycle. They should add up to N*Euler_characteristic of the mesh.
  //  cycleCurvature: #c the original curvature for each basis cycle.
  //  solver: The Simplicial LDLT solver used to solver the problem. It will only prefactor the matrix once upon the first call to the function; the state of  the solver solely depends on the basisCycles, therefore it only needs to be reset if the basisCycles matrix changed.
  //  N: the degree of the field.
  // Output:
  //  rotationAngles: #iE rotation angles (difference from parallel transport) per inner dual edge
  //  linfError: l_infinity error of the computation. If this is not approximately 0, the prescribed indices are likely inconsistent (don't add up to the correct sum).
  IGL_INLINE void index_prescription(directional::CartesianField& field,
                                     const Eigen::VectorXi& cycleIndices,
                                     Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >& ldltSolver,
                                     const int N,
                                     const double globalRotation,
                                     Eigen::VectorXd& rotationAngles,
                                     double &linfError)
  {
    using namespace Eigen;
    using namespace std;
    
    VectorXd cycleNewCurvature = cycleIndices.cast<double>()*(2.0*igl::PI/(double)N);
    
    //Initialize solver if never before
    if (!ldltSolver.rows())
    {
      SparseMatrix<double> AAt = field.dualCycles*field.dualCycles.transpose();
      ldltSolver.compute(AAt);
    }
    
    VectorXd innerRotationAngles = field.dualCycles.transpose()*ldltSolver.solve(-field.cycleCurvatures + cycleNewCurvature);
    rotationAngles.conservativeResize(field.adjSpaces.rows());
    rotationAngles.setZero();
    for (int i=0;i<field.innerAdjacencies.rows();i++)
      rotationAngles(field.innerAdjacencies(i))=innerRotationAngles(i);
    
    linfError = (field.dualCycles*innerRotationAngles - (-field.cycleCurvatures + cycleNewCurvature)).template lpNorm<Infinity>();
    
    Eigen::MatrixXd representative;
    directional::rotation_to_raw(rotationAngles,N,globalRotation, field);
  }
  
  //Minimal version: no provided solver
  IGL_INLINE void index_prescription(directional::CartesianField& field,
                                     const Eigen::VectorXi& cycleIndices,
                                     const int N,
                                     const double globalRotation,
                                     Eigen::VectorXd& rotationAngles,
                                     double &error)
  {
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > ldltSolver;
    index_prescription(field, ldltSolver, N, rotationAngles, error);
  }
}



#endif


