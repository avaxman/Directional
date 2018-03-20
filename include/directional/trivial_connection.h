// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef TRIVIAL_CONNECTIONS_H
#define TRIVIAL_CONNECTIONS_H
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <igl/parallel_transport_angles.h>
#include <igl/edge_topology.h>
#include <igl/boundary_loop.h>
#include <directional/cycle_curvature.h>
#include <Eigen/Core>
#include <vector>
#include <cmath>


namespace directional
{    
  // Computes the rotation angles to form a trivial connection according to given cone curvatures (or singularity indices) around basis cycles.
  // In case the sum of curvature is not consistent with the topology, the system is solved in least squares and unexpected singularities may appear elsewhere.
  // The output is the modification to the parallel transport.
  // Inputs:
  //  V: #V X 3 vertex coordinates
  //  F: #F by 3 face vertex indices
  //  basisCycles: #basisCycles X #E the oriented (according to EF) basis cycles around which the curvatures are measured
  //               the basis cycles must be arranged so that the first |V| are the vertex cycles (for instance, the result of igl::basis_cycles())
  //  indices: #basisCycles the index around each cycle. They should add up to N*Euler_characteristic of the mesh.
  //  cycleCurvature: the angle defect for each basis cycle.
  //  solver: The Simplicial LDLT solver used to calculate the trivial connections. If initialized the solve step will be skipped when calculating the field.
  //			The state of  the solver solely depends on the basisCycles, therefore it only needs to be reset if the basisCycles matrix changed.
  //			If the solver is not yet set the solver will be called to prepare the basisCycles.
  //  N: the degree of the field. The curvature of a cycle is measured by (singIndex/N)*(2*pi) (can be negative)
  // Outputs:
  //  rotationAngles: the difference between the parallel transport and the modified one.
  //  error: gives the total error of the field. If this is not approximately 0 your singularities probably don't add up properly.
  IGL_INLINE void trivial_connection(const Eigen::MatrixXd& V,
                                     const Eigen::MatrixXi& F,
                                     const Eigen::SparseMatrix<double>& basisCycles,
                                     const Eigen::VectorXd& cycleCurvature,
                                     const Eigen::VectorXi& cycleIndices,
                                     Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >& ldltSolver,
                                     const int N,
                                     Eigen::VectorXd& rotationAngles,
                                     double &linfError)
  {
    using namespace Eigen;
    using namespace std;
    
    VectorXd cycleNewCurvature = cycleIndices.cast<double>()*(2.0*igl::PI/(double)N);
    
    //Initialize solver if never before
    if (!ldltSolver.rows())
    {
      SparseMatrix<double> AAt = basisCycles*basisCycles.transpose();
      ldltSolver.compute(AAt);
    }
    
    rotationAngles = basisCycles.transpose()*ldltSolver.solve(-cycleCurvature + cycleNewCurvature);
    
    linfError = (basisCycles*rotationAngles - (-cycleCurvature + cycleNewCurvature)).lpNorm<Infinity>();
  }
  
  //Minimal version: no solver or pre-computed cycle curvature
  IGL_INLINE void trivial_connection(const Eigen::MatrixXd& V,
                                     const Eigen::MatrixXi& F,
                                     const Eigen::SparseMatrix<double>& basisCycles,
                                     const Eigen::VectorXi& cycleIndices,
                                     const int N,
                                     Eigen::VectorXd& rotationAngles,
                                     double &error)
  {
    Eigen::MatrixXi EV, x, EF;
    igl::edge_topology(V, F, EV, x, EF);
    Eigen::MatrixXd B1, B2, B3;
    igl::local_basis(V, F, B1, B2, B3);
    Eigen::VectorXd cycleCurvature;
    cycle_curvature(V, F, EV, EF, B1, B2, basisCycles, cycleCurvature);
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > ldltSolver;
    trivial_connection(V, F, basisCycles,cycleCurvature, cycleIndices, ldltSolver, N, rotationAngles, error);
  }
}




#endif


