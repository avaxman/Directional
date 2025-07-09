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
#include <directional/definitions.h>
#include <directional/CartesianField.h>
#include <directional/rotation_to_raw.h>


namespace directional
{
// Computes the rotation angles that are required to reproduce a prescribed set of indices on the dual cycles of the mesh.
// In case the sum of curvature is not consistent with the topology, the system is solved in least squares and unexpected singularities may appear elsewhere. linfError will mostl like be far from zero.
// Input:
//  cycleIndices:   a prescribed index per cycle (either local, generator, or boundary). This must it the cycles in the type of the field.
//  N:              degree of the field
//  globalRotation: the orientation of the directional in the first tangent space (mostly arbitrary)
//  ldltSolver:     Since index prescription can benefit from prefactoring, this is an option to give the already-factored solver.
//  field:          Cartesian field object.
// Output:
//  rotationAngles: #adjSpaces rotation angles (difference from parallel transport) per inner space adjacency relation
//  linfError:      l_infinity error of the computation. If this is not approximately 0, the prescribed indices are likely inconsistent (don't add up to the correct sum).
inline void index_prescription(const Eigen::VectorXi& cycleIndices,
                               const int N,
                               const double globalRotation,
                               Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> >& ldltSolver,
                               directional::CartesianField& field,
                               Eigen::VectorXd& rotationAngles,
                               double &linfError)
{
    using namespace Eigen;
    using namespace std;
    
    VectorXd cycleNewCurvature = cycleIndices.cast<double>()*(2.0*directional::PI/(double)N);
    
    //Initialize solver if never before
    if (!ldltSolver.rows())
    {
        SparseMatrix<double> AAt = field.tb->cycles*field.tb->cycles.transpose();
        ldltSolver.compute(AAt);
    }
    
    VectorXd innerRotationAngles = field.tb->cycles.transpose()*ldltSolver.solve(-field.tb->cycleCurvatures + cycleNewCurvature);
    rotationAngles.conservativeResize(field.tb->adjSpaces.rows());
    rotationAngles.setZero();
    for (int i=0;i<field.tb->innerAdjacencies.rows();i++)
        rotationAngles(field.tb->innerAdjacencies(i))=innerRotationAngles(i);
    
    linfError = (field.tb->cycles*innerRotationAngles - (-field.tb->cycleCurvatures + cycleNewCurvature)).template lpNorm<Infinity>();
    
    directional::rotation_to_raw(*(field.tb), rotationAngles,N,globalRotation,field);
}

//Minimal version: without a provided solver
inline void index_prescription(const Eigen::VectorXi& cycleIndices,
                               const int N,
                               const double globalRotation,
                               directional::CartesianField& field,
                               Eigen::VectorXd& rotationAngles,
                               double &error)
{
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > ldltSolver;
    index_prescription(cycleIndices, N, globalRotation,ldltSolver,  field, rotationAngles, error);
}
}



#endif


