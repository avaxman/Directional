// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2025 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POLYVECTOR_DATA_H
#define DIRECTIONAL_POLYVECTOR_DATA_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace directional
{


//Data for the precomputation of the PolyVector algorithm
struct PolyVectorData{
public:
    
    //User parameters
    Eigen::VectorXi constSpaces;    // List of tangent spaces where there are (partial) constraints. The faces can repeat to constrain more vectors
    Eigen::MatrixXd constVectors;   // Corresponding to constSpaces.
    
    int N;                          // Degree of field
    TangentBundle* tb;              //The tangent bundle on which the field is defined
    bool verbose;                   //whether to output anything
    bool signSymmetry;              // Whether field enforces a ssign symmetry (only when N is even, otherwise by default set to false)
    bool perfectRoSy;               // Whether the field must be perfect rotationally-symmetric (but not unit).
    double wSmooth;                 // Weight of smoothness
    double wRoSy;                   // Weight of rotational-symmetry. "-1" means a perfect RoSy field (power field)
    Eigen::VectorXd wAlignment;     // Weight of alignment per each of the constfaces. "-1" means a fixed vector
    double initImplicitFactor;          // Implicit smoothing factor
    double currImplicitCoeff;
    double implicitScheduler;        //How much to attenuate implicit factor by
    int numIterations;              //  Iterate energy reduction->(possibly)normalize->(possibly)project curl
    int currIteration;
    
    Eigen::SparseMatrix<std::complex<double>> smoothMat;    //Smoothness energy
    Eigen::SparseMatrix<std::complex<double>> roSyMat;      //Rotational-symmetry energy
    Eigen::SparseMatrix<std::complex<double>> alignMat;     //(soft) alignment energy.
    Eigen::SparseMatrix<std::complex<double>> reducMat;     //reducing the fixed dofs (for instance with sign symmetry or fixed partial constraints)
    Eigen::VectorXcd reducRhs;                              //The uncompressed PV coeffs are reducMat*true_dofs+reducRhs
    Eigen::VectorXcd alignRhs;                              //encoding the soft constraints
    
    //Mass and stiffness matrices
    Eigen::SparseMatrix<std::complex<double>> WSmooth, WAlign, WRoSy, M;
    double totalRoSyWeight, totalConstrainedWeight, totalSmoothWeight;    //for co-scaling energies
    
    PolyVectorData():signSymmetry(true),  tb(NULL), verbose(false), wSmooth(1.0), wRoSy(0.0), numIterations(0), currIteration(0), currImplicitCoeff(0.0), initImplicitFactor(0.5), implicitScheduler(0.8) {wAlignment.resize(0); constSpaces.resize(0); constVectors.resize(0,3);}
    ~PolyVectorData(){}
};

}

#endif
