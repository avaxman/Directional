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


//Data structure for the extended PolyVector algorithm
struct PolyVectorData{
public:
    
    //User parameters
    Eigen::VectorXi constSpaces;    // List of tangent spaces where there are (partial) constraints. The faces can repeat to constrain more vectors
    Eigen::MatrixXd constVectors;   // Corresponding to constSpaces.
    
    int N;                          // Degree of field
    const TangentBundle* tb;        //The tangent bundle on which the field is defined
    bool verbose;                   //whether to output anything
    bool signSymmetry;              // Whether field enforces a ssign symmetry (only when N is even, otherwise by default set to false)
    bool perfectRoSy;               // Whether the field must be perfect rotationally-symmetric (but not unit).
    double wSmooth;                 // Weight of smoothness
    double wRoSy;                   // Weight of rotational-symmetry. "-1" means a perfect RoSy field (power field)
    Eigen::VectorXd wAlignment;     // Weight of alignment per each of the constfaces. "-1" means a fixed vector
    Eigen::VectorXd confidence;     // The confidence in each output of the iteration functions in the previous step (how much we want to keep it). Should be between [0,1], default is "1" (so it's transparent)
    double initImplicitFactor;      // Implicit smoothing factor
    double currImplicitFactor;      // The current implicit coeff in the given iteration
    double implicitScheduler;       // How much to attenuate the implicit factor with in each iteration
    int iterationMode;              // Making it possible to iterate energy reduction -> some custom projection function
    int currIteration;              // The current iteration
    bool implicitFirst;             // Whether to do the implicit step first or last
 
    Eigen::SparseMatrix<std::complex<double>> smoothMat;    //Smoothness energy
    Eigen::SparseMatrix<std::complex<double>> roSyMat;      //Rotational-symmetry energy
    Eigen::SparseMatrix<std::complex<double>> alignMat;     //(soft) alignment energy.
    Eigen::SparseMatrix<std::complex<double>> reducMat;     //reducing the fixed dofs (for instance with sign symmetry or fixed partial constraints)
    Eigen::VectorXcd reducRhs;                              //The uncompressed PV coeffs are reducMat*true_dofs+reducRhs
    Eigen::VectorXcd alignRhs;                              //encoding the soft constraints
    
    //Mass and stiffness matrices
    Eigen::SparseMatrix<std::complex<double>> WSmooth, WAlign, WRoSy, M;
    double totalRoSyWeight, totalConstrainedWeight, totalSmoothWeight;    //for co-scaling energies
    
    //state-machine solvers and vectors
    Eigen::SparseMatrix<std::complex<double>> totalLhs;             //The total left-hand-side of the system
    Eigen::VectorXcd totalRhs;                                      //The total right-hand side
    Eigen::VectorXcd reducedDofs;                                   //The net solution in the minimal independent degrees of freedom
    Eigen::SparseMatrix<std::complex<double>> implicitLhs;          //The implicit system left-hand-side
    Eigen::VectorXcd implicitRhs;                                   //The implicit right-hand side
    Eigen::SparseMatrix<std::complex<double>> confidenceMat;        //The confidence matrix for implicit iterations
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> reducProjSolver;           //The solver for the reduced dofs from the full dofs (in Least squares)
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>> implicitSolver;            //The solver for the implicit system
    
    PolyVectorData():signSymmetry(true),  tb(NULL), verbose(false), wSmooth(1.0), wRoSy(0.0), iterationMode(false), currIteration(0), currImplicitFactor(0.0), initImplicitFactor(0.5), implicitScheduler(0.8), implicitFirst(true) {wAlignment.resize(0); constSpaces.resize(0); constVectors.resize(0,3);}
    ~PolyVectorData(){}
};

}

#endif
