// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_COCHAIN_COMPLEX_H
#define DIRECTIONAL_COCHAIN_COMPLEX_H

#include <cassert>
#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <directional/sparse_block.h>

/***
 This is the interface class for cochain complexes, which are a sequence of commutative groups related by operators such
 that the image of one is contained in the kernel of the other, and where the remainders are the cohomology of the complex.
 
 The class is quite abstract, and assumes 1) the rings are vector spaces, 2) the operators are sparse and linear. In addition,
 the class contain inner-product matrices (metrics) that implement the norms in each vector space. This allows for a generalized
 Hodge decomposition, and solving Poisson systems.
 ***/

namespace directional{

//GaugeFixing = false means we are solving a generic Poisson equation prevCochain = argmin |d*prevCochain - cochain|^2 with the metric of M, outputting exactCochain = di*prevCochain, which doesn't care about d null space components in prevCochain.
//This has a precondition that d0 is injective, that is that d'*d is invertible. It will likely work softly (for instance for d0), but no guarantee.
//GaugeFixing = true means we want a minimum 2-norm solution (orthogonal to the null space of d), so we solve instead for:
// prevCochain = argmin |prevCochain|^2 s.t. d*prevCochain = cochain, with the metric of M.  The precondition is that d is surjective in the space where cochain is defined, meaning that the constraint is satisfiable.
//The mass matrix if either for the cochain (gaugeFixing true) or the prevChain (gaugeFixing false)
template<typename NumberType>
void project_exact(const Eigen::SparseMatrix<NumberType> d,
                   const Eigen::SparseMatrix<NumberType> M,
                   const Eigen::Vector<NumberType, Eigen::Dynamic>& cochain,
                   Eigen::Vector<NumberType, Eigen::Dynamic>& prevCochain,
                   Eigen::Vector<NumberType, Eigen::Dynamic>& exactCochain,
                   const bool gaugeFixing = false){
    
    if (!gaugeFixing) {
        Eigen::SparseMatrix<NumberType> A = d.adjoint() * M * d;
        Eigen::Vector<NumberType, Eigen::Dynamic> b = d.adjoint() * M * cochain;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<NumberType>> solver;
        solver.compute(A);
        assert(solver.info() == Eigen::Success && "project_exact(): Decomposition failed!");
        prevCochain = solver.solve(b);
        assert(solver.info() == Eigen::Success && "project_exact(): Solver failed!");
        
    } else {
        Eigen::Matrix2i orderMat; orderMat<<0,1,2,3;
        std::vector<Eigen::SparseMatrix<NumberType>> matVec;
        matVec.push_back(M);
        matVec.push_back(d.adjoint());
        matVec.push_back(d);
        matVec.push_back(Eigen::SparseMatrix<NumberType>(d.rows(), d.rows()));
        Eigen::SparseMatrix<NumberType> A;
        directional::sparse_block(orderMat, matVec, A);
        Eigen::Vector<NumberType, Eigen::Dynamic> b(M.rows()+d.rows());
        b.head(M.rows()).setZero();
        b.tail(d.rows())=cochain;
        
        Eigen::SparseLU<Eigen::SparseMatrix<NumberType>> solver;
        solver.compute(A);
        
        if (solver.info() != Eigen::Success) {
            std::cerr << "project_exact(): Decomposition failed!" << std::endl;
            std::abort();
        }
        Eigen::Vector<NumberType, Eigen::Dynamic> x = solver.solve(b);
        if (solver.info() != Eigen::Success) {
            std::cerr << "project_exact(): Solving failed!" << std::endl;
            std::abort();
        }
        prevCochain = x.head(M.rows());
    }
    exactCochain = d * prevCochain;
}


//We solve for the closest coexact field and its generator nextChain
// coexactCochain = argmin |coexactCochain|^2 s.t. di*cochain = di*coexactCochain with the metric of Mi (i = cochainNum)
// nextCochain is the (i+1) form such that cocexactCochain = -M^{-1}*d1^T*MNext*nextCochain (codifferential)
// k is the order of the differential k-form (k = 1 by default, meaning a 1-form)
// Note the sign! the is the canonical definition so that the we have codiff*nextCochain = coexactCochain, where codiff = (-1)^k*invHodge*d^T*MNext
template<typename NumberType>
void project_coexact(const Eigen::SparseMatrix<NumberType> d,
                     const Eigen::SparseMatrix<NumberType> M,
                     const Eigen::SparseMatrix<NumberType> MNext,
                     const Eigen::Vector<NumberType, Eigen::Dynamic>& cochain,
                     Eigen::Vector<NumberType, Eigen::Dynamic>& coexactCochain,
                     Eigen::Vector<NumberType, Eigen::Dynamic>& nextCochain
                     double k = 1){
    
    Eigen::Matrix2i orderMat; orderMat<<0,1,2,3;
    std::vector<Eigen::SparseMatrix<NumberType>> matVec;
    matVec.push_back(M);
    matVec.push_back(d.adjoint()*MNext.adjoint());
    matVec.push_back(Mext*d);
    matVec.push_back(Eigen::SparseMatrix<NumberType>(MNext.rows(), d.rows()));
    Eigen::SparseMatrix<NumberType> A;
    directional::sparse_block(orderMat, matVec, A);
    Eigen::Vector<NumberType, Eigen::Dynamic> b(M.rows()+d.rows());
    b.head(M.rows()).setZero();
    b.tail(d.rows())=d*cochain;
    
    Eigen::SparseLU<Eigen::SparseMatrix<NumberType>> solver;
    solver.compute(A);
    assert("project_next(): Decomposition failed!" && solver.info() == Eigen::Success);
    Eigen::Vector<NumberType, Eigen::Dynamic> x = solver.solve(b);
    assert("project_next(): Solver failed!" && solver.info() == Eigen::Success);
    coexactCochain = x.head(M.rows());
    double sign = (k%2==0 ? 1 : -1);  //it will actually be a -1 for odd k since x is in the left-hand side
    nextCochain = sign*x.tail(d.rows());
}

// Solve A = H^T H (sparse), M mass matrix (dense or sparse), k eigenvectors near zero
// maxIters and tol control convergence
template<typename Scalar>
void computeSmallestEigenvectors(
                                 const Eigen::SparseMatrix<Scalar>& H,
                                 const Eigen::SparseMatrix<Scalar>& M,
                                 int k,
                                 Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic>& eigenvectors,
                                 int maxIters = 1000,
                                 double tol = 1e-8)
{
    using Vec = Eigen::Matrix<Scalar, Eigen::Dynamic, 1>;
    int n = H.cols();
    
    // Compute A = H^T * H
    Eigen::SparseMatrix<Scalar> A = H.transpose() * H;
    
    // Solver for (A + sigma I) x = b
    Scalar sigma = 1e-4; // small positive shift to stabilize inverse
    Eigen::SparseMatrix<Scalar> Ashift = A;
    for (int i = 0; i < Ashift.outerSize(); ++i)
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(Ashift, i); it; ++it)
            if (it.row() == it.col())
                it.valueRef() += sigma;
    
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<Scalar>> solver;
    solver.compute(Ashift);
    if (solver.info() != Eigen::Success) {
        std::cerr << "Decomposition failed\n";
        return;
    }
    
    eigenvectors.resize(n, k);
    std::vector<Vec> vecs(k);
    
    // Initialize random vectors
    for (int i = 0; i < k; ++i) {
        vecs[i] = Vec::Random(n);
        // Optional: M-orthonormalize initial vecs[i] here
    }
    
    for (int iter = 0; iter < maxIters; ++iter) {
        bool converged = true;
        for (int i = 0; i < k; ++i) {
            // Solve (A + sigma I) y = x
            Vec y = solver.solve(vecs[i]);
            if (solver.info() != Eigen::Success) {
                std::cerr << "Solve failed at iteration " << iter << "\n";
                return;
            }
            
            // Orthogonalize y against previous vectors with respect to M
            for (int j = 0; j < i; ++j) {
                Scalar dot = vecs[j].transpose() * M * y;
                y -= dot * vecs[j];
            }
            
            // Normalize y with respect to M-inner product
            Scalar norm = std::sqrt(y.transpose() * M * y);
            if (norm < tol) {
                std::cerr << "Near zero vector detected, reinitializing\n";
                y = Vec::Random(n);
                norm = std::sqrt(y.transpose() * M * y);
            }
            y /= norm;
            
            // Check convergence
            Scalar diff = (y - vecs[i]).norm();
            if (diff > tol)
                converged = false;
            
            vecs[i] = y;
        }
        if (converged) break;
    }
    
    // Copy results to output matrix
    for (int i = 0; i < k; ++i) {
        eigenvectors.col(i) = vecs[i];
    }
}


//computes the hodge decomposition of a cochain, where
//cochain = exactCochain + coexactCochain + harmCochain and where
// exactCochain = d*prevCochain,
// coexactCochain = M^{-1}*dNext*MNext*nextCochain
// and harmCochain is the remainder.
template<typename NumberType>
void hodge_decomposition(const Eigen::SparseMatrix<NumberType> d,
                         const Eigen::SparseMatrix<NumberType> dNext,
                         const Eigen::SparseMatrix<NumberType> M,
                         const Eigen::Vector<NumberType, Eigen::Dynamic>& cochain,
                         Eigen::Vector<NumberType, Eigen::Dynamic>& prevCochain,
                         Eigen::Vector<NumberType, Eigen::Dynamic>& exactCochain,
                         Eigen::Vector<NumberType, Eigen::Dynamic>& coexactCochain,
                         Eigen::Vector<NumberType, Eigen::Dynamic>& nextCochain,
                         Eigen::Vector<NumberType, Eigen::Dynamic>& harmCochain){
    
    project_exact(d, M, cochain, prevCochain, exactCochain);
    project_coexact(dNext, M, cochain, coexactCochain, nextCochain);
    harmCochain = cochain - exactCochain - coexactCochain;
}

//Computes a basis (orthogonal under the metric) of harmonic fields for a k-form where
//dNext is the k->k+1 differential, and d is the k-1->k differential.
//This constructs the harmonic operator [dNext; d.adjoint()*M] and takes its (approximate) kernels as the lowest bettiNumber eigenvalues
template<typename NumberType>
void cohomology_basis(const Eigen::SparseMatrix<NumberType> d,
                      const Eigen::SparseMatrix<NumberType> dNext,
                      const Eigen::SparseMatrix<NumberType> M,
                      const int bettiNumber,
                      Eigen::Matrix<NumberType, Eigen::Dynamic, Eigen::Dynamic>& harmFields) {
    
    std::vector<Eigen::SparseMatrix<NumberType>> matVec;
    Eigen::MatrixXi matOrder(2, 1);
    matOrder << 0, 1;
    matVec.push_back(dNext);
    matVec.push_back(d.adjoint() * M);
    
    Eigen::SparseMatrix<NumberType> H;
    directional::sparse_block(matOrder, matVec, H);
    
    /*std::cout<<"Decomposition "<<std::endl;
     Eigen::SparseQR<Eigen::SparseMatrix<NumberType>, Eigen::COLAMDOrdering<int>> qr;
     qr.compute(H.adjoint());
     assert(qr.info() == Eigen::Success && "harmonic_field(): Decomposition failed!");
     
     bettiNumber = H.cols() - qr.rank();
     //std::cout << "bettiNumber: " << bettiNumber << std::endl;
     
     const int Qcols = qr.matrixQ().cols();
     harmFields.resize(Qcols, bettiNumber);
     
     std::cout<<"Extracting field "<<std::endl;
     for (int i = 0; i < bettiNumber; ++i) {
     // Construct a standard basis vector e_{Qcols - bettiNumber + i}
     Eigen::VectorXd ei = Eigen::VectorXd::Zero(Qcols);
     ei(Qcols - bettiNumber + i) = 1.0;
     
     // Apply Q to ei to get the i-th harmonic field
     Eigen::VectorXd col = qr.matrixQ() * ei;
     //Normalizing col according to the mass matrix
     col.array()/=sqrt((col.adjoint()*M*col).coeff(0,0));
     harmFields.col(i) = col;
     }*/
    
    computeSmallestEigenvectors(H, M, bettiNumber, harmFields);
    
    
}


#endif
