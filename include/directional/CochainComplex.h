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

    //TODO: should be solving also a minimal problem since we need the solution to be orthogonal to the null space in the prevM metric...
    //Solving the Poisson equation prevCochain = argmin |di*prevCochain - cochain|^2 with the metric of Mi, and where i = cochainNum, and also outputting exactCochain = di*prevCochain
    template<typename NumberType>
    void project_exact(const Eigen::SparseMatrix<NumberType> d,
                       const Eigen::SparseMatrix<NumberType> M,
                       const Eigen::Vector<NumberType, Eigen::Dynamic>& cochain,
                       Eigen::Vector<NumberType, Eigen::Dynamic>& prevCochain,
                       Eigen::Vector<NumberType, Eigen::Dynamic>& exactCochain){

        Eigen::SparseMatrix<NumberType> A = d.adjoint() * M * d;
        Eigen::Vector<NumberType, Eigen::Dynamic> b = d.adjoint() * M * cochain;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        assert("project_prev(): Decomposition failed!" && solver.info() == Eigen::Success);
        prevCochain = solver.solve(b);
        assert("project_prev(): Solver failed!" && solver.info() == Eigen::Success);
        exactCochain = d*prevCochain;
    }

    //Solving the equation nextCochain = argmin |codiff_i*nextCochain - cochain|^2 with the metric of Mi, and where i = cochain, and also outputting coexactCochain = codiff_i*nextCochain
    //However, since the codifferential is not assumed to be given explicitly, we instead solve for:
    // coexactCochain = argmin |coexactCochain|^2 s.t. M_(i+1)*di*cochain = M_(i+1)*di*coexactCochain with the metric of Mi (i = cochainNum)
    // nextCochain is the (i+1) form such that cocexactCochain = codiff*nextCochain (codifferential)
    template<typename NumberType>
    void project_coexact(const Eigen::SparseMatrix<NumberType> d,
                         const Eigen::SparseMatrix<NumberType> M,
                         const Eigen::SparseMatrix<NumberType> MNext,
                         const Eigen::Vector<NumberType, Eigen::Dynamic>& cochain,
                         Eigen::Vector<NumberType, Eigen::Dynamic>& coexactCochain,
                         Eigen::Vector<NumberType, Eigen::Dynamic>& nextCochain){

        Eigen::Matrix2i orderMat; orderMat<<0,1,2,3;
        std::cout<<"orderMat: "<<orderMat<<std::endl;
        std::vector<Eigen::SparseMatrix<NumberType>> matVec;
        matVec.push_back(M);
        matVec.push_back(d.adjoint()*MNext);
        matVec.push_back(MNext*d);
        matVec.push_back(Eigen::SparseMatrix<NumberType>(MNext.rows(), d.rows()));
        Eigen::SparseMatrix<NumberType> A;
        directional::sparse_block(orderMat, matVec, A);
        Eigen::Vector<NumberType, Eigen::Dynamic> b(M.rows()+d.rows());
        b.head(M.rows()).setZero();
        b.tail(d.rows())=MNext*d*cochain;

        Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
        solver.compute(A);
        assert("project_next(): Decomposition failed!" && solver.info() == Eigen::Success);
        Eigen::Vector<NumberType, Eigen::Dynamic> x = solver.solve(b);
        assert("project_next(): Solver failed!" && solver.info() == Eigen::Success);
        coexactCochain = x.head(M.rows());
        nextCochain = -x.tail(d.rows());
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
                             const Eigen::SparseMatrix<NumberType> MNext,
                             const Eigen::Vector<NumberType, Eigen::Dynamic>& cochain,
                             Eigen::Vector<NumberType, Eigen::Dynamic>& prevCochain,
                             Eigen::Vector<NumberType, Eigen::Dynamic>& exactCochain,
                             Eigen::Vector<NumberType, Eigen::Dynamic>& coexactCochain,
                             Eigen::Vector<NumberType, Eigen::Dynamic>& nextCochain,
                             Eigen::Vector<NumberType, Eigen::Dynamic>& harmCochain){

        std::cout<<"before project exact"<<std::endl;
        project_exact(d, M, cochain, prevCochain, exactCochain);
        std::cout<<"before project coexact"<<std::endl;
        project_coexact(dNext, M, MNext, cochain, coexactCochain, nextCochain);
        harmCochain = cochain - exactCochain - coexactCochain;
    }

    //Computes a basis (orthogonal under the metric) of harmonic fields for a k-form where
    //dNext is the k->k+1 differential, and d is the k-1->k differential.
    //This constructs the harmonic operator [dNext; d.adjoint()*M] and takes its kernel.
    //if orthogonalize = true, the kernel is made so it's M-orthogonal.
    //To get the proper boundary conditions, the operators need to be appropriate
    template<typename NumberType>
    void cohomology_basis(const Eigen::SparseMatrix<NumberType> d,
                         const Eigen::SparseMatrix<NumberType> dNext,
                         const Eigen::SparseMatrix<NumberType> M,
                         Eigen::Matrix<NumberType, Eigen::Dynamic, Eigen::Dynamic>& harmFields,
                         int& bettiNumber,
                         const bool orthogonalize = true){

        std::vector<Eigen::SparseMatrix<NumberType>> matVec;
        Eigen::MatrixXi matOrder(2,1);
        matOrder<<0,1;
        matVec.push_back(dNext);
        matVec.push_back(d.adjoint()*M);
        Eigen::SparseMatrix<NumberType> H;
        //std::cout<<"Before sparse block"<<std::endl;
        directional::sparse_block(matOrder, matVec, H);
        //std::cout<<"After sparse block"<<std::endl;
        Eigen::SparseQR<Eigen::SparseMatrix<NumberType>, Eigen::COLAMDOrdering<int>> qr;
        qr.compute(H.adjoint());
        assert("harmonic_field(): Decomposition failed!" && qr.info() == Eigen::Success);
        bettiNumber = H.cols() - qr.rank();
        //std::cout<<"BettiNumber: "<<bettiNumber<<std::endl;
        //Eigen::VectorXi PIndices = qr.colsPermutation().indices();
        Eigen::SparseMatrix<NumberType> I(qr.matrixQ().cols(), bettiNumber);
        std::vector<Eigen::Triplet<NumberType>> ITris;
        for (int i=0;i<bettiNumber;i++) {
            ITris.push_back(Eigen::Triplet<NumberType>(qr.matrixQ().cols()-i-1, bettiNumber - i - 1, 1));
            //std::cout<<PIndices(PIndices.size() - i - 1)<<","<<bettiNumber - i - 1<<std::endl;
        }
        I.setFromTriplets(ITris.begin(), ITris.end());
        //std::cout<<"Eigen::MatrixXd(I): "<<Eigen::MatrixXd(I)<<std::endl;
        harmFields = Eigen::MatrixXd(qr.matrixQ() * Eigen::MatrixXd(I)); //TODO: inefficient
        //std::cout<<"H*Q: "<<(H*Eigen::MatrixXd(qr.matrixQ()))<<std::endl;
        //std::cout<<"harmonicity: "<<(H*harmFields).cwiseAbs().maxCoeff()<<std::endl;

    }

        /*void codifferential(const int cochainNum, const Eigen::Vector<NumberType, Eigen::Dynamic>& cochain, Eigen::Vector<NumberType, Eigen::Dynamic>& result){
            assert("The first cochain in the sequence doesn't have a codifferential" && cochainNum>0);
            if (invMetrics[cochainNum].size()!=0) //inverse matrices are given
                result =  invMetrics[cochainNum]*differentials[cochainNum-1].adjoint()*metrics[cochainNum]*cochain;
            else{
                Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
                solver.compute(metrics[cochainNum]);
                result = solver.solve(differentials[cochainNum-1].adjoint()*metrics[cochainNum]*cochain);
            }

        }*/

        /*void dual_differential(const int cochainNum, const Eigen::Vector<NumberType, Eigen::Dynamic>& dualCochain, Eigen::Vector<NumberType, Eigen::Dynamic>& result){
            assert("The first cochain in the sequence doesn't have a codifferential" && cochainNum>0);
            if (invMetrics.size()!=0) //inverse matrices are given
                result =  invMetrics[cochainNum]*differentials[cochainNum-1].adjoint()*dualCochain;
            else{
                Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
                solver.compute(metrics[cochainNum]);
                result = solver.solve(differentials[cochainNum-1].adjoint()*dualCochain);
            }
        }*/

        /*void codifferential_matrix(const int cochainNum){
            assert("codifferential_matrix(): Inverse metrics have not been defined!" && invMetrics[cochainNum].size()==0);
            return invMetrics[cochainNum]*differentials[cochainNum-1].adjoint()*metrics[cochainNum];
        }*/

   // };

}


#endif
