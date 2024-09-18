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

    //an interface to a single element in the class
    /*template<typename NumberType>
    class Cochain{
    public:

        //An unrolled version of the space that can works with a metric matrix (for instance)
        virtual Eigen::Vector<NumberType, Eigen::Dynamic> to_vector() = 0;

        //Updating back to the object
        virtual void from_vector(const Eigen::Vector<NumberType, Eigen::Dynamic>& vec) = 0;

        //The differential operator that returns an unrolled different cochain (of a different type) down the sequence
        virtual Eigen::SparseMatrix<NumberType> differential() = 0;

        //Inner product matrix so that inner product on the cochain is to_vector().adjoint() * metric() * t_vector()
        //The metric is optional, but necessary for several operators in a cochain sequence
        virtual Eigen::SparseMatrix<NumberType> metric(){return Eigen::SparseMatrix<NumberType>();}

        //This is another optional function that is only useful for explicit codifferential or dual cochains
        virtual Eigen::SparseMatrix<NumberType> inverse_metric(){return Eigen::SparseMatrix<NumberType>();}
    };

    template<typename NumberType, typename PrimalCochain>
    class DualCochain: public Cochain<NumberType>{
        static_assert(std::is_base_of<Cochain<NumberType>, PrimalCochain>::value, "The primal cochain must be of type class Cochain");

        Eigen::SparseMatrix<NumberType> differential(){
            return
        }

    };*/

    template<typename NumberType>
    class CochainComplex {
    public:

        std::vector<Eigen::SparseMatrix<NumberType>> differentials;
        std::vector<Eigen::SparseMatrix<NumberType>> metrics;

        //optional: used if exist, otherwise an equation is solved
        std::vector<Eigen::SparseMatrix<NumberType>> invMetrics;

        CochainComplex(){}
        ~CochainComplex(){}

        //checking that the cochain is propoer, that is d^2 = 0 and d^{N-1} = 0
        bool check_cochain_validity(double tolerance){
            for (int i=0;i<differentials.size()-1;i++){
                Eigen::SparseMatrix<NumberType> d2 = differentials[i+1]*differentials[i];
                double maxAbsValue=-32767000.0;
                for (int k = 0; k < d2.outerSize(); ++k) {
                    for (typename Eigen::SparseMatrix<NumberType>::InnerIterator it(d2, k); it; ++it) {
                        double absValue = std::abs(it.value());
                        if (absValue > maxAbsValue) {
                            maxAbsValue = absValue;
                        }
                    }
                }
                if (maxAbsValue > tolerance)
                    return false;
            }
            return true;
        }

        //Solving the Poisson equation prevCochain = argmin |di*prevCochain - cochain|^2 with the metric of Mi, and where i = cochainNum, and also outputting exactCochain = di*prevCochain
        void project_exact(const int cochainNum,
                          const Eigen::Vector<NumberType, Eigen::Dynamic>& cochain,
                          Eigen::Vector<NumberType, Eigen::Dynamic>& prevCochain,
                          Eigen::Vector<NumberType, Eigen::Dynamic>& exactCochain){
            if (cochainNum==0)
                return Eigen::Vector<NumberType, Eigen::Dynamic>();  //there is no prev

            assert("cochain is out of bounds!" && cochainNum > 0 && cochainNum < differentials.size()-1);
            Eigen::SparseMatrix<NumberType> Mi = metrics[cochainNum];
            Eigen::SparseMatrix<NumberType> di = differentials[cochainNum-1];
            assert("the cochain did not implement the metric function " && Mi.nonZeros()!=0);

            Eigen::SparseMatrix<NumberType> A = di.adjoint() * Mi * di;
            Eigen::Vector<NumberType, Eigen::Dynamic> b = di.adjoint() * Mi * cochains[cochainNum];
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
            solver.compute(A);
            assert("project_prev(): Decomposition failed!" && solver.info() == Eigen::Success);
            prevCochain = solver.solve(b);
            assert("project_prev(): Solver failed!" & solver.info() == Eigen::Success);
            exactCochain = di*fim1;
        }

        //Solving the equation coexactCochain = argmin |coexactCochain|^2 s.t. di*cochain = di*coexactCochain with the metric of Mi
        // the nextDualCochain returned is the lagrange multipliers, also the dual coexact generator of fi (so that Mi*coexactCochain = di^T*nextDualCochain)
        void project_coexact(const int cochainNum,
                             const Eigen::Vector<NumberType, Eigen::Dynamic>& cochain,
                             Eigen::Vector<NumberType, Eigen::Dynamic>& coexactCochain,
                             Eigen::Vector<NumberType, Eigen::Dynamic>& nextDualCochain){
            if (cochainNum==cochains.size())
                return Eigen::Vector<NumberType, Eigen::Dynamic>();  //there is no next

            assert("cochain is out of bounds!" && cochainNum > 0 && cochainNum < differentials.size()-1);
            Eigen::SparseMatrix<NumberType> Mi = metrics[cochainNum];
            Eigen::SparseMatrix<NumberType> di = differentials[cochainNum];
            assert("the cochain did not implement the metric function " && Mi.nonZeros()!=0);
            Eigen::Matrix2i orderMat; orderMat<<1,2,3,4;
            std::vector<Eigen::SparseMatrix<NumberType>> matVec;
            matVec.push_back(Mi);
            matVec.push_back(di.adjoint());
            matVec.push_back(di);
            matVec.push_back(Eigen::SparseMatrix<NumberType>(Mi.rows(), di.adjoint()));
            Eigen::SparseMatrix<NumberType> A = sparse_block(orderMat, matVec);
            Eigen::Vector<NumberType, Eigen::Dynamic> b(Mi.rows()+di.rows());
            b.head(Mi.rows()).setZero();
            b.tail(di.rows())=differentials[cochainNum]*cochain;

            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(A);
            assert("project_next(): Decomposition failed!" && solver.info() == Eigen::Success);
            Eigen::Vector<NumberType, Eigen::Dynamic> x = solver.solve(b);
            assert("project_next(): Solver failed!" && solver.info() == Eigen::Success);
            coexactCochain = x.head(Mi.rows());
            nextDualCochain = -x.tail(di.rows());
        }

        //computes the hodge decomposition of a cochain i(=cochainNum)
        void hodge_decomposition(const int cochainNum,
                                 const Eigen::Vector<NumberType, Eigen::Dynamic>& cochain,
                                 Eigen::Vector<NumberType, Eigen::Dynamic>& prevCochain,
                                 Eigen::Vector<NumberType, Eigen::Dynamic>& exactCochain
                                 Eigen::Vector<NumberType, Eigen::Dynamic>& coexactCochain,
                                 Eigen::Vector<NumberType, Eigen::Dynamic>& nextDualCochainת
                                 Eigen::Vector<NumberType, Eigen::Dynamic>& harmCurr){

            project_exact(cochainNum, cochain, prevCochain, exactCochain);
            project_coexact(cochainNum, cochain, coexactCochain, nextDualCochain);
            harmCurr = cochain - exactCurr - coexactCurr;
        }

        //Computes a basis (orthogonal under the metric) of harmonic fields for a given level i(=cochainNum)
        //If the cochain has a metric, the fields are diagonalized through it.
        void harmonic_fields(const int cochainNum,
                             Eigen::Matrix<NumberType, Eigen::Dynamic, Eigen::Dynamic>& harmFields,
                             int& bettiNumber,
                             const bool orthogonalize = true){

            std::vector<Eigen::SparseMatrix<NumberType>> matVec;
            Eigen::MatrixXi matOrder(2,1);
            matOrder<<1,2;
            matVec.push_back(differentials[cochainNum-1]);
            matVec.push_back(differentials[cochainNum].adjoint());
            Eigen::SparseMatrix<NumberType> H = sparse_block(matOrder, matVec);
            Eigen::SparseQR<Eigen::SparseMatrix<NumberType>, Eigen::COLAMDOrdering<int>> qr;
            qr.compute(H);
            assert("harmonic_field(): Decomposition failed!" && qr.info() == Eigen::Success);
            harmFields = qr.matrixQ().rightCols(H.cols() - qr.rank());
            bettiNumber = harmFields.cols();

            //TODO: orthogonalize according to Mi
        }

        void codifferential(const int cochainNum, const Eigen::Vector<NumberType, Eigen::Dynamic>& cochain, Eigen::Vector<NumberType, Eigen::Dynamic>& result){
            assert("The first cochain in the sequence doesn't have a codifferential" && cochainNum>0);
            if (invMetrics.size()!=0) //inverse matrices are given
                result =  invMetrics[cochainNum]*differentials[cochainNum-1].adjoint()*metrics[cochainNum]*cochain;
            else{
                Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
                solver.compute(metrics[cochainNum]);
                result = solver.solve(differentials[cochainNum-1].adjoint()*metrics[cochainNum]*cochain);
            }

        }

        void dual_differential(const int cochainNum, const Eigen::Vector<NumberType, Eigen::Dynamic>& dualCochain, Eigen::Vector<NumberType, Eigen::Dynamic>& result){
            assert("The first cochain in the sequence doesn't have a codifferential" && cochainNum>0);
            if (invMetrics.size()!=0) //inverse matrices are given
                result =  invMetrics[cochainNum]*differentials[cochainNum-1].adjoint()*dualCochain;
            else{
                Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
                solver.compute(metrics[cochainNum]);
                result = solver.solve(differentials[cochainNum-1].adjoint()*dualCochain);
            }
        }

        void codifferential_matrix(const int cochainNum){
            assert("codifferential_matrix(): Inverse metrics have not been defined!" && invMetrics.size()==0);
            return invMetrics[cochainNum]*differentials[cochainNum-1].adjoint()*metrics[cochainNum];
        }

        //creating a cochain complex that is dual to the current one
       CochainComplex dual_complex(){
            assert("dual_complex(): Inverse metrics have not been defined!" && invMetrics.size()==0);
            CochainComplex dualComplex;
            dualComplex.differentials.resize(differentials.size());
            dualComplex.metrics = invMetrics;
            dualComplex.invMetrics = metrics;
            for (int i=0;i<differentials.size();i++)
                dualComplex.differentials[differentials.size()-1-i] = invMetrics[i]*differentials[i].adjoint();

        }
    };

}


#endif
