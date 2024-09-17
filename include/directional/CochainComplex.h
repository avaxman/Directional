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

        std::vector<Eigen::Vector<NumberType, Eigen::Dynamic>> cochains;
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

        //Solving the Poisson equation fi-1 = argmin |di*fi-1-1 - fi|^2 with the metric of Mi, and where i = cochainNum
        // (so one needs to feed sequence[cochainNum]=fi)
        void project_exact(const int cochainNum,
                          Eigen::Vector<NumberType, Eigen::Dynamic>& fim1,
                          Eigen::Vector<NumberType, Eigen::Dynamic>& fi){
            if (cochainNum==0)
                return Eigen::Vector<NumberType, Eigen::Dynamic>();  //there is no prev

            assert("cochain is out of bounds!" && cochainNum > 0 && cochainNum < cochains.size()-1);
            Eigen::SparseMatrix<NumberType> Mi = metrics[cochainNum];
            Eigen::SparseMatrix<NumberType> di = differentials[cochainNum-1];
            assert("the cochain did not implement the metric function " && Mi.nonZeros()!=0);

            Eigen::SparseMatrix<NumberType> A = di.adjoint() * Mi * di;
            Eigen::Vector<NumberType, Eigen::Dynamic> b = di.adjoint() * Mi * cochains[cochainNum];
            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver;
            solver.compute(A);
            assert("project_prev(): Decomposition failed!" && solver.info() == Eigen::Success);
            fim1 = solver.solve(b);
            assert("project_prev(): Solver failed!" & solver.info() == Eigen::Success);
            fi = di*fim1;
        }

        //Solving the equation fi = argmin |fi|^2 s.t. di*fi = fi+1 with the metric of Mi
        //and where i = cochainNum. The input is then sequence[cochainNum].
        // Note this is conceptually different than project_next!
        // the gip1 returned is the lagrange multipliers, also the dual coexact generator of fi (so that Mi*fi = di^T*g)
        void project_coexact(const int cochainNum,
                             Eigen::Vector<NumberType, Eigen::Dynamic>& gip1,
                             Eigen::Vector<NumberType, Eigen::Dynamic>& fi){
            if (cochainNum==cochains.size())
                return Eigen::Vector<NumberType, Eigen::Dynamic>();  //there is no next

            assert("cochain is out of bounds!" && cochainNum > 0 && cochainNum < sequence.size()-1);
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
            b.tail(di.rows())=cochains[cochainNum+1];

            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.compute(A);
            assert("project_next(): Decomposition failed!" && solver.info() == Eigen::Success);
            Eigen::Vector<NumberType, Eigen::Dynamic> x = solver.solve(b);
            assert("project_next(): Solver failed!" && solver.info() == Eigen::Success);
            fi = x.head(Mi.rows());
            gip1 = -x.tail(di.rows());
        }

        //computes the hodge decomposition of a cochain i(=cochainNum)
        void hodge_decomposition(int cochainNum,
                                 Eigen::Vector<NumberType, Eigen::Dynamic>& exactPrev,
                                 Eigen::Vector<NumberType, Eigen::Dynamic>& exactCurr,
                                 Eigen::Vector<NumberType, Eigen::Dynamic>& coexactNext,
                                 Eigen::Vector<NumberType, Eigen::Dynamic>& coexactCurr,
                                 Eigen::Vector<NumberType, Eigen::Dynamic>& harmCurr){

            project_exact(cochainNum, exactPrev, exactCurr);
            project_coexact(cochainNum, coexactNext, coexactCurr);
            //TODO: perhaps extracting "to_vector()" is inefficient
            harmCurr = cochains[cochainNum] - exactCurr - coexactCurr;
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

        void codifferential(const int cochainNum, Eigen::Vector<NumberType, Eigen::Dynamic>& result){
            //TODO
        }

        void codifferential_matrix(const int cochainNum){
            assert("Inverse metrics have not been defined!" && invMetrics.size()==0);
            //TODO
        }

       CochainComplex dual_complex(){
            //TODO: return dual complex
        }

    };

}


#endif
