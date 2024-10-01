// This file is part of SaddlePoint, a simple library for Eigen-based sparse nonlinear optimization
//
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_SPARSE_DIAGONAL_H
#define DIRECTIONAL_SPARSE_DIAGONAL_H
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>


namespace directional
{

    template <typename Scalar>
    Eigen::SparseMatrix<Scalar> sparse_diagonal(const Eigen::Vector<Scalar, Eigen::Dynamic>& diagValues){

        Eigen::SparseMatrix<Scalar> diagMatrix(diagValues.size(), diagValues.size());
        std::vector<Eigen::Triplet<Scalar>> diagMatTris;
        for (int i=0;i<diagValues.size();i++)
            diagMatTris.push_back(Eigen::Triplet<Scalar>(i,i,diagValues(i)));

        diagMatrix.setFromTriplets(diagMatTris.begin(), diagMatTris.end());
        return diagMatrix;

    }



}


#endif
