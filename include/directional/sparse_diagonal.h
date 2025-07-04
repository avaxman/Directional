// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2025 Amir Vaxman <avaxman@gmail.com>
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

//Computing a diagonal sparse matrix from a given vector of values.
template <typename Scalar>
Eigen::SparseMatrix<Scalar> sparse_diagonal(const Eigen::Vector<Scalar, Eigen::Dynamic>& diagValues){
    
    Eigen::SparseMatrix<Scalar> diagMatrix(diagValues.size(), diagValues.size());
    std::vector<Eigen::Triplet<Scalar>> diagMatTris;
    for (int i=0;i<diagValues.size();i++)
        diagMatTris.push_back(Eigen::Triplet<Scalar>(i,i,diagValues(i)));
    
    diagMatrix.setFromTriplets(diagMatTris.begin(), diagMatTris.end());
    return diagMatrix;
    
}

//A version that returns the matrix as a parameter
template <typename Scalar>
inline void sparse_diagonal(const std::vector<Eigen::SparseMatrix<Scalar>>& diagValues,
                            Eigen::SparseMatrix<Scalar>& diagMatrix){
    
    int numRows=0, numCols=0;
    Eigen::MatrixXi offsets(diagValues.size(),2);
    for (int i=0;i<diagValues.size();i++) {
        offsets.row(i)<<numRows, numCols;
        numRows += diagValues[i].rows();
        numCols += diagValues[i].cols();
    }
    diagMatrix.resize(numRows, numCols);
    std::vector<Eigen::Triplet<Scalar>> diagMatTris;
    for (int i=0;i<diagValues.size();i++) {
        for (int k = 0; k < diagValues[i].outerSize(); ++k) {
            for (Eigen::SparseMatrix<double>::InnerIterator it(diagValues[i], k); it; ++it) {
                diagMatTris.push_back(Eigen::Triplet<Scalar>(offsets(i,0)+it.row(), offsets(i,1)+it.col(), it.value()));
            }
        }
    }
    
    diagMatrix.setFromTriplets(diagMatTris.begin(), diagMatTris.end());
    //return diagMatrix;
    
}



}


#endif
