// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2025 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_SINGLE_TO_N_H
#define DIRECTIONAL_SINGLE_TO_N_H
#include <Eigen/Core>
#include <Eigen/Sparse>
#include <vector>


namespace directional
{

/***Computing a kronecker product of a sparse matrix with Id(N), in order to make a matrix applicable to a raw N-field or N-function
 Input:
 singMat:           the sparse matrix
 N:                     the degree for multiplication
 dRows, dCols:  the unit block in the sparse matrix that gets copied, to allow for vector sizes. For instnace, for (3,3), the matrix creates copies of size (3N, 3N)
 Output:
 The resulting sparse matrix.
 ***/
template <typename Scalar>
Eigen::SparseMatrix<Scalar> single_to_N_matrix(const Eigen::SparseMatrix<Scalar>& singMat,
                                               const int N,
                                               const int dRows,
                                               const int dCols){
    Eigen::SparseMatrix<Scalar> NMat(singMat.rows()*N, singMat.cols()*N);
    std::vector<Eigen::Triplet<Scalar>> NMatTris;
    for (int k=0; k<singMat.outerSize(); ++k)
        for (typename Eigen::SparseMatrix<Scalar>::InnerIterator it(singMat,k); it; ++it)
        {
            int rowPack = int(it.row() / dRows);
            int colPack = int(it.col() / dCols);
            int rowPos = it.row() % dRows;
            int colPos = it.col() % dCols;
            for (int j=0;j<N;j++)
                NMatTris.push_back(Eigen::Triplet<Scalar>(N*dRows*rowPack+j*dRows+rowPos, N*dCols*colPack+j*dCols+colPos, it.value()));
        }
    
    NMat.setFromTriplets(NMatTris.begin(), NMatTris.end());
    return NMat;
    
}



}


#endif
