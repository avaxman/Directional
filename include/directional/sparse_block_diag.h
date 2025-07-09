// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2025 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_BLOCKDIAG_H
#define DIRECTIONAL_BLOCKDIAG_H

#include <Eigen/Eigen>
#include <numeric>


namespace directional{

/***Computing a block diagonal matrix
 Input:
 matrices:      the blocks to be put on the diagonal of the matrix
 Output:
 C:                The resulting matrix
 ***/

template <typename Scalar>
void sparse_block_diag(
                       const std::vector<Eigen::SparseMatrix<Scalar>*>& matrices,
                       Eigen::SparseMatrix<Scalar> & C)
{
    using namespace Eigen;
    using SMat = Eigen::SparseMatrix<Scalar>;
    
    const size_t newRows = std::accumulate(matrices.begin(), matrices.end(), 0, [](const size_t& s1, const SMat* s2)-> size_t
                                           {
        return s1 + s2->rows();
    });
    const size_t newCols = std::accumulate(matrices.begin(), matrices.end(), 0, [](const size_t& s1, const SMat* s2)-> size_t
                                           {
        return s1 + s2->cols();
    });
    C = SMat(newRows, newCols);
    
    VectorXi colCount;
    colCount.resize(newCols);
    int curr = 0;
    for(int i = 0; i < matrices.size(); i++)
    {
        for (int k = 0; k < matrices[i]->outerSize(); ++k, ++curr)
        {
            int cnt = 0;
            for (typename SMat::InnerIterator it(*matrices[i], k); it; ++it)
            {
                ++cnt;
            }
            colCount(curr) = cnt;
        }
    }
    C.reserve(colCount);
    
    int rOff = 0, cOff = 0;
    for (int i = 0; i < matrices.size(); i++)
    {
        for (int k = 0; k < matrices[i]->outerSize(); ++k)
        {
            for (typename SMat::InnerIterator it(*matrices[i], k); it; ++it)
            {
                C.insert(rOff + it.row(), cOff + k) = it.value();
            }
        }
        rOff += matrices[i]->rows();
        cOff += matrices[i]->cols();
    }
    C.makeCompressed();
}
}

#endif
