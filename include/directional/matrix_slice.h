// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2025 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_MATRIX_SLICE_H
#define DIRECTIONAL_MATRIX_SLICE_H

#include <vector>
#include <set>
#include <Eigen/Core>

namespace directional{

//Slicing a matrix by rowIndices and colIndices, and outputting the results in B.
template<typename T>
inline void matrix_slice(const Eigen::Matrix<T, Eigen::Dynamic,  Eigen::Dynamic> A,
                         const Eigen::VectorXi& rowIndices,
                         const Eigen::VectorXi& colIndices,
                         Eigen::Matrix<T, Eigen::Dynamic,  Eigen::Dynamic>& B){
    
    B.resize(rowIndices.size(), colIndices.size());
    for (int i=0;i<rowIndices.size();i++)
        for (int j=0;j<colIndices.size();j++)
            B(i,j) = A(rowIndices(i), colIndices(j));
    
}

}

#endif //DIRECTIONAL_TUTORIALS_INDEX_OPERATIONS_H
