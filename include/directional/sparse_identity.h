// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2025 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_SPARSE_IDENTITY_H
#define DIRECTIONAL_SPARSE_IDENTITY_H

#include <vector>
#include <set>
#include <Eigen/Core>

namespace directional{

//Generating a sparse identity matrix of size (numRows, numCols), in the main diagonal.
template<typename T>
inline void sparse_identity(const int numRows,
                            const int numCols,
                            Eigen::SparseMatrix<T>& eyeMat){
    eyeMat.resize(numRows, numCols);
    std::vector<Eigen::Triplet<T>> eyeMatTris;
    for (int i=0;i<(numRows < numCols ? numRows : numCols);i++)
        eyeMatTris.push_back(Eigen::Triplet<T>(i, i, 1.0));
    eyeMat.setFromTriplets(eyeMatTris.begin(), eyeMatTris.end());
}

}

#endif
