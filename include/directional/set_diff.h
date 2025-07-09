// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2025 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_SET_DIFF_H
#define DIRECTIONAL_SET_DIFF_H

#include <vector>
#include <set>
#include <Eigen/Core>

namespace directional{

//Set Difference between two vectors A/B, putting the result in C, and outputting in IA the indices into A of all elements in B.
template<typename T>
inline void set_diff(const Eigen::Matrix<T, Eigen::Dynamic,  1>& A,
                     const Eigen::Matrix<T, Eigen::Dynamic,  1>& B,
                     Eigen::Matrix<T, Eigen::Dynamic,  1>& C,
                     Eigen::VectorXi& IA) {
    
    
    std::vector<T> CList;
    std::vector<int> IAList;
    std::set<T> BSet(B.begin(), B.end());
    for (int i = 0; i < A.size(); i++) {
        if (BSet.find(A[i]) != BSet.end())  //found in both A and B and thus not included
            continue;
        
        CList.push_back(A[i]);
        IAList.push_back(i);
    }
    C.resize(CList.size());
    IA.resize(IAList.size());
    std::copy(CList.begin(), CList.end(), C.data());
    std::copy(IAList.begin(), IAList.end(), IA.data());
}
}

#endif
