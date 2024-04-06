//
// Created by Amir Vaxman on 03.04.24.
//

#ifndef DIRECTIONAL_TUTORIALS_INDEX_OPERATIONS_H
#define DIRECTIONAL_TUTORIALS_INDEX_OPERATIONS_H

#include <vector>
#include <set>
#include <Eigen/Core>

namespace directional{

    /// Set difference of elements of matrices
    ///
    /// @param[in] A  m-long vector of indices
    /// @param[in] B  n-long vector of indices
    /// @param[out] C  (k<=m)-long vector of unique elements appearing in A but not in B
    /// @param[out] IA  (k<=m)-long list of indices into A so that C = A(IA)

    template<typename T>
    inline void setdiff(const Eigen::Matrix<T, Eigen::Dynamic,  1>& A,
                        const Eigen::Matrix<T, Eigen::Dynamic,  1>& B,
                        Eigen::Matrix<T, Eigen::Dynamic,  1>& C,
                        Eigen::VectorXi& IA) {


        std::vector<T> CList;
        std::vector<int> IAList;
        std::set<T> BSet(B.begin(), B.end());
        for (int i = 0; i < A.size(); i++) {
            if (BSet.find(A[i]) == BSet.end())
                continue;

            CList.push_back(A[i]);
            IAList.push_back(i);
        }
        C.resize(CList.size());
        IA.resize(IAList.size());
        std::copy(CList.begin(), CList.end(), C.data());
        std::copy(IAList.begin(), IAList.end(), IA.data());
    }


    template<typename T>
    inline void slice(const Eigen::Matrix<T, Eigen::Dynamic,  Eigen::Dynamic> A,
                      const Eigen::VectorXi& rowIndices,
                      const Eigen::VectorXi& colIndices,
                      Eigen::Matrix<T, Eigen::Dynamic,  Eigen::Dynamic>& B){

        B.resize(rowIndices.size(), colIndices.size());
        for (int i=0;i<rowIndices.size();i++)
            for (int j=0;j<colIndices.size();j++)
                B(i,j) = A(rowIndices(i), colIndices(j));

    }

    template<typename T>
    inline void speye(const int numRows,
                      const int numCols,
                      Eigen::SparseMatrix<T>& eyeMat){
        eyeMat.resize(numRows, numCols);
        std::vector<Eigen::Triplet<T>> eyeMatTris;
        for (int i=0;i<(numRows < numCols ? numRows : numCols);i++)
            eyeMatTris.push_back(Eigen::Triplet<T>(i, i, 1.0));
        eyeMat.setFromTriplets(eyeMatTris.begin(), eyeMatTris.end());
    }

}

#endif //DIRECTIONAL_TUTORIALS_INDEX_OPERATIONS_H
