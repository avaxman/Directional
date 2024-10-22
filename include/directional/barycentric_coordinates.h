// This file is part of libigl, a simple c++ geometry processing library.
// 
// Copyright (C) 2014 Alec Jacobson <alecjacobson@gmail.com>
// 
// This Source Code Form is subject to the terms of the Mozilla Public License 
// v. 2.0. If a copy of the MPL was not distributed with this file, You can 
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_BARYCENTRIC_COORDINATES_H
#define IGL_BARYCENTRIC_COORDINATES_H


#include <Eigen/Core>
namespace directional {

    template<
            typename DerivedP,
            typename DerivedA,
            typename DerivedB,
            typename DerivedC,
            typename DerivedL>
    inline void barycentric_coordinates(
            const Eigen::MatrixBase <DerivedP> &P,
            const Eigen::MatrixBase <DerivedA> &A,
            const Eigen::MatrixBase <DerivedB> &B,
            const Eigen::MatrixBase <DerivedC> &C,
            Eigen::PlainObjectBase <DerivedL> &L) {
        using namespace Eigen;
#ifndef NDEBUG
        const int DIM = P.cols();
        assert(A.cols() == DIM && "corners must be in same dimension as query");
        assert(B.cols() == DIM && "corners must be in same dimension as query");
        assert(C.cols() == DIM && "corners must be in same dimension as query");
        assert(P.rows() == A.rows() && "Must have same number of queries as corners");
        assert(A.rows() == B.rows() && "Corners must be same size");
        assert(A.rows() == C.rows() && "Corners must be same size");
#endif

        // http://gamedev.stackexchange.com/a/23745
        typedef
        Eigen::Array<
                typename DerivedP::Scalar,
                DerivedP::RowsAtCompileTime,
                DerivedP::ColsAtCompileTime>
                ArrayS;
        typedef
        Eigen::Array<
                typename DerivedP::Scalar,
                DerivedP::RowsAtCompileTime,
                1>
                VectorS;

        const ArrayS v0 = B.array() - A.array();
        const ArrayS v1 = C.array() - A.array();
        const ArrayS v2 = P.array() - A.array();
        VectorS d00 = (v0 * v0).rowwise().sum();
        VectorS d01 = (v0 * v1).rowwise().sum();
        VectorS d11 = (v1 * v1).rowwise().sum();
        VectorS d20 = (v2 * v0).rowwise().sum();
        VectorS d21 = (v2 * v1).rowwise().sum();
        VectorS denom = d00 * d11 - d01 * d01;
        L.resize(P.rows(), 3);
        L.col(1) = (d11 * d20 - d01 * d21) / denom;
        L.col(2) = (d00 * d21 - d01 * d20) / denom;
        L.col(0) = 1.0f - (L.col(1) + L.col(2)).array();
    }
}


#endif
