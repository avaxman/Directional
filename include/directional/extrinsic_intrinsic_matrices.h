// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_EXTRINSIC_INTRINSIC_MATRICES_H
#define DIRECTIONAL_EXTRINSIC_INTRINSIC_MATRICES_H

#include <Eigen/Sparse>
#include <directional/TriMesh.h>
#include <directional/single_to_N_matrix.h>

namespace directional {

//Computing a matrix that transforms an xyzxyz extrinsic field field in raw format to an xyxyxy intrinsic vector.
//Note that this operation is lossy: any normal components in the extrinsic field will be lost.
//Input:
//  Mesh:     a triangle mesh
//  N:        the degree of the field
//  d:        the polynomial degree (inactive)
//Output:
//  A 2N x 3N conversion matrix.

template<typename NumberType>
Eigen::SparseMatrix<NumberType> face_extrinsic_to_intrinsic_matrix_2D(const TriMesh& mesh,
                                                                      const int N = 1,
                                                                      const int d = 1) {
    
    assert(d == 1 && "This method is currently defined only for d==1");
    Eigen::SparseMatrix<NumberType> EI(2 * mesh.F.rows(), 3 * mesh.F.rows());
    std::vector<Eigen::Triplet<NumberType>> EITris;
    for (int i = 0; i < mesh.F.rows(); i++) {
        EITris.push_back(Eigen::Triplet<NumberType>(2 * i, 3 * i, mesh.FBx(i, 0)));
        EITris.push_back(Eigen::Triplet<NumberType>(2 * i, 3 * i + 1, mesh.FBx(i, 1)));
        EITris.push_back(Eigen::Triplet<NumberType>(2 * i, 3 * i + 2, mesh.FBx(i, 2)));
        EITris.push_back(Eigen::Triplet<NumberType>(2 * i + 1, 3 * i, mesh.FBy(i, 0)));
        EITris.push_back(Eigen::Triplet<NumberType>(2 * i + 1, 3 * i + 1, mesh.FBy(i, 1)));
        EITris.push_back(Eigen::Triplet<NumberType>(2 * i + 1, 3 * i + 2, mesh.FBy(i, 2)));
    }
    EI.setFromTriplets(EITris.begin(), EITris.end());
    return (N == 1 ? EI : directional::single_to_N_matrix(EI, N, 2, 3));
}

//This is the opposite matrix. Note that since extrinsic to intrinsic is lossy, it's just the adjoint of the above matrix.
template<typename NumberType>
Eigen::SparseMatrix<NumberType> face_intrinsic_to_extrinsic_matrix_2D(const TriMesh& mesh,
                                                                      const int N = 1,
                                                                      const int d = 1) {
    assert(d == 1 && "This method is currently defined only for d==1");
    return face_extrinsic_to_intrinsic_matrix_2D<NumberType>(mesh, N, d).adjoint();
}
}

#endif
