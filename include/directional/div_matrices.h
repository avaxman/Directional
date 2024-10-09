// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_DIV_MATRICES_H
#define DIRECTIONAL_DIV_MATRICES_H

#include <eigen/sparse>
#include <directional/TriMesh.h>
#include <directional/single_to_N_matrix.h>
#include <directional/mass_matrices.h>
#include <directional/gradient_matrices.h>
#include <directional/mass_matrices.h>

namespace directional {

    template<typename NumberType>
    Eigen::SparseMatrix<NumberType> div_matrix_2D(const TriMesh* mesh,
                                                  const bool isIntrinsic,
                                                  const int N,
                                                  const int d){

        assert("This method is currently defined only for d==1" && d==1);
        Eigen::SparseMatrix<NumberType> G = directional::conf_gradient_matrix_2D<NumberType>(mesh,isIntrinsic,N,d);
        Eigen::SparseMatrix<NumberType> Mx = directional::face_vectors_mass_matrix_2D<NumberType>(mesh, isIntrinsic, N, d);
        return (G.adjoint()*Mx);
    }
}

#endif