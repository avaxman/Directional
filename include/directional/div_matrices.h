// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_DIV_MATRICES_H
#define DIRECTIONAL_DIV_MATRICES_H

#include <Eigen/Sparse>
#include <directional/TriMesh.h>
#include <directional/single_to_N_matrix.h>
#include <directional/mass_matrices.h>
#include <directional/gradient_matrices.h>
#include <directional/mass_matrices.h>

namespace directional {

//Produces the divergence operator on vector fields, as the adjoint of the conforming face-based gradient matrix.
//Input:
//mesh:         a triangle mesh
//isIntrinsic:  whether the input field is intrinsic (2D) or extrinsic (3D) in every face
//N:            the order of the field
//d:            the polynomial order (inactive)
//Output:
//A |V|x2N if intrinsic or |V|x3N if extrinsic divergence matrix.
template<typename NumberType>
Eigen::SparseMatrix<NumberType> div_matrix_2D(const TriMesh& mesh,
                                              const bool isIntrinsic = false,
                                              const int N = 1,
                                              const int d = 1){
    
    assert(d==1 && "This method is currently defined only for d==1");
    Eigen::SparseMatrix<NumberType> G = directional::conf_gradient_matrix_2D<NumberType>(mesh,isIntrinsic,N,d);
    Eigen::SparseMatrix<NumberType> Mx = directional::face_mass_matrix_2D<NumberType>(mesh, false, (isIntrinsic ? 2*N : 3*N), d);
    return (-G.adjoint()*Mx);  //Note the sign! True divergence is like that. See tutorial for explanation.
}
}

#endif
