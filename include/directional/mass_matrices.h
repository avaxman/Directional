// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_MASS_MATRICES_H
#define DIRECTIONAL_MASS_MATRICES_H

#include <Eigen/Sparse>
#include <directional/TriMesh.h>
#include <directional/single_to_N_matrix.h>
#include <directional/sparse_diagonal.h>

namespace directional {

//This file concentrates several methods that return a sparse mass matrix targeting a specific FEM or DEC function space inner product.
//The common inputs:
//mesh:     The input mesh
//dim:        The dimension of the quantity in question (1 = scalar, 3 = 3D vector etc.)
//deg:        The polynomial degree/order of the space. Currently still TBD, so only accepting d = 1 (linear elements).
//Output:   a symmatrix sparse mass matrix of the relevant sizing.

//The PL conforming elements with vertex dofs mass matrix (non-lumped). This is inner products of the PL elements
template<typename NumberType>
Eigen::SparseMatrix<NumberType> conf_mass_matrix_2D(const TriMesh& mesh,
                                                    const int dim = 1,
                                                    const int deg = 1) {
    assert("Currently only implemented for d=1" && d == 1);
    Eigen::SparseMatrix<NumberType> M1(mesh.V.rows(), mesh.V.rows());
    std::vector<Eigen::Triplet<NumberType>> MTris;
    for (int i = 0; i < mesh.F.rows(); i++) {
        for (int j = 0; j < 3; j++)
            for (int k = 0; k < 3; k++)
                MTris.push_back(Eigen::Triplet<NumberType>(mesh.F(i, j), mesh.F(i, k),
                                                           mesh.faceAreas(i) * (j == k ? 1.0 / 6.0 : 1.0 / 12.0)));
    }
    M1.setFromTriplets(MTris.begin(), MTris.end());
    if (N==1)
        return M1;
    return single_to_N_matrix(M1, dim, 1, 1);
}


//The (popular) lumped version of the PL conforming mass matrix (for d=1) matrix. Basically just voronoi area masses
template<typename NumberType>
Eigen::SparseMatrix<NumberType> lumped_voronoi_mass_matrix_2D(const TriMesh& mesh,
                                                              const int dim = 1){
    
    Eigen::SparseMatrix<NumberType> M1(mesh.V.rows(), mesh.V.rows());
    std::vector<Eigen::Triplet<NumberType>> MTris;
    for (int i = 0; i < mesh.F.rows(); i++) {
        //adding the 1/3 of each face's area to the vertex
        for (int j = 0; j < 3; j++)
            MTris.push_back(Eigen::Triplet<NumberType>(mesh.F(i, j), mesh.F(i, j), mesh.faceAreas(i) / 3.0));
    }
    M1.setFromTriplets(MTris.begin(), MTris.end());
    if (N==1)
        return M1;
    return single_to_N_matrix(M1, dim, 1, 1);
}

//The mass matrix for piecewise face-based vector quantities, which is basically face areas
//The inverse allows to use inverse face areas in case the quantity is an integrated one.
//the default dim is 3 because this is usually used for face-based vectors
template<typename NumberType>
Eigen::SparseMatrix<NumberType> face_mass_matrix_2D(const TriMesh& mesh,
                                                    const bool isInverse = false,
                                                    const int dim=3,
                                                    const int deg=0){
    assert("Currently only works for deg==0" && deg==0);
    Eigen::SparseMatrix<NumberType> M1 = sparse_diagonal((isInverse ? 1.0/(mesh.faceAreas.array()) : mesh.faceAreas));
    return single_to_N_matrix(M1, dim, 1, 1);
}

//The mass matrix for face-based scalar quantities, which are *inverse* area by default
/*template<typename NumberType>
 Eigen::SparseMatrix<NumberType> face_scalar_mass_matrix_2D(const TriMesh& mesh,
 const bool isInverse = false,
 const int N=1,
 const int d=1){
 assert("Currently only works for d==1" && d==1);
 Eigen::SparseMatrix<NumberType> M1 = sparse_diagonal((isInverse ? mesh.faceAreas : 1.0/(mesh.faceAreas.array())));
 return single_to_N_matrix(M1, N, 1, 1);
 }*/


//Mass matrix for quantities on edge diamond regions, which is the area of the diamond
template<typename NumberType>
Eigen::SparseMatrix<NumberType> edge_diamond_mass_matrix_2D(const TriMesh& mesh,
                                                            const bool isInverse = false;
                                                            const int deg=1,
                                                            const int dim=1){
    assert(dim==1 && "Currently only works for dim==1");
    Eigen::SparseMatrix<NumberType> M1(mesh.EV.rows(), mesh.EV.rows());
    std::vector<Eigen::Triplet<NumberType>> MTris;
    for (int i = 0; i < mesh.EF.rows(); i++) {
        //adding the 1/3 of each face's area to the edge, and taking an inverse
        NumberType mass=0.0;
        if (mesh.EF(i,0)!=-1) mass+=mesh.faceAreas(mesh.EF(i,0))/3.0;
        if (mesh.EF(i,1)!=-1) mass+=mesh.faceAreas(mesh.EF(i,1))/3.0;
        MTris.push_back(Eigen::Triplet<NumberType>(i,i, (isInverse ? 1.0/mass : mass)));
    }
    M1.setFromTriplets(MTris.begin(), MTris.end());
    if (deg==1)
        return M1;
    return single_to_N_matrix(M1, deg, 1, 1);
}

//This is not a mass matrix, but rather the "J" rotation operator around the normal of each face
//The input here is N, the usual degree of the field (number of vectors), because the choices are fundamentally different
template<typename NumberType>
Eigen::SparseMatrix<NumberType> face_vector_rotation_matrix_2D(const TriMesh& mesh,
                                                               const bool isIntrinsic = false,
                                                               const int N = 1,
                                                               const int d = 1){
    assert("Currently only works for d==1" && d==1);
    Eigen::SparseMatrix<NumberType> J(2*mesh.F.rows(), 2*mesh.F.rows());
    std::vector<Eigen::Triplet<NumberType>> JTris;
    for (int i=0;i<mesh.F.rows();i++){
        JTris.push_back(Eigen::Triplet<NumberType>(2*i,2*i, 0.0));
        JTris.push_back(Eigen::Triplet<NumberType>(2*i+1,2*i+1, 0.0));
        JTris.push_back(Eigen::Triplet<NumberType>(2*i+1,2*i, 1.0));
        JTris.push_back(Eigen::Triplet<NumberType>(2*i,2*i+1, -1.0));
    }
    J.setFromTriplets(JTris.begin(), JTris.end());
    if (!isIntrinsic){
        Eigen::SparseMatrix<NumberType> EI = face_extrinsic_to_intrinsic_matrix_2D<NumberType>(mesh, N, d);
        if (N==1)
            return EI.adjoint()*J*EI;
        return EI.adjoint()*single_to_N_matrix(J, N, 2, 2)*EI;
    } else {
        if (N==1)
            return J;
        return single_to_N_matrix(J, N, 2, 2);
    }
}

}

#endif //DIRECTIONAL_MASS_MATRICES_H
