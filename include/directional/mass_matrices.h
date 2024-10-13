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

    template<typename NumberType>
    Eigen::SparseMatrix<NumberType> conf_mass_matrix_2D(const TriMesh* mesh,
                                                        const int N = 1,
                                                        const int d = 1) {
        assert("Currently only implemented for d=1" && d == 1);
        Eigen::SparseMatrix<NumberType> M1(mesh->V.rows(), mesh->V.rows());
        std::vector<Eigen::Triplet<NumberType>> MTris;
        for (int i = 0; i < mesh->F.rows(); i++) {
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    MTris.push_back(Eigen::Triplet<NumberType>(mesh->F(i, j), mesh->F(i, k),
                                                               mesh->faceAreas(i) * (j == k ? 1.0 / 6.0 : 1.0 / 12.0)));
        }
        M1.setFromTriplets(MTris.begin(), MTris.end());
        if (N==1)
            return M1;
        return single_to_N_matrix(M1, N, 1, 1);
    }

    //The popular lumped version of the M0 (for d=0) matrix. Just voronoi area masses
    template<typename NumberType>
    Eigen::SparseMatrix<NumberType> lumped_voronoi_mass_matrix_2D(const TriMesh* mesh,
                                                                  const int N = 1){

        Eigen::SparseMatrix<NumberType> M1(mesh->V.rows(), mesh->V.rows());
        std::vector<Eigen::Triplet<NumberType>> MTris;
        for (int i = 0; i < mesh->F.rows(); i++) {
            //adding the 1/3 of each face's area to the vertex
            for (int j = 0; j < 3; j++)
                MTris.push_back(Eigen::Triplet<NumberType>(mesh->F(i, j), mesh->F(i, j), mesh->faceAreas(i) / 3.0));
        }
        M1.setFromTriplets(MTris.begin(), MTris.end());
        if (N==1)
            return M1;
        return single_to_N_matrix(M1, N, 1, 1);
    }

    //The mass matrix for face-based vector quantities, which can be vector
    template<typename NumberType>
    Eigen::SparseMatrix<NumberType> face_vectors_mass_matrix_2D(const TriMesh* mesh,
                                                                const bool isIntrinsic = false,
                                                                const bool isInverse = false,
                                                                const int N=1,
                                                                const int d=1){
        assert("Currently only works for d==1" && d==1);
        Eigen::SparseMatrix<NumberType> M1 = sparse_diagonal((isInverse ? mesh->faceAreas : 1.0/(mesh->faceAreas.array())));
        return single_to_N_matrix(M1, (isIntrinsic ? 2 : 3)*N, 1, 1);
    }

    //Mass matrix for edge diamond regions, which is *inverse* areas
    template<typename NumberType>
    Eigen::SparseMatrix<NumberType> edge_diamond_mass_matrix_2D(const TriMesh* mesh,
                                                                const int N=1,
                                                                const int d=1){
        assert("Currently only works for d==1" && d==1);
        Eigen::SparseMatrix<NumberType> M1(mesh->EV.rows(), mesh->EV.rows());
        std::vector<Eigen::Triplet<NumberType>> MTris;
        for (int i = 0; i < mesh->EF.rows(); i++) {
            //adding the 1/3 of each face's area to the edge, and taking an inverse
            NumberType mass=0.0;
            if (mesh->EF(i,0)!=-1) mass+=mesh->faceAreas(mesh->EF(i,0))/3.0;
            if (mesh->EF(i,1)!=-1) mass+=mesh->faceAreas(mesh->EF(i,1))/3.0;
            MTris.push_back(Eigen::Triplet<NumberType>(i,i, 1.0/mass));
        }
        M1.setFromTriplets(MTris.begin(), MTris.end());
        if (N==1)
            return M1;
        return single_to_N_matrix(M1, N, 1, 1);
    }

    template<typename NumberType>
    Eigen::SparseMatrix<NumberType> face_vector_rotation_matrix_2D(const TriMesh* mesh,
                                                                   const bool isIntrinsic = false,
                                                                   const int N = 1,
                                                                   const int d = 1){
        assert("Currently only works for d==1" && d==1);
        Eigen::SparseMatrix<NumberType> J(2*mesh->F.rows(), 2*mesh->F.rows());
        std::vector<Eigen::Triplet<NumberType>> JTris;
        Eigen::SparseMatrix<NumberType> EI = face_extrinsic_to_intrinsic_matrix_2D<NumberType>(mesh, N, d);
        for (int i=0;i<mesh->F.rows();i++){
            JTris.push_back(Eigen::Triplet<NumberType>(2*i,2*i, 0.0));
            JTris.push_back(Eigen::Triplet<NumberType>(2*i+1,2*i+1, 0.0));
            JTris.push_back(Eigen::Triplet<NumberType>(2*i+1,2*i, 1.0));
            JTris.push_back(Eigen::Triplet<NumberType>(2*i,2*i+1, -1.0));
        }
        J.setFromTriplets(JTris.begin(), JTris.end());
        if (N==1)
            return EI.adjoint()*J*EI;
        return EI.adjoint()*single_to_N_matrix(J, N, 2, 2)*EI;
    }

}

#endif //DIRECTIONAL_MASS_MATRICES_H
