// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_GRADIENT_MATRIX_H
#define DIRECTIONAL_GRADIENT_MATRIX_H

#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <directional/TriMesh.h>
#include <directional/sparse_diagonal.h>
#include <directional/single_to_N_matrix.h>
#include <directional/extrinsic_intrinsic_matrices.h>


namespace directional{

    //Computues the classical conforming gradient matrix from P^d polynomial scalar conforming fields to X^{d-1} piecewise polynomial fields
    //this allows for N-functions as well going to N-fields (of the appropriate polynomial degree).
    //the result is always extrinsic (i.e., does not take into account any tangent bundle intrinsic parameterization).
    template<typename NumberType>
    Eigen::SparseMatrix<NumberType> conf_gradient_matrix_2D(const TriMesh* mesh,
                                                         const bool isIntrinsic,
                                                         const int N,
                                                         const int d){

        assert("Currently only implemented for d=1" && d==1);

        Eigen::SparseMatrix<NumberType> G1(3*mesh->F.rows(), mesh->V.rows());
        std::vector<Eigen::Triplet<NumberType>> GTris;
        for (int i=0;i<mesh->F.rows();i++){
            Eigen::RowVector3d n  = mesh->faceNormals.row(i);
            Eigen::RowVector3d e01 = mesh->V.row(mesh->F(i, 1)) - mesh->V.row(mesh->F(i, 0));
            Eigen::RowVector3d e12 = mesh->V.row(mesh->F(i, 2)) - mesh->V.row(mesh->F(i, 1));
            Eigen::RowVector3d e20 = mesh->V.row(mesh->F(i, 0)) - mesh->V.row(mesh->F(i, 2));

            Eigen::Matrix3d ep; ep<<n.cross(e12), n.cross(e20), n.cross(e01);
            double faceArea = mesh->faceAreas(i);
            for (int j=0;j<3;j++)
                for (int k=0;k<3;k++) {
                    //std::cout<<"("<<3 * i + k<<","<<mesh->F(i, j)<<","<<ep(j, k) / (2.0 * faceArea)<<")"<<std::endl;
                    GTris.push_back(Eigen::Triplet<NumberType>(3 * i + k, mesh->F(i, j), ep(j, k) / (2.0 * faceArea)));
                }

        }
        G1.setFromTriplets(GTris.begin(), GTris.end());
        if (isIntrinsic){
            Eigen::SparseMatrix<NumberType> EI = directional::face_extrinsic_to_intrinsic_matrix_2D<NumberType>(mesh, N, d);
            return (N==1 ? G1 : single_to_N_matrix(G1, N, 3, 1) )*EI;
        } else return (N==1 ? G1 : single_to_N_matrix(G1, N, 3, 1));
    }

    //The non-conforming Crouzeix-Raviart gradient matrix
    template<typename NumberType>
    Eigen::SparseMatrix<NumberType> non_conf_gradient_matrix_2D(const TriMesh* mesh,
                                                             const int N,
                                                             const int d){

        assert("Currently only implemented for d=1" && d==1);
        Eigen::SparseMatrix<NumberType> G1(3*N*mesh->F.rows(), mesh->EV.rows());
        std::vector<Eigen::Triplet<NumberType>> GTris;
        for (int i=0;i<mesh->F.rows();i++){
            Eigen::RowVector3d n  = mesh->faceNormals.row(i);
            Eigen::RowVector3d e01 = mesh->V.row(mesh->F(i, 1)) - mesh->V.row(mesh->F(i, 0));
            Eigen::RowVector3d e12 = mesh->V.row(mesh->F(i, 2)) - mesh->V.row(mesh->F(i, 1));
            Eigen::RowVector3d e20 = mesh->V.row(mesh->F(i, 0)) - mesh->V.row(mesh->F(i, 2));

            Eigen::Matrix3d ep; ep<<n.cross(e01), n.cross(e12), n.cross(e20);
            double faceArea = mesh->faceAreas(i);
            for (int j=0;j<3;j++)
                for (int k=0;k<3;k++)
                    GTris.push_back(Eigen::Triplet<NumberType>(3*i+k, mesh->FE(i,j), -ep(j,k)/(faceArea)));

        }
        G1.setFromTriplets(GTris.begin(), GTris.end());
        if (N==1)
            return G1;
        else return single_to_N_matrix(G1, N, 3, 1);

    }


}

#endif //DIRECTIONAL_GRADIENT_MATRIX_H
