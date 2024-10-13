// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_CURL_MATRICES_H
#define DIRECTIONAL_CURL_MATRICES_H

#include <eigen/sparse>
#include <directional/TriMesh.h>
#include <directional/single_to_N_matrix.h>


namespace directional {

    template<typename NumberType>
    Eigen::SparseMatrix<NumberType> curl_matrix_2D(const TriMesh* mesh,
                                                   const bool isIntrinsic=false,
                                                   const int N = 1,
                                                   const int d = 1){

        assert("This method is currently defined only for d==1" && d==1);
        Eigen::SparseMatrix<double> singleCurlMatrix(mesh->innerEdges.size(), (isIntrinsic ? 2*mesh->F.rows() : 3*mesh->F.rows()));
        std::vector<Eigen::Triplet<double>> singleCurlMatTris;
        for (int i=0;i<mesh->innerEdges.size();i++){
            Eigen::RowVector3d e = mesh->V.row(mesh->EV(mesh->innerEdges(i),1))-mesh->V.row(mesh->EV(mesh->innerEdges(i),0));
            //curl is <right_face - left_face , e>
            if (isIntrinsic) {
                Eigen::RowVector2d einLeft;
                einLeft << e.dot(mesh->FBx.row(mesh->EF(mesh->innerEdges(i), 0))),
                        e.dot(mesh->FBy.row(mesh->EF(mesh->innerEdges(i), 0)));

                Eigen::RowVector2d einRight;
                einRight << e.dot(mesh->FBx.row(mesh->EF(mesh->innerEdges(i), 1))),
                        e.dot(mesh->FBy.row(mesh->EF(mesh->innerEdges(i), 1)));

                singleCurlMatTris.push_back(
                        Eigen::Triplet<double>(i, 2 * mesh->EF(mesh->innerEdges(i), 0), -einLeft(0)));
                singleCurlMatTris.push_back(
                        Eigen::Triplet<double>(i, 2 * mesh->EF(mesh->innerEdges(i), 0) + 1, -einLeft(1)));
                singleCurlMatTris.push_back(
                        Eigen::Triplet<double>(i, 2 * mesh->EF(mesh->innerEdges(i), 1), einRight(0)));
                singleCurlMatTris.push_back(
                        Eigen::Triplet<double>(i, 2 * mesh->EF(mesh->innerEdges(i), 1) + 1, einRight(1)));
            } else {
                singleCurlMatTris.push_back(
                        Eigen::Triplet<double>(i, 3 * mesh->EF(mesh->innerEdges(i), 0), -e(0)));
                singleCurlMatTris.push_back(
                        Eigen::Triplet<double>(i, 3 * mesh->EF(mesh->innerEdges(i), 0) + 1, -e(1)));
                singleCurlMatTris.push_back(
                        Eigen::Triplet<double>(i, 3 * mesh->EF(mesh->innerEdges(i), 0) + 2, -e(2)));
                singleCurlMatTris.push_back(
                        Eigen::Triplet<double>(i, 3 * mesh->EF(mesh->innerEdges(i), 1), e(0)));
                singleCurlMatTris.push_back(
                        Eigen::Triplet<double>(i, 3 * mesh->EF(mesh->innerEdges(i), 1) + 1, e(1)));
                singleCurlMatTris.push_back(
                        Eigen::Triplet<double>(i, 3 * mesh->EF(mesh->innerEdges(i), 1) + 2, e(2)));
            }
        }

        singleCurlMatrix.setFromTriplets(singleCurlMatTris.begin(), singleCurlMatTris.end());

        if (N==1)
            return singleCurlMatrix;
        else return single_to_N_matrix(singleCurlMatrix, N, 1, (isIntrinsic ? 2 : 3));
    }
}

#endif