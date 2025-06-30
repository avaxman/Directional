// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_CURL_MATRICES_H
#define DIRECTIONAL_CURL_MATRICES_H

#include <cassert>
#include <Eigen/Sparse>
#include <directional/TriMesh.h>
#include <directional/single_to_N_matrix.h>


namespace directional {

template<typename NumberType>
Eigen::SparseMatrix<NumberType> curl_matrix_2D(const TriMesh& mesh,
                                               const Eigen::VectorXi& _matching=Eigen::VectorXi(),
                                               const bool isIntrinsic=false,
                                               const int N = 1,
                                               const int d = 1){
    
    assert("This method is currently defined only for d==1" && d==1);
    Eigen::VectorXi matching;
    if (_matching.size()==0)
        matching = Eigen::VectorXi::Zero(mesh.innerEdges.size());
    else
        matching = _matching;
    
    Eigen::SparseMatrix<double> singleCurlMatrix(N*mesh.innerEdges.size(), (isIntrinsic ? 2*N*mesh.F.rows() : 3*N*mesh.F.rows()));
    std::vector<Eigen::Triplet<double>> singleCurlMatTris;
    for (int i=0;i<mesh.innerEdges.size();i++){
        Eigen::RowVector3d e = mesh.V.row(mesh.EV(mesh.innerEdges(i),1))-mesh.V.row(mesh.EV(mesh.innerEdges(i),0));
        //curl is <right_face - left_face , e>
        for (int n=0;n<N;n++){
            if (isIntrinsic) {
                Eigen::RowVector2d einLeft;
                einLeft << e.dot(mesh.FBx.row(mesh.EF(mesh.innerEdges(i), 0))),
                e.dot(mesh.FBy.row(mesh.EF(mesh.innerEdges(i), 0)));
                
                Eigen::RowVector2d einRight;
                einRight << e.dot(mesh.FBx.row(mesh.EF(mesh.innerEdges(i), 1))),
                e.dot(mesh.FBy.row(mesh.EF(mesh.innerEdges(i), 1)));
                
                singleCurlMatTris.push_back(
                                            Eigen::Triplet<double>(N * i + n, 2 * N * mesh.EF(mesh.innerEdges(i), 0) + 2*n, -einLeft(0)));
                singleCurlMatTris.push_back(
                                            Eigen::Triplet<double>(N * i + n, 2 * N * mesh.EF(mesh.innerEdges(i), 0) + 2*n + 1, -einLeft(1)));
                singleCurlMatTris.push_back(
                                            Eigen::Triplet<double>(N * i + n, 2 * N * mesh.EF(mesh.innerEdges(i), 1) + 2*((n+matching(i)+N)%N), einRight(0)));
                singleCurlMatTris.push_back(
                                            Eigen::Triplet<double>(N * i + n, 2 * N * mesh.EF(mesh.innerEdges(i), 1) + 2*((n+matching(i)+N)%N) + 1, einRight(1)));
            } else {
                singleCurlMatTris.push_back(
                                            Eigen::Triplet<double>(N * i + n, 3 * N * mesh.EF(mesh.innerEdges(i), 0) + 3*n, -e(0)));
                singleCurlMatTris.push_back(
                                            Eigen::Triplet<double>(N * i + n, 3 * N * mesh.EF(mesh.innerEdges(i), 0) + 3*n + 1, -e(1)));
                singleCurlMatTris.push_back(
                                            Eigen::Triplet<double>(N * i + n, 3 * N * mesh.EF(mesh.innerEdges(i), 0) + 3*n + 2, -e(2)));
                singleCurlMatTris.push_back(
                                            Eigen::Triplet<double>(N * i + n, 3 * N * mesh.EF(mesh.innerEdges(i), 1) + 3*((n+matching(i)+N)%N), e(0)));
                singleCurlMatTris.push_back(
                                            Eigen::Triplet<double>(N * i + n, 3 * N * mesh.EF(mesh.innerEdges(i), 1) + 3*((n+matching(i)+N)%N) + 1, e(1)));
                singleCurlMatTris.push_back(
                                            Eigen::Triplet<double>(N * i + n, 3 * N * mesh.EF(mesh.innerEdges(i), 1) + 3*((n+matching(i)+N)%N) + 2, e(2)));
            }
        }
    }
    
    singleCurlMatrix.setFromTriplets(singleCurlMatTris.begin(), singleCurlMatTris.end());
    
    //if (N==1)
    return singleCurlMatrix;
    //else return single_to_N_matrix(singleCurlMatrix, N, 1, (isIntrinsic ? 2 : 3));
}
}

#endif
