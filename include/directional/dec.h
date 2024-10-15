// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_DEC_H
#define DIRECTIONAL_DEC_H

#include <eigen/sparse>
#include <directional/TriMesh.h>

namespace directional {

    template<typename NumberType>
    Eigen::SparseMatrix<NumberType> d0_matrix(const TriMesh *mesh,
                                              const bool isDirichlet=false) {
        Eigen::SparseMatrix<NumberType> d0(mesh->EV.rows(), mesh->V.rows());
        std::vector<Eigen::Triplet<NumberType>> d0Tris;
        for (int i=0;i<mesh->EV.rows();i++) {
            d0Tris.push_back(Eigen::Triplet<NumberType>(i, mesh->EV(i, 0), -1.0));
            d0Tris.push_back(Eigen::Triplet<NumberType>(i, mesh->EV(i, 1), 1.0));
        }
        d0.setFromTriplets(d0Tris.begin(), d0Tris.end());
        return d0;
    }

    template<typename NumberType>
    Eigen::SparseMatrix<NumberType> d1_matrix(const TriMesh *mesh,
                                              const bool isDirichlet=false) {
        Eigen::SparseMatrix<NumberType> d1(mesh->F.rows(), mesh->EF.rows());
        std::vector<Eigen::Triplet<NumberType>> d1Tris;
        for (int i=0;i<mesh->FE.rows();i++)
            for (int j=0;j<3;j++)
                d1Tris.push_back(Eigen::Triplet<NumberType>(i, mesh->FE(i, j), mesh->FEs(i,j)));

        d1.setFromTriplets(d1Tris.begin(), d1Tris.end());
        return d1;
    }

    //The primal/dual diagonal hodge star (with a choice of center so as to make it positive)
    template<typename NumberType>
    void hodge_star_1_matrix(const TriMesh *mesh,
                             Eigen::SparseMatrix<NumberType>& hodgeStar,
                             Eigen::SparseMatrix<NumberType>& invHodgeStar,
                             const bool circumcenter = true,
                             const bool isDirichlet = false) {

        /*Eigen::VectorXd primalEdgeLengths(mesh->EV.rows());
        Eigen::VectorXd dualEdgeLengths(mesh->EV.rows());*/
        Eigen::VectorXd M1Weights(mesh->EV.rows());

        for (int i=0;i<mesh->EV.rows();i++){
            if (!circumcenter) {
                double primalEdgeLengths = (mesh->V.row(mesh->EV(i, 0)) - mesh->V.row(mesh->EV(i, 1))).norm();
                double dualEdgeLength=0.0;
                Eigen::RowVector3d midEdge = (mesh->V.row(mesh->EV(i, 0)) + mesh->V.row(mesh->EV(i, 1))) / 2.0;
                bool isBoundary = false;
                if (mesh->EF(i,0)!=-1) {
                    isBoundary = true;
                    Eigen::RowVector3d leftBarycenter = mesh->barycenters.row(mesh->EF(i, 0));
                    dualEdgeLength += (leftBarycenter - midEdge).norm();
                }
                if (mesh->EF(i,1)!=-1) {
                    isBoundary = true;
                    Eigen::RowVector3d rightBarycenter = mesh->barycenters.row(mesh->EF(i, 1));
                    dualEdgeLength += (rightBarycenter - midEdge).norm();
                }
                M1Weights(i)=(isBoundary ? 0.5 : 1.0)*dualEdgeLength/primalEdgeLengths;
            } else {//cot weights
                M1Weights(i)=0.0;
                for (int j=0;j<2;j++){
                    if (mesh->EF(i,j)==-1)
                        continue;
                    int nextHe = mesh->nextH(mesh->EH(i,j));
                    int prevHe = mesh->prevH(mesh->EH(i,j));

                    Eigen::RowVector3d e1 = -(mesh->V.row(mesh->HV(mesh->nextH(nextHe)))-mesh->V.row(mesh->HV(nextHe)));
                    Eigen::RowVector3d e2 = mesh->V.row(mesh->HV(mesh->nextH(prevHe)))-mesh->V.row(mesh->HV(prevHe));
                    double sinAngle = e1.cross(e2).norm();
                    double cosAngle = e1.dot(e2);
                    if (std::abs(sinAngle)<10e-7)
                        continue;  //using 0 weight

                    M1Weights(i)+=0.5*cosAngle/sinAngle;
                }
            }
        }

        hodgeStar = directional::sparse_diagonal(M1Weights);
        //Note! inverse hodge has a minus sign (see here, Fig. 2: https://dl.acm.org/doi/pdf/10.1145/2897824.2925880)
        Eigen::VectorXd invWeights = M1Weights.unaryExpr([](double v) {
            return (std::abs(v)>10e-7 ? -1.0/v : 0.0);
        });
        invHodgeStar = directional::sparse_diagonal(invWeights);
    }


}

#endif //DIRECTIONAL_DEC_H
