// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_GAUSSIAN_CURVATURE_H
#define DIRECTIONAL_GAUSSIAN_CURVATURE_H

#include <directional/definitions.h>
#include <Eigen/Core>


namespace directional
{
    // Computes boundary-aware discrete Gaussian curvature on vertices (angle defect).
    // Input:
    //  V:                  #V by 3 vertices.
    //  F:                  #F by 3 triangles.
    //  isBoundaryVertex:   #V boolean indicating if vertex is a boundary.
    //output:
    //  G:                  #V discrete Gaussian curvature. sum(G) = eulerChar of mesh.
    inline void gaussian_curvature(const Eigen::MatrixXd& V,
                                       const Eigen::MatrixXi& F,
                                       const Eigen::VectorXi& isBoundaryVertex,
                                       Eigen::VectorXd& G){

        G.resize(V.rows());
        for (int i=0;i<V.rows();i++)
            G(i)=(isBoundaryVertex(i) ? directional::PI : 2.0*directional::PI);

        for (int i=0;i<F.rows();i++) {
            for (int j = 0; j < 3; j++) {
                Eigen::RowVector3d v1 = V.row(F(i, (j + 1) % 3)) - V.row(F(i, j));
                Eigen::RowVector3d v2 = V.row(F(i, (j + 2) % 3)) - V.row(F(i, j));
                double currAngle = std::acos(v1.dot(v2) / (v1.norm() * v2.norm()));
                G(F(i, j)) -= currAngle;
            }
        }
    }
}

#endif //DIRECTIONAL_TUTORIALS_GAUSSIAN_CURVATURE_H
