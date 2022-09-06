// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_GAUSSIAN_CURVATURE_H
#define DIRECTIONAL_GAUSSIAN_CURVATURE_H

#include <igl/PI.h>
#include <Eigen/Core>


namespace directional
{
    // Creates the set of independent dual cycles (closed loops of connected faces that cannot be morphed to each other) on a mesh. Primarily used for index prescription.
    // The basis cycle matrix first contains #V-#b cycles for every inner vertex (by order), then #b boundary cycles, and finally 2*g generator cycles around all handles. Total #c cycles.The cycle matrix sums information on the dual edges between the faces, and is indexed into the inner edges alone (excluding boundary)
    //input:
    //  V: #V by 3 vertices.
    //  F: #F by 3 triangles.
    //  EV: #E by 2 matrix of edges (vertex indices)
    //  EF: #E by 2 matrix of oriented adjacent faces
    //output:
    //  basisCycles:    #c by #iE basis cycles
    //  cycleCurvature:   #c by 1 curvatures of each cycle (for inner-vertex cycles, simply the Gaussian curvature.
    //  vertex2cycle:     #v by 1 map between vertex and corresponding cycle (for comfort of input from the user's side; inner vertices map to their cycles, boundary vertices to the bigger boundary cycle.
    //  innerEdges:       #iE by 1 the subset of #EV that are inner edges, and with the same ordering as the columns of basisCycles.

    IGL_INLINE void gaussian_curvature(const Eigen::MatrixXd& V,
                                       const Eigen::MatrixXi& F,
                                       const Eigen::VectorXi& isBoundaryVertex,
                                       Eigen::VectorXd& G){

        G.resize(V.rows());
        for (int i=0;i<V.rows();i++)
            G(i)=(isBoundaryVertex(i) ? igl::PI : 2.0*igl::PI);

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
