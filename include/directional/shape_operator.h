// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_SHAPE_OPERATOR_H
#define DIRECTIONAL_SHAPE_OPERATOR_H


#include <Eigen/Core>
#include <vector>
#include <directional/TriMesh.h>

/**********Computing principal curvatures and directions by estimating per-vertex shape operators
 Based on "Smooth Feature Lines on Surface Meshes" by Hildebrandt et al. 2005
 ***********************************************************************************/

//TODO: find out what boundaries should actually do.
namespace directional
{
    inline void shape_operator(const TriMesh& mesh,
                               std::vector<Eigen::Matrix3d>& Se,
                               std::vector<Eigen::Matrix3d>& Sv,
                               std::vector<Eigen::Matrix3d>& Sf)

    {
        using namespace Eigen;
        MatrixXd e(mesh.EV.rows(),3);
        MatrixXd Ne(mesh.EV.rows(),3);
        VectorXd cosHalfTheta(mesh.EV.rows());
        Se.resize(mesh.EV.rows());
        Sv.resize(mesh.V.rows());
        Sf.resize(mesh.F.rows());
        VectorXd He(mesh.EV.rows());

        for (int i=0;i<mesh.V.rows();i++)
            Sv[i] = Matrix3d::Zero();

        for (int i=0;i<mesh.EV.rows();i++) {
            e.row(i) = mesh.V.row(mesh.EV(i, 1)) - mesh.V.row(mesh.EV(i, 0));
            if (mesh.EF(i,1)!=-1)
                Ne.row(i) = (mesh.faceNormals.row(mesh.EF(i,0))+mesh.faceNormals.row(mesh.EF(i,1))).normalized();
            else
                Ne.row(i) = mesh.faceNormals.row(mesh.EF(i,0));

            //Assuming it's always positive?
            cosHalfTheta(i) = (Ne.row(i).cross(mesh.faceNormals.row(mesh.EF(i,0)))).norm();
            He(i) = 2* e.row(i).norm() * cosHalfTheta(i);
            RowVector3d eNe = e.row(i).normalized().cross(Ne.row(i));
            Se[i] = He * eNe.transpose() * eNe;
            Sv[mesh.EV(i,0)] += 0.5 * Se[i] * Ne.row(i).dot(mesh.vertexNormals.row(mesh.EV(i,0)));
            Sv[mesh.EV(i,1)] += 0.5 * Se[i] * Ne.row(i).dot(mesh.vertexNormals.row(mesh.EV(i,1)));

            Sf[mesh.EF(i,0)] += 0.5 * Se[i] * Ne.row(i).dot(mesh.faceNormals.row(mesh.EV(i,0)));
            Sf[mesh.EF(i,1)] += 0.5 * Se[i] * Ne.row(i).dot(mesh.faceNormals.row(mesh.EV(i,1)));

        }

    }
}




#endif


