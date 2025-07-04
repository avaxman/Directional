// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2025 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_SHAPE_OPERATOR_H
#define DIRECTIONAL_SHAPE_OPERATOR_H


#include <Eigen/Core>
#include <vector>

/***Computing principal curvatures and directions by estimating per-vertex shape operators
 Based on "Smooth Feature Lines on Surface Meshes" by Hildebrandt et al. 2005
 Input:
 V:          matrix |V|x3 of vertices
 F:          matrix of vertex indices indicating faces
 EV:         |E|x2 of edges (vertex indices)
 EF:         |E|x2 left and right faces of each edge (-1 where boundary).
 faceNormals:        |F|x3 normalized face normals
 vertexNormals:      |V|x3 normalized vertex normals
 Output:
 Se:         vector of 3x3 shape operators per edge
 Sf:           the same for vertices
 Sv:         the same for vertices
 ***/

namespace directional
{
inline void shape_operator(const Eigen::MatrixXd& V,
                           const Eigen::MatrixXi& F,
                           const Eigen::MatrixXi& EV,
                           const Eigen::MatrixXi& EF,
                           const Eigen::MatrixXd& faceNormals,
                           const Eigen::MatrixXd& vertexNormals,
                           std::vector<Eigen::Matrix3d>& Se,
                           std::vector<Eigen::Matrix3d>& Sv,
                           std::vector<Eigen::Matrix3d>& Sf)

{
    using namespace Eigen;
    MatrixXd e(EV.rows(),3);
    MatrixXd Ne(EV.rows(),3);
    Se.resize(EV.rows());
    Sv.resize(V.rows());
    Sf.resize(F.rows());
  
    for (int i=0;i<V.rows();i++)
        Sv[i] = Matrix3d::Zero();
    
    for (int i=0;i<EV.rows();i++) {
        e.row(i) = V.row(EV(i, 1)) - V.row(EV(i, 0));
        if (EF(i,1)!=-1)
            Ne.row(i) = (faceNormals.row(EF(i,0))+faceNormals.row(EF(i,1))).normalized();
        else
            Ne.row(i) = faceNormals.row(EF(i,0));
        
        //Assuming it's always positive?
        RowVector3d n = faceNormals.row(EF(i,0)).eval();
        RowVector3d ne = Ne.row(i);
        double He = 2*  (n.cross(ne)).dot(e.row(i)); //e.row(i).norm() * cosHalfTheta(i);
        RowVector3d evec = e.row(i).normalized();
        RowVector3d eNe = evec.cross(ne);
        Se[i] = He * eNe.transpose() * eNe;
        Sv[EV(i,0)] += 0.5 * Se[i] * (ne.dot(vertexNormals.row(EV(i,0))));
        Sv[EV(i,1)] += 0.5 * Se[i] * (ne.dot(vertexNormals.row(EV(i,1))));
    }
    
    //averaging operator to faces
    for (int i=0;i<F.rows();i++){
        SelfAdjointEigenSolver<Matrix3d> esA(Sv[F(i,0)]), esB(Sv[F(i,1)]), esC(Sv[F(i,2)]);
        Matrix3d Qsum = esA.eigenvectors() + esB.eigenvectors() + esC.eigenvectors();
        JacobiSVD<Matrix3d> svd(Qsum, ComputeFullU | ComputeFullV);
        Matrix3d Qavg = svd.matrixU() * svd.matrixV().transpose();
        Vector3d eval_avg = (esA.eigenvalues() + esB.eigenvalues() + esC.eigenvalues()) / 3.0;
        Sf[i] = Qavg * eval_avg.asDiagonal() * Qavg.transpose();
    }
}
}




#endif


