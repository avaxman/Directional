// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_DISCRETE_EXTERIOR_CALCULUS_H
#define DIRECTIONAL_DISCRETE_EXTERIOR_CALCULUS_H

#include <eigen/sparse>
#include <directional/TriMesh.h>


namespace directional {

//This header provide a suite of operations for low-order discrete exterior calculus. The function names are mostly self-explanating

//The output is the |E|x|V| d0 matrix (in the default Neumann boundary conditions). In case of Dirichlet boundary conditions (current inactive), the matrix is of size |E|x|Vi| where "Vi" indicates only interior vertices.
template<typename NumberType>
Eigen::SparseMatrix<NumberType> d0_matrix(const TriMesh& mesh,
                                          const bool isDirichlet=false) {
    
    assert("d0_matrix(): only Neumann boundary conditions are currently implemented"  && !isDirichlet);
    Eigen::SparseMatrix<NumberType> d0(mesh.EV.rows(), mesh.V.rows());
    std::vector<Eigen::Triplet<NumberType>> d0Tris;
    for (int i=0;i<mesh.EV.rows();i++) {
        d0Tris.push_back(Eigen::Triplet<NumberType>(i, mesh.EV(i, 0), -1.0));
        d0Tris.push_back(Eigen::Triplet<NumberType>(i, mesh.EV(i, 1), 1.0));
    }
    d0.setFromTriplets(d0Tris.begin(), d0Tris.end());
    return d0;
}

//The output is the |F|x|E| d1 matrix in the default Neumann boundary conditions. For Dirichlet boundary conditions (current inactive) the matrix would be sized |F+Fb|x|E| accepting values beyond the boundary edges
template<typename NumberType>
Eigen::SparseMatrix<NumberType> d1_matrix(const TriMesh& mesh,
                                          const bool isDirichlet=false) {
    assert("d1_matrix(): only Neumann boundary conditions are currently implemented"  && !isDirichlet);
    Eigen::SparseMatrix<NumberType> d1(mesh.F.rows(), mesh.EF.rows());
    std::vector<Eigen::Triplet<NumberType>> d1Tris;
    for (int i=0;i<mesh.FE.rows();i++)
        for (int j=0;j<3;j++)
            d1Tris.push_back(Eigen::Triplet<NumberType>(i, mesh.FE(i, j), mesh.FEs(i,j)));
    
    d1.setFromTriplets(d1Tris.begin(), d1Tris.end());
    return d1;
}

//The is M1: |E|x|E|, the linear Whitney mass matrix, which is the original mass matrix rather than the lumped one.
template<typename NumberType>
void linear_whitney_mass_matrix(const TriMesh& mesh,
                                Eigen::SparseMatrix<NumberType>& M1){
    
    std::vector<Eigen::Triplet<NumberType>> M1Tris;
    M1.resize(3*mesh.F.rows(), 3*mesh.F.rows());
    for (int i=0;i<mesh.F.rows();i++){
        Eigen::RowVector3d n  = mesh.faceNormals.row(i);
        Eigen::RowVector3d e01 = mesh.V.row(mesh.F(i, 1)) - mesh.V.row(mesh.F(i, 0));
        Eigen::RowVector3d e12 = mesh.V.row(mesh.F(i, 2)) - mesh.V.row(mesh.F(i, 1));
        Eigen::RowVector3d e20 = mesh.V.row(mesh.F(i, 0)) - mesh.V.row(mesh.F(i, 2));
        
        Eigen::Matrix3d ep; ep<<n.cross(e12), n.cross(e20), n.cross(e01);
        double faceArea = mesh.faceAreas(i);
        ep.array()/=(2*faceArea);
        //diagonal elements
        for (int j=0;j<3;j++)
            M1Tris.push_back(Eigen::Triplet<NumberType>(3*i+j, 3*i+j,
                                                        (ep.row(j).squaredNorm()+ep.row((j+1)%3).squaredNorm() -
                                                         ep.row(j).dot(ep.row((j+1)%3)))*(faceArea/6.0)));
        
        
        //off diagonal elements
        for (int j=0;j<3;j++) {
            M1Tris.push_back(Eigen::Triplet<NumberType>(3*i+j, 3*i + (j + 1) % 3,
                                                        (ep.row((j + 1) % 3).squaredNorm() +
                                                         ep.row(j).dot(ep.row((j + 2) % 3))) * (-faceArea / 6.0)));
            M1Tris.push_back(Eigen::Triplet<NumberType>(3*i + (j + 1) % 3, 3*i+j,
                                                        (ep.row((j + 1) % 3).squaredNorm() +
                                                         ep.row(j).dot(ep.row((j + 2) % 3))) * (-faceArea / 6.0)));
        }
    }
    M1.setFromTriplets(M1Tris.begin(), M1Tris.end());
    
    //conversion matrix
    std::vector<Eigen::Triplet<NumberType>> FEMatTris;
    Eigen::SparseMatrix<NumberType> FEMat(3*mesh.F.rows(), mesh.EV.rows());
    for (int i=0;i<mesh.F.rows();i++)
        for (int j=0;j<3;j++)
            FEMatTris.push_back(Eigen::Triplet<NumberType>(3*i+j, mesh.FE(i,j), mesh.FEs(i,j)));
    
    FEMat.setFromTriplets(FEMatTris.begin(), FEMatTris.end());
    M1 = FEMat.adjoint() * M1 * FEMat;
}

//The primal/dual diagonal hodge star of size |E|x|E| (with a choice of center between circle circumcenter, leading to cot weights, and barycenter, which is not linear reproducing, but is always positive), which also returns the dual hodge star that is common in DEC. Note that Hodge star is an involution and thus dual*primal = -Identity. As such we put inv = -dual.
template<typename NumberType>
void hodge_star_1_matrix(const TriMesh& mesh,
                         Eigen::SparseMatrix<NumberType>& hodgeStar,
                         Eigen::SparseMatrix<NumberType>& invHodgeStar,
                         const bool circumcenter = true,
                         const bool isDirichlet = false) {
    
    Eigen::VectorXd M1Weights(mesh.EV.rows());
    
    for (int i=0;i<mesh.EV.rows();i++){
        if (!circumcenter) {
            double primalEdgeLengths = (mesh.V.row(mesh.EV(i, 0)) - mesh.V.row(mesh.EV(i, 1))).norm();
            double dualEdgeLength=0.0;
            Eigen::RowVector3d midEdge = (mesh.V.row(mesh.EV(i, 0)) + mesh.V.row(mesh.EV(i, 1))) / 2.0;
            bool isBoundary = false;
            if (mesh.EF(i,0)!=-1) {
                isBoundary = true;
                Eigen::RowVector3d leftBarycenter = mesh.barycenters.row(mesh.EF(i, 0));
                dualEdgeLength += (leftBarycenter - midEdge).norm();
            }
            if (mesh.EF(i,1)!=-1) {
                isBoundary = true;
                Eigen::RowVector3d rightBarycenter = mesh.barycenters.row(mesh.EF(i, 1));
                dualEdgeLength += (rightBarycenter - midEdge).norm();
            }
            M1Weights(i)=(isBoundary ? 0.5 : 1.0)*dualEdgeLength/primalEdgeLengths;
        } else {//cot weights
            M1Weights(i)=0.0;
            for (int j=0;j<2;j++){
                if (mesh.EF(i,j)==-1)
                    continue;
                int nextHe = mesh.nextH(mesh.EH(i,j));
                int prevHe = mesh.prevH(mesh.EH(i,j));
                
                Eigen::RowVector3d e1 = -(mesh.V.row(mesh.HV(mesh.nextH(nextHe)))-mesh.V.row(mesh.HV(nextHe)));
                Eigen::RowVector3d e2 = mesh.V.row(mesh.HV(mesh.nextH(prevHe)))-mesh.V.row(mesh.HV(prevHe));
                double sinAngle = e1.cross(e2).norm();
                double cosAngle = e1.dot(e2);
                if (std::abs(sinAngle)<10e-7)
                    continue;  //using 0 weight
                
                M1Weights(i)+=0.5*cosAngle/sinAngle;
            }
        }
    }
    
    hodgeStar = directional::sparse_diagonal(M1Weights);
    Eigen::VectorXd invWeights = M1Weights.unaryExpr([](double v) {
        return (std::abs(v)>10e-7 ? 1.0/v : 0.0);
    });
    invHodgeStar = directional::sparse_diagonal(invWeights);
}


}

#endif
