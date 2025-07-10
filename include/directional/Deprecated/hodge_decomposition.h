// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_HODGE_DECOMPOSITION_H
#define DIRECTIONAL_HODGE_DECOMPOSITION_H

#include <queue>
#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <igl/igl_inline.h>
#include <igl/diag.h>
#include <igl/local_basis.h>
#include <igl/edge_topology.h>
#include <directional/FEM_masses.h>
#include <directional/FEM_suite.h>
#include <igl/per_face_normals.h>


namespace directional
{


// Creating non-conforming mid-edge mesh, where the faces are between the midedges of each original face. This is generally only for visualization
// Input:
//  VMesh:      #V x 3 conforming mesh vertices
//  FMesh:      #F x 3 conforming mesh faces
//  EV:         #E x 2 edges to vertices indices
//  EF:         #E x 2 edges to faces indices
//  FE:         #F x 3 faces to edges indices
// Output:
//  Gv:    #3f x V Conforming gradient matrix, returning vector of xyzxyz per face gradient vectors
//  Ge:    #3f x V Non-conforming gradient of the same style, but for mid-edge functions
//  J:    #3F x 3F rotation operator [Nx] per face
//  C:  Curl operator which is basically (JGe)^T * Mchi
//  C:  Divergence operator which is basically Gv^T * Mchi

IGL_INLINE void hodge_decomposition(const Eigen::MatrixXd& V,
                                    const Eigen::MatrixXi& F,
                                    const Eigen::MatrixXi& EV,
                                    const Eigen::MatrixXi& FE,
                                    const Eigen::MatrixXi& EF,
                                    const Eigen::MatrixXd& rawField,
                                    Eigen::VectorXd& exactFunc,
                                    Eigen::VectorXd& coexactFunc,
                                    Eigen::MatrixXd& harmField)
{
    
    using namespace Eigen;
    using namespace std;
    
    SparseMatrix<double> Gv, Ge, J, Mv, Mchi, Mf, Me, C, D;
    VectorXd MvVec, MeVec, MfVec, MchiVec;
    
    VectorXd rawFieldVec(3*F.rows(),1);
    for (int i=0;i<F.rows();i++)
        rawFieldVec.segment(3*i,3)=rawField.row(i).transpose();
    
    directional::FEM_suite(V, F, EV, FE, EF, Gv, Ge, J, C, D);
    directional::FEM_masses(V, F, EV, FE, EF, MvVec, MeVec, MfVec, MchiVec);
    
    igl::diag(MvVec,Mv);
    igl::diag(MeVec,Me);
    igl::diag(MfVec,Mf);
    igl::diag(MchiVec,Mchi);
    
    SparseMatrix<double> Lv = D*Gv;   //Gv^T * Mchi * Gv
    SparseMatrix<double> Le = C*J*Ge; //(JGe)^T * Mchi * JGe
    
    //solving for exact part
    VectorXd scalarFunc;
    igl::min_quad_with_fixed_data<double> mqwfExact;
    VectorXd B = -D*rawFieldVec;
    VectorXd Beq;
    SparseMatrix<double> Aeq;
    Eigen::VectorXi b(1); b(0)=0;
    Eigen::VectorXd bc(1); bc(0)=0;
    
    igl::min_quad_with_fixed_precompute(Lv,b,Aeq,true,mqwfExact);
    igl::min_quad_with_fixed_solve(mqwfExact,B,bc,Beq,exactFunc);
    
    //solving for coexact part
    igl::min_quad_with_fixed_data<double> mqwfCoexact;
    // Linear term is 0
    B = -C*rawFieldVec;
    
    igl::min_quad_with_fixed_precompute(Le,b,Aeq,true,mqwfCoexact);
    igl::min_quad_with_fixed_solve(mqwfCoexact,B,bc,Beq,coexactFunc);
    
    Eigen::VectorXd gradFieldVec = Gv*exactFunc;
    Eigen::VectorXd rotCogradFieldVec = J*Ge*coexactFunc;
    Eigen::VectorXd harmFieldVec = rawFieldVec - gradFieldVec -rotCogradFieldVec;
    
    harmField.resize(F.rows(),3);
    for (int i=0;i<F.rows();i++)
        for (int j=0;j<3;j++)
            harmField(i,j)=harmFieldVec(3*i+j);
}
}

#endif


