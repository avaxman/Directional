// This file is part of libhedra, a library for polyhedral mesh processing
//
// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_VERTEX_AREA_MESH_H
#define DIRECTIONAL_VERTEX_AREA_MESH_H
#include <igl/igl_inline.h>
#include <igl/barycenter.h>
#include <igl/parula.h>
#include <Eigen/Core>
#include <string>
#include <vector>
#include <cstdio>

namespace directional
{
    // returns a triangle mesh s.t. every face is tesselated by a central vertex, and then every edge is supported by two triangles. The purpose is to visualize edge-based quantities)
    // Inputs:
    //  V  eigen double matrix     #V by 3 - vertex coordinates
    //  D  eigen int vector        #F by 1 - face degrees
    //  F  eigen int matrix        #F by max(D) - vertex indices in face
    //  EV eigen int matrix     #E by 2 - map from edges to end vertices
    //  EF eigen int matrix     #E by 2 - map from edges to adjacent faces
    //
    // Outputs:
    //  edgeV  eigen double matrix  #F+#V by 3 - new vertices
    //  edgeT  eigen int matrix    2*#E-#Boundary - new edge-based triangles
    //  edgeTE eigen int vector     #edgeT edgeT -> original edge in EV.
    IGL_INLINE bool vertex_area_mesh(const Eigen::MatrixXd& V,
                                     const Eigen::MatrixXi& F,
                                     const Eigen::MatrixXi& EV,
                                     const Eigen::MatrixXi& FE,
                                     const Eigen::MatrixXi& EF,
                                     const Eigen::VectorXd& vertexAreaFunc,
                                     Eigen::MatrixXd& VVAMesh,
                                     Eigen::MatrixXi& FVAMesh,
                                     Eigen::MatrixXd& CVAMesh)
    {
        using namespace Eigen;
        MatrixXd faceCenters;
        igl::barycenter(V,F,faceCenters);

        MatrixXd midEdges(EV.rows(),3);
        for (int i=0;i<EV.rows();i++)
            midEdges.row(i) = (V.row(EV(i,0))+V.row(EV(i,1)))/2.0;

        VVAMesh.conservativeResize(V.rows()+EV.rows()+F.rows(),3);
        VVAMesh<<V, midEdges, faceCenters;

        FVAMesh.conservativeResize(6*F.rows(),3);
        VectorXd funcScalars(6*F.rows());

        for (int i=0;i<F.rows();i++){
            FVAMesh.block(6*i,0,6,3)<<F(i,0), V.rows()+FE(i,0), V.rows()+EV.rows()+i,
                    V.rows()+EV.rows()+i, V.rows()+FE(i,2), F(i,0),
                    F(i,1), V.rows()+FE(i,1), V.rows()+EV.rows()+i,
                    V.rows()+EV.rows()+i, V.rows()+FE(i,0), F(i,1),
                    F(i,2), V.rows()+FE(i,2), V.rows()+EV.rows()+i,
                    V.rows()+EV.rows()+i, V.rows()+FE(i,1), F(i,2);

            funcScalars.segment(6*i,6)<<vertexAreaFunc(F(i,0)), vertexAreaFunc(F(i,0)),
                    vertexAreaFunc(F(i,1)),vertexAreaFunc(F(i,1)),
                    vertexAreaFunc(F(i,2)),vertexAreaFunc(F(i,2));

        }

        igl::parula(funcScalars, funcScalars.minCoeff(), funcScalars.maxCoeff(), CVAMesh);

        return true;
    }
}


#endif


