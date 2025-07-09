// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POLY_TO_CONSTANT_SUBDIVISION_H
#define DIRECTIONAL_POLY_TO_CONSTANT_SUBDIVISION_H

#include <Eigen/Core>
#include <vector>
#include <directional/pointwise_vectors.h>
#include <directional/TriMesh.h>


namespace directional{



//subdividing a mesh with a polynomial field (of order K) into a mesh with a piecewise constant sampled field. Used for visualization alone!
//subdividing 4-1 subdDepth
void inline poly_to_constant_subdivision(const TriMesh& mesh,
                                         const Eigen::MatrixXd& rawField,
                                         const int K,
                                         const int subdDepth,
                                         TriMesh& subdMesh,
                                         Eigen::VectorXi& FFK,
                                         Eigen::MatrixXd& rawFieldK)
{
    using namespace Eigen;
    
    subdMesh = mesh;
    
    int degree = rawField.cols()/3;
    //iteratively doing 4-1 subdivisions
    for (int i=0;i<subdDepth;i++){
        MatrixXd subdVK;
        MatrixXi subdFK;
        subdVK.resize(subdMesh.V.rows()+subdMesh.EV.rows(),3);
        MatrixXd midEdges(subdMesh.EV.rows(),3);
        for (int j=0;j<subdMesh.EV.rows();j++)
            midEdges.row(j)<<(subdMesh.V.row(subdMesh.EV(j,0))+subdMesh.V.row(subdMesh.EV(j,1)))/2.0;
        subdVK<<subdMesh.V, midEdges;
        subdFK.resize(subdMesh.F.rows()*4,3);
        for (int j=0;j<subdMesh.F.rows();j++){
            subdFK.row(4*j)<<subdMesh.F(j,0), subdMesh.V.rows()+subdMesh.FE(j,0), subdMesh.V.rows()+subdMesh.FE(j,2);
            subdFK.row(4*j+1)<<subdMesh.F(j,1), subdMesh.V.rows()+subdMesh.FE(j,1), subdMesh.V.rows()+subdMesh.FE(j,0);
            subdFK.row(4*j+2)<<subdMesh.F(j,2), subdMesh.V.rows()+subdMesh.FE(j,2), subdMesh.V.rows()+subdMesh.FE(j,1);
            subdFK.row(4*j+3)<<subdMesh.V.rows()+subdMesh.FE(j,0), subdMesh.V.rows()+subdMesh.FE(j,1), subdMesh.V.rows()+subdMesh.FE(j,2);
        }
        
        TriMesh newMesh;
        newMesh.set_mesh(subdVK, subdFK);
        
        subdMesh = newMesh;
        
        /*FK = subdFK;
         VK = subdVK;
         
         igl::edge_topology(VK, FK, EVK, FEK, EFK);*/
    }
    
    int numSubdFaces = pow(4,subdDepth);
    FFK.resize(mesh.F.rows()*numSubdFaces);
    for (int i=0;i<mesh.F.rows();i++)
        FFK.segment(i*numSubdFaces, numSubdFaces).array()=i;
    
    /*Eigen::MatrixXd barycentersK(FFK.rows(),3);
     for (int i=0;i<FFK.rows();i++)
     barycentersK.row(i) = mesh.barycenters.row(FFK(i));*/
    
    directional::pointwise_vectors(mesh, FFK, subdMesh.barycenters, rawField,K,rawFieldK);
    
}




}


#endif /* poly_to_constant_subdivision_h */
