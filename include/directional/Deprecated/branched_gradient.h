// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


#ifndef DIRECTIONAL_BRANCHED_GRADIENT_H
#define DIRECTIONAL_BRANCHED_GRADIENT_H


namespace directional{


//assuming corner functions are arranged packets of N per vertex.
inline void branched_gradient(const TriMesh& mesh,
                              const int N,
                              Eigen::SparseMatrix<double>& G)
{
    
    using namespace Eigen;
    using namespace std;
    
    VectorXd dblA = mesh.faceAreas*2.0;
    Eigen::MatrixXd normals = mesh.faceNormals;
    vector<Triplet<double>> GTri;
    
    for (int i=0;i<mesh.F.rows();i++){
        RowVector3d currNormal=normals.row(i);
        for (int k=0;k<N;k++){
            RowVector3d localGradient(0.0,0.0,0.0);
            for (int j=0;j<3;j++){
                RowVector3d eVec = mesh.V.row(mesh.F(i,(j+2)%3))-mesh.V.row(mesh.F(i,(j+1)%3));
                RowVector3d gradComp = currNormal.cross(eVec)/dblA(i);
                for (int l=0;l<3;l++)
                    GTri.push_back(Triplet<double>(3*N*i+k*3+l, N*mesh.F(i,j)+k,gradComp(l)));
                
            }
        }
    }
    G.conservativeResize(3*N*mesh.F.rows(), N*mesh.V.rows());
    G.setFromTriplets(GTri.begin(), GTri.end());
}
}




#endif /* branched_gradient_h */
