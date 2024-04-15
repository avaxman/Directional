// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_FEM_SUITE_H
#define DIRECTIONAL_FEM_SUITE_H

#include <queue>
#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <directional/FEM_masses.h>


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
  
  inline void FEM_suite(const TriMesh& mesh,
                            Eigen::SparseMatrix<double>& Gv,
                            Eigen::SparseMatrix<double>& Ge,
                            Eigen::SparseMatrix<double>& J,
                            Eigen::SparseMatrix<double>& C,
                            Eigen::SparseMatrix<double>& D)
  {
    
    using namespace Eigen;
    using namespace std;
    
    VectorXd dblA = mesh.faceAreas*2.0;
    SparseMatrix<double> Mchi;
    VectorXd MvVec, MeVec, MfVec, MchiVec;
    directional::FEM_masses(mesh, MvVec, MeVec, MfVec, MchiVec);
    Mchi = MchiVec.asDiagonal();
    Eigen::MatrixXd N = mesh.faceNormals;
    //igl::per_face_normals(V, F, N);
    
    vector<Triplet<double> > GvTriplets, GeTriplets, JTriplets;
    for (int i=0;i<mesh.F.rows();i++){
      RowVector3d currNormal=N.row(i);
      for (int j=0;j<3;j++){
        RowVector3d eVec = mesh.V.row(mesh.F(i,(j+1)%3))-mesh.V.row(mesh.F(i,j));
        RowVector3d eVecRot = currNormal.cross(eVec);
        
        //TODO: I cannot count on FE(i,j) to be the correct edge, need to search for it
        int currEdge=-1;
        for (int k=0;k<3;k++){
          if (((mesh.F(i,j) == mesh.EV(mesh.FE(i,k),0))&&(mesh.F(i,(j+1)%3) == mesh.EV(mesh.FE(i,k),1)))||
            ((mesh.F(i,(j+1)%3) == mesh.EV(mesh.FE(i,k),0))&&(mesh.F(i,j) == mesh.EV(mesh.FE(i,k),1))))
            currEdge=mesh.FE(i,k);
        }
        assert (currEdge!=-1 && "Something wrong with edge topology!");
        for (int k=0;k<3;k++){
          GvTriplets.push_back(Triplet<double>(3*i+k,mesh.F(i,(j+2)%3),eVecRot(k)/dblA(i)));
          GeTriplets.push_back(Triplet<double>(3*i+k,currEdge,-2*eVecRot(k)/dblA(i)));
        }
      }
      
      JTriplets.push_back(Triplet<double>(3*i, 3*i+1, -N(i,2)));
      JTriplets.push_back(Triplet<double>(3*i+1, 3*i, N(i,2)));
      JTriplets.push_back(Triplet<double>(3*i, 3*i+2, N(i,1)));
      JTriplets.push_back(Triplet<double>(3*i+2, 3*i, -N(i,1)));
      JTriplets.push_back(Triplet<double>(3*i+1, 3*i+2, -N(i,0)));
      JTriplets.push_back(Triplet<double>(3*i+2, 3*i+1, N(i,0)));
    }
    
    Gv.resize(3*mesh.F.rows(), mesh.V.rows());
    Gv.setFromTriplets(GvTriplets.begin(), GvTriplets.end());
    Ge.resize(3*mesh.F.rows(), mesh.EV.rows());
    Ge.setFromTriplets(GeTriplets.begin(), GeTriplets.end());
    J.resize(3*mesh.F.rows(), 3*mesh.F.rows());
    J.setFromTriplets(JTriplets.begin(), JTriplets.end());
    
    C = (J*Ge).transpose()*Mchi;
    D = Gv.transpose()*Mchi;
  }
}


#endif


