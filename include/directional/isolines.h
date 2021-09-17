// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Oded Stein <oded.stein@columbia.edu>, 2021 Amir Vaxman <avaxman@gmail.com> (modifications)
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


#ifndef IGL_ISOLINES_H
#define IGL_ISOLINES_H
#include <igl/igl_inline.h>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <vector>
#include <igl/remove_duplicate_vertices.h>
#include <igl/per_face_normals.h>


namespace igl
  {
  // Constructs isolines for a function z given on a mesh (V,F)
  //
  //
  // Inputs:
  //   V  #V by dim list of mesh vertex positions
  //   F  #F by 3 list of mesh faces (must be triangles)
  //   z  #V by 1 list of function values evaluated at vertices
  //   n  the number of desired isolines
  // Outputs:
  //   isoV  #isoV by dim list of isoline vertex positions
  //   isoE  #isoE by 2 list of isoline edge positions
  //   isoN  $isoE by 3 list of normals to the edges
  //
  
  template <typename DerivedV,
  typename DerivedF,
  typename DerivedZ,
  typename DerivedIsoV,
  typename DerivedIsoE>
  IGL_INLINE void isolines(
                           const Eigen::MatrixBase<DerivedV>& V,
                           const Eigen::MatrixBase<DerivedF>& F,
                           const Eigen::MatrixBase<DerivedZ>& z,
                           const int n,
                           Eigen::PlainObjectBase<DerivedIsoV>& isoV,
                           Eigen::PlainObjectBase<DerivedIsoE>& isoE,
                           Eigen::PlainObjectBase<DerivedIsoV>& isoN){
    //Constants
    const int dim = V.cols();
    assert(dim==2 || dim==3);
    //const int nVerts = V.rows();
    assert(z.rows() == nVerts &&
           "There must be as many function entries as vertices");
    const int nFaces = F.rows();
    int np1 = n+1;
    const double min = z.minCoeff(), max = z.maxCoeff();
    
    
    //Following http://www.alecjacobson.com/weblog/?p=2529
    typedef typename DerivedZ::Scalar Scalar;
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vec;
    
    //Amir's corrections to integers
    int roundmin = floor(min);
    int roundmax = ceil(max);
    np1 = roundmax - roundmin + 1;
    Vec iso(np1);
    
    
    /*std::cout<<"roundmin: "<<roundmin<<std::endl;
     std::cout<<"roundmax: "<<roundmax<<std::endl;*/
    
    for (int i=roundmin;i<=roundmax;i++)
      iso(i-roundmin)=i;
    
    
    /*for(int i=0; i<np1; ++i)
     iso(i) = Scalar(i)/Scalar(np1-1)*(roundmax-roundmin) + roundmin;*/
    
    //std::cout<<"iso: "<<iso<<std::endl;
    
    Eigen::Matrix<Scalar, Eigen::Dynamic, 3> normals;
    igl::per_face_normals(V, F, normals);
    
    typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
    std::array<Matrix,3> t{{Matrix(nFaces, np1),
      Matrix(nFaces, np1), Matrix(nFaces, np1)}};
    for(int i=0; i<nFaces; ++i) {
      for(int k=0; k<3; ++k) {
        const Scalar z1=z(F(i,k)), z2=z(F(i,(k+1)%3));
        for(int j=0; j<np1; ++j) {
          t[k](i,j) = (iso(j)-z1) / (z2-z1);
          if(t[k](i,j)<0 || t[k](i,j)>1)
            t[k](i,j) = std::numeric_limits<Scalar>::quiet_NaN();
        }
      }
    }
    
    
    std::array<std::vector<int>,3> Fij, Iij;
    for(int i=0; i<nFaces; ++i) {
      for(int j=0; j<np1; ++j) {
        for(int k=0; k<3; ++k) {
          const int kp1=(k+1)%3, kp2=(k+2)%3;
          if(std::isfinite(t[kp1](i,j)) && std::isfinite(t[kp2](i,j))) {
            Fij[k].push_back(i);
            Iij[k].push_back(j);
          }
        }
      }
    }
    
    const int K = Fij[0].size()+Fij[1].size()+Fij[2].size();
    isoV.resize(2*K, dim);
    if (dim==3)
      isoN.resize(K,3);
    
    int b = 0;
    for(int k=0; k<3; ++k) {
      const int kp1=(k+1)%3, kp2=(k+2)%3;
      for(int i=0; i<Fij[k].size(); ++i) {
        isoV.row(b+i) = (1.-t[kp1](Fij[k][i],Iij[k][i]))*
        V.row(F(Fij[k][i],kp1)) +
        t[kp1](Fij[k][i],Iij[k][i])*V.row(F(Fij[k][i],kp2));
        isoV.row(K+b+i) = (1.-t[kp2](Fij[k][i],Iij[k][i]))*
        V.row(F(Fij[k][i],kp2)) +
        t[kp2](Fij[k][i],Iij[k][i])*V.row(F(Fij[k][i],k));
        if (dim==3)
          isoN.row(b+i)=normals.row(Fij[k][i]);
      }
      b += Fij[k].size();
    }
    
    isoE.resize(K,2);
    for(int i=0; i<K; ++i)
      isoE.row(i) << i, K+i;
    
    
    //Remove double entries
    typedef typename DerivedIsoV::Scalar LScalar;
    typedef typename DerivedIsoE::Scalar LInt;
    typedef Eigen::Matrix<LInt, Eigen::Dynamic, 1> LIVec;
    typedef Eigen::Matrix<LScalar, Eigen::Dynamic, Eigen::Dynamic> LMat;
    typedef Eigen::Matrix<LInt, Eigen::Dynamic, Eigen::Dynamic> LIMat;
    LIVec dummy1, dummy2;
    igl::remove_duplicate_vertices(LMat(isoV), LIMat(isoE),
                                   2.2204e-15, isoV, dummy1, dummy2, isoE);
    
  }
  }


#endif
