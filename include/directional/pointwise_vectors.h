// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POINTWISE_VECTORS_H
#define DIRECTIONAL_POINTWISE_VECTORS_H

#include <Eigen/Core>
#include <vector>
#include <directional/TriMesh.h>
#include <directional/polyvector_to_raw.h>
#include <directional/barycentric_coordinates.h>


namespace directional{

//This function interpolates vectors inside each face from the polynomial dofs


/*IGL_INLINE void multiply_polynomials(Eigen::RowVectorXcd& p1,
                                     const Eigen::RowVectorXcd& p2)
{
  Eigen::RowVectorXcd newp = Eigen::RowVectorXcd::Zero(p1.size()+p2.size());
  for (int i=0;i<p1.size();i++)
    for (int j=0;j<p2.size();j++)
      newp(i+j)+=p1(i)*p2(j);
  
  p1 =newp;
}*/


inline void raw_to_polyvector_polynomial(const Eigen::MatrixXd& V,
                                  const Eigen::MatrixXi& F,
                                  const Eigen::MatrixXd& B1,
                                  const Eigen::MatrixXd& B2,
                                  const int N,
                                  const int K,
                                  const Eigen::MatrixXd& rawField,
                                  Eigen::MatrixXcd& pvField){
  
  using namespace std;
  int dofPerFace = (K+1)*(K+2)/2;
  pvField.conservativeResize(dofPerFace*F.rows(), N);
  for (int i=0;i<dofPerFace*F.rows();i++){
    int currFace = (i - i%dofPerFace)/dofPerFace;
    Eigen::RowVectorXcd complexField(N);
    for (int j=0;j<N;j++){
      Eigen::RowVector3d currVector = rawField.block(i,j*3, 1,3);
      complexField(j) = std::complex<double>(currVector.dot(B1.row(currFace)),currVector.dot(B2.row(currFace)));
    }
    //std::cout<<"complexField: "<<complexField<<endl;
    Eigen::RowVectorXcd pvCoeffs=Eigen::RowVectorXcd::Zero(2);
    pvCoeffs(0)=-complexField(0);
    pvCoeffs(1)=1;
    for (int j=1;j<N;j++){
      Eigen::RowVectorXcd newRoot=Eigen::RowVectorXcd::Zero(2);
      newRoot(0)=-complexField(j);
      newRoot(1)=1;
      multiply_polynomials(pvCoeffs, newRoot);
    }
    pvField.row(i)=pvCoeffs.head(N);
    if (N%2==0)  //force pvField to have sign symmetry, meaning that odd powers are strictly zero
      for (int n=1;n<N;n+=2)
        pvField.col(n).setZero();
    //std::cout<<"pvField.row(i): "<<pvField.row(i)<<endl;
  }
}


    void inline pointwise_vectors(const TriMesh& mesh,
                                    const Eigen::VectorXi& faces,
                                    const Eigen::MatrixXd& locations,
                                    const Eigen::MatrixXd& rawField,
                                    const int order,
                                    Eigen::MatrixXd& interpField)


  {
    using namespace Eigen;
    using namespace std;
    
    //only going up to quadratic vectors for now (integrating to cubic functions)
    assert(order<=2 && order >=0);
    
    int N = rawField.cols()/3;
    //cout<<"vecIndices: "<<vecIndices<<endl;
    //cout<<"locations: "<<locations<<endl;
    int dofInFace = (order+1)*(order+2)/2;
    
    interpField.resize(faces.rows(), rawField.cols());
    
    MatrixXd BCoords;  //barycentric coordinate of locations inside respective triangles
    MatrixXd cornersi(locations.rows(),3), cornersj(locations.rows(),3), cornersk(locations.rows(),3);
    for (int i=0;i<faces.rows();i++){
      cornersi.row(i)=mesh.V.row(mesh.F(faces(i),0));
      cornersj.row(i)=mesh.V.row(mesh.F(faces(i),1));
      cornersk.row(i)=mesh.V.row(mesh.F(faces(i),2));
    }
    
    directional::barycentric_coordinates(locations, cornersi, cornersj, cornersk, BCoords);

    MatrixXcd pvField;
    directional::raw_to_polyvector_polynomial(mesh.V,mesh.F,mesh.FBx,mesh.FBy,N,order, rawField,pvField);
    MatrixXcd interpPVField(faces.rows(),N);

    MatrixXd B1interp(faces.rows(),3);
    MatrixXd B2interp(faces.rows(),3);
    /*isInside.resize(faces.rows());
    for (int i=0;i<faces.rows();i++)
      isInside(i) = (BCoords(i,0)>=0.0)&&(BCoords(i,1)>=0.0)&&(BCoords(i,2)>=0.0)&&(BCoords(i,0)<=1.0)&&(BCoords(i,1)<=1.0)&&(BCoords(i,2)<=1.0);*/
    for (int i=0;i<faces.rows();i++){
      /*if (!isInside(i))
        continue;  //the result is not defined*/
      
      B1interp.row(i)=mesh.FBx.row(faces(i));
      B2interp.row(i)=mesh.FBy.row(faces(i));
      
      if (order==0){  //the field is constant
        interpPVField.row(i) = pvField.row(faces(i));
        continue;
      }
      
      RowVectorXd basisFunctions(dofInFace);
      if (order==1){
        basisFunctions = BCoords.row(i);
      }
      
      if (order==2){
        //order: i,j,k,ij,jk,ki
        basisFunctions(0) = 2.0*BCoords(i,0)*BCoords(i,0) - BCoords(i,0);
        basisFunctions(1) = 2.0*BCoords(i,1)*BCoords(i,1) - BCoords(i,1);
        basisFunctions(2) = 2.0*BCoords(i,2)*BCoords(i,2) - BCoords(i,2);
        basisFunctions(3) = 4.0*BCoords(i,0)*BCoords(i,1);
        basisFunctions(4) = 4.0*BCoords(i,1)*BCoords(i,2);
        basisFunctions(5) = 4.0*BCoords(i,2)*BCoords(i,0);
      }

      interpPVField.row(i) = basisFunctions*pvField.block(dofInFace*faces(i), 0, dofInFace,N);
      //cout<<"pvField.block(dofInFace*faces(i), 0, dofInFace,N): "<<pvField.block(dofInFace*faces(i), 0, dofInFace,N)<<endl;
      //cout<<"interpPVField.row(i): "<<interpPVField.row(i)<<endl;
    }

    MatrixXcd interpComplexField;
    //hack that should be fixed
    MatrixXd realInterpPVField(interpPVField.rows(), interpPVField.cols()*2);
    for (int i=0;i<interpPVField.cols();i+=2){
        realInterpPVField.col(i) = interpPVField.col(i/2).real();
        realInterpPVField.col(i+1) = interpPVField.col(i/2).imag();
    }
    directional::polyvector_to_raw(realInterpPVField, N, interpComplexField);
    //converting to raw field
    //interpField.resize(interpComplexField.rows(), interpComplexField.cols()*3);
    for (int i=0;i<interpComplexField.rows();i++)
        for (int j=0;j<N;j++)
            interpField.block(i,3*j, 1, 3) = B1interp.row(i)*interpComplexField(i,j).real()+B2interp.row(i)*interpComplexField(i,j).imag();
    //cout<<"interpField: "<<interpField<<endl;

    
  }
}



#endif
