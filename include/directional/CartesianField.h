// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_CARTESIAN_FIELD_H
#define DIRECTIONAL_CARTESIAN_FIELD_H

#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <directional/TriMesh.h>

namespace directional{

#define RAW_FIELD 0
#define POWER_FIELD 1
#define POLYVECTOR_FIELD 2


//This is the interface class for any directional fields represented in cartesian coordinates, of any order N.
template<class _TangentBundle>
class CartesianField{
public:
  
  const _TangentBundle* tb;
  
  int N;  //degree of field (how many vectors are in each point);
  int fieldType;       //What the field actually represents  (for instance, either a raw field or a power/polyvector field)
  
  Eigen::MatrixXd intField;  //the field in intrinsic representation (depending on the local basis of the face). Size #T x 2N
  Eigen::MatrixXd extField;  //the field in ambient coordinates. Size Size #T x 3N
  
  Eigen::VectorXi matching;         //matching(i)=j when vector k in adjSpaces(i,0) matches to vector (k+j)%N in adjSpaces(i,1)
  Eigen::VectorXd effort;           //The effort of the entire matching (sum of deviations from parallel transport)
  Eigen::VectorXi singLocalCycles;     //The singular (dual elements). Only the local cycles! not the generators or boundary cycles
  Eigen::VectorXi singIndices;      //Corresponding indices (this is the numerator where the true fractional index is singIndices/N);
  
  CartesianField(){}
  CartesianField(const _TangentBundle& _tb):tb(&_tb){}
  ~CartesianField(){}
  
  //Initializing the field with the proper tangent spaces
  void IGL_INLINE init(const _TangentBundle& _tb, const int _fieldType, const int _N){
    tb = &_tb;
    fieldType = _fieldType;
    N=_N;
  };
  
  void IGL_INLINE set_intrinsic_field(const Eigen::MatrixXd& _intField){
    assert (!(fieldType==POWER_FIELD) || (_intField.cols()==2));
    assert ((_intField.cols()==2*N) || !(fieldType==POLYVECTOR_FIELD || fieldType==RAW_FIELD));
    intField = _intField;
    
    extField = tb->project_to_extrinsic(Eigen::VectorXi(), intField);
  }
  
  //The same, just with complex coordinates
  void virtual IGL_INLINE set_intrinsic_field(const Eigen::MatrixXcd& _intField){
    intField.resize(_intField.rows(),_intField.cols()*2);
    for (int i=0;i<N;i++){
      intField.col(2*i)=_intField.col(i).real();
      intField.col(2*i+1)=_intField.col(i).imag();
    }
    set_intrinsic_field(intField);
  }
  
  //Setting the field by the extrinsic ambient field, which will get projected to the intrinsic tangent spaces.
  void IGL_INLINE set_extrinsic_field(const Eigen::MatrixXd& _extField){
    assert(_extField.cols()==3*N);
    extField=_extField;
    extField = tb->project_to_intrinsic(Eigen::VectorXi(), extField);
  }
  
  
  //Directly setting the singularities of the the field (only at the local dual elements; not at generator or boundary cycles).
  void IGL_INLINE set_singularities(const Eigen::VectorXi& _singLocalCycles,
                                    const Eigen::VectorXi& _singIndices){
    
    //TODO: remove boundary elements
    singLocalCycles = _singLocalCycles;
    singIndices = _singIndices;
  }
};

}



#endif 
