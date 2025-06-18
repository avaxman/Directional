// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2025 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POLYVECTOR_ITERATION_FUNCTIONS_H
#define DIRECTIONAL_POLYVECTOR_ITERATION_FUNCTIONS_H

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <iostream>
#include <directional/TangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/polyvector_to_raw.h>
#include <directional/raw_to_polyvector.h>
#include <directional/project_curl.h>
#include <directional/principal_matching.h>
#include <directional/sparse_diagonal.h>
#include <directional/PolyVectorData.h>

namespace directional{

typedef std::function<CartesianField(const CartesianField&, const PolyVectorData&)> PvIterationFunction;


//Normalizes the field to have unit length on all its vectors
CartesianField hard_normalization(const CartesianField& pvField, const PolyVectorData& pvData){
    CartesianField rawField, normPvField;
    polyvector_to_raw(pvField, rawField, pvData.N%2==0, true);
    directional::raw_to_polyvector(rawField, normPvField);
    return normPvField;
}


CartesianField hard_rosy(const CartesianField& pvField, const PolyVectorData& pvData){
    Eigen::MatrixXcd rosyField = pvField.get_complex_intrinsic_field();
    rosyField.block(0,1,rosyField.rows(), rosyField.cols()-1).setZero();
    rosyField.col(0) = rosyField.col(0).array() / rosyField.col(0).array().abs();
    //std::cout<<"rosyField: "<<rosyField<<std::endl;
    CartesianField hardRosyField = pvField;
    hardRosyField.set_intrinsic_field(rosyField);
    return hardRosyField;
}



CartesianField soft_rosy(const CartesianField& pvField, const PolyVectorData& pvData){
    Eigen::MatrixXcd rosyField = pvField.get_complex_intrinsic_field();
    rosyField.block(0,1,rosyField.rows(), rosyField.cols()-1).setZero();
    rosyField.col(0) = rosyField.col(0).array() / rosyField.col(0).array().abs();
    //doing closed-form implicit step
    CartesianField softRosyField = pvField;
    Eigen::MatrixXcd interpField = (pvField.get_complex_intrinsic_field().array() + 2.0*pvData.currImplicitCoeff*rosyField.array())/(1+2.0*pvData.currImplicitCoeff);
    softRosyField.set_intrinsic_field(interpField);
    return softRosyField;
}

CartesianField curl_projection(const CartesianField& pvField, const PolyVectorData& pvData){
    assert(pvData.tb->discTangType()==directional::discTangTypeEnum::FACE_SPACES && "Projecting curl only works for face-based fields for now!");
    CartesianField rawField, curlFreeFieldRaw, curlFreeFieldPv;
    polyvector_to_raw(pvField, rawField, pvData.N%2==0, false);
    directional::principal_matching(rawField);
    project_curl(rawField, Eigen::VectorXi(), Eigen::MatrixXd(), curlFreeFieldRaw);
    directional::raw_to_polyvector(curlFreeFieldRaw,  curlFreeFieldPv);
    return curlFreeFieldPv;
}

}


#endif
