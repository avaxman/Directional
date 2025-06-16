// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2021 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POLYVECTOR_ITERATION_FUNCTIONS_H
#define DIRECTIONAL_POLYVECTOR_ITERATION_FUNCTIONS_H

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/Polynomials>
#include <iostream>
#include <directional/TangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/polyvector_to_raw.h>
#include <directional/raw_to_polyvector.h>
#include <directional/project_curl.h>
#include <directional/principal_matching.h>
#include <directional/sparse_diagonal.h>
#include <directional/PolyVectorData.h>


//Normalizes the field to have unit length on all its vectors
CartesianField& hard_normalization(const CartesianField& pvField, const& PolyVectorData& pvData){
    CartesianField rawField, normPvField;
    polyvector_to_raw(pvField, rawField, pvData.N%2==0, true);
    directional::raw_to_polyvector(rawField, normPvField);
    return normPvField;
}

CartesianField& soft_orthonormalization(const CartesianField& pvField, , const& PolyVectorData& pvData){
    
}


#endif
