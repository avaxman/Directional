// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POWER_FIELD_H
#define DIRECTIONAL_POWER_FIELD_H

#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <directional/polyvector_field.h>
#include <directional/TangentBundle.h>
#include <directional/CartesianField.h>


namespace directional
{
    // Computes a power field on the entire mesh from given values at the prescribed indices.
    // If no constraints are given the lowest-eigenvalue (of smoothness energy) field will be returned.
    // Input:
    //  tb: underlying tangent bundle.
    //  constFaces: the faces on which the polyvector is prescribed. If a face is repeated and the alignment is hard then all but the first vector in the face will be ignored.
    //  constVectors: #F by 3 in representative form of the N-RoSy's on the tangent spaces.
    //  alignWeights: #constFaces x 1 soft weights for alignment (negative values = fixed faces).
    //  N: The degree of the field.
    // Output:
    //  powerField: a cartesian power-field object.
    inline void power_field(const TangentBundle& tb,
                            const Eigen::VectorXi& constSpaces,
                            const Eigen::MatrixXd& constVectors,
                            const Eigen::VectorXd& alignWeights,
                            const int N,
                            directional::CartesianField& field,
                            const bool normalizeField = false)
    {
        
        PolyVectorData pvData;
        pvData.N = N;
        pvData.tb = &tb;
        if (constSpaces.size()!=0) {
            pvData.constSpaces = constSpaces;
            pvData.constVectors = constVectors;
            pvData.wAlignment = alignWeights;
        }else{
            pvData.constSpaces.resize(1); pvData.constSpaces(0)=0;
            Eigen::RowVector2d intConstVector; intConstVector<<1.0,0.0;
            pvData.constVectors = tb.project_to_extrinsic(pvData.constSpaces, intConstVector);
            pvData.wAlignment = Eigen::VectorXd::Constant(pvData.constSpaces.size(),-1.0);
        }
        pvData.wSmooth = 1.0;
        pvData.wRoSy = -1.0;  //Perfect RoSy
        polyvector_field(pvData, field);
        field.fieldType = fieldTypeEnum::POWER_FIELD;
        //getting rid of the redundant zeros, in case they were allocated.
        field.intField.conservativeResize(field.intField.rows(),2);
        field.extField.conservativeResize(field.extField.rows(),3);
        if (normalizeField){
            Eigen::MatrixXd intField = field.intField;
            intField = intField.array().colwise() / intField.rowwise().norm().array();
            field.set_intrinsic_field(intField);
        }
        //powerField=-pvField.col(0);  //powerfield is represented positively
    }
}


#endif
