// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_ROTATION_TO_RAW_H
#define DIRECTIONAL_ROTATION_TO_RAW_H

#include <directional/definitions.h>
#include <directional/CartesianField.h>

namespace directional
{
    // Converts the rotation angle representation + global rotation to raw format
    // Inputs:
    //  tb              The tangent bundle on which the field should be defined. The field is initialized with this tangent bundle (by reference!)
    //  rotationAngles: #E angles that encode deviation from parallel transport EF(i,0)->EF(i,1)
    //  N:              The degree of the field.
    //  globalRotation: The angle between the vector on the first face and its basis in radians.
    // Outputs:
    //  field:          The raw Cartesian field.
    inline void rotation_to_raw(const TangentBundle& tb,
                                const Eigen::VectorXd& rotationAngles,
                                const int N,
                                const double globalRotation,
                                directional::CartesianField& field)
    {
        typedef std::complex<double> Complex;
        using namespace Eigen;
        using namespace std;

        Complex globalRot = exp(Complex(0, globalRotation));
        field.init(tb, fieldTypeEnum::RAW_FIELD, N);
        SparseMatrix<Complex> aP1Full(tb.adjSpaces.rows(), field.intField.rows());
        SparseMatrix<Complex> aP1(tb.adjSpaces.rows(), field.intField.rows() - 1);
        vector<Triplet<Complex> > aP1Triplets, aP1FullTriplets;
        for (int i = 0; i<tb.adjSpaces.rows(); i++) {
            if (tb.adjSpaces(i, 0) == -1 || tb.adjSpaces(i, 1) == -1)
                continue;

            aP1FullTriplets.push_back(Triplet<Complex>(i, tb.adjSpaces(i, 0), pow(tb.connection(i),(double)N)*exp(Complex(0, (double)N*rotationAngles(i)))));
            aP1FullTriplets.push_back(Triplet<Complex>(i, tb.adjSpaces(i, 1), -1.0));
            if (tb.adjSpaces(i, 0) != 0)
                aP1Triplets.push_back(Triplet<Complex>(i, tb.adjSpaces(i, 0)-1, pow(tb.connection(i),(double)N)*exp(Complex(0, (double)N*rotationAngles(i)))));
            if (tb.adjSpaces(i, 1) != 0)
                aP1Triplets.push_back(Triplet<Complex>(i, tb.adjSpaces(i, 1)-1, -1.0));
        }
        aP1Full.setFromTriplets(aP1FullTriplets.begin(), aP1FullTriplets.end());
        aP1.setFromTriplets(aP1Triplets.begin(), aP1Triplets.end());
        VectorXcd torhs = VectorXcd::Zero(field.intField.rows()); torhs(0) = globalRot;  //global rotation
        VectorXcd rhs = -aP1Full*torhs;

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<Complex> > solver;
        solver.compute(aP1.adjoint()*aP1);
        assert(solver.info() == Success);
        VectorXcd complexPowerField(field.intField.rows());
        complexPowerField(0) = globalRot;
        complexPowerField.tail(field.intField.rows() - 1) = solver.solve(aP1.adjoint()*rhs);
        assert(solver.info() == Success);

        VectorXcd complexField = pow(complexPowerField.array(), 1.0 / (double)N);

        MatrixXd intField(complexField.rows(),2*N);
        for (int j=0;j<N;j++){
            VectorXcd currComplexField = complexField.array()*exp(Complex(0,2*j*directional::PI/N));
            intField.block(0,2*j,intField.rows(),2)<<currComplexField.array().real(),currComplexField.array().imag();
        }
        //constructing raw intField
        field.set_intrinsic_field(intField);
    }
}

#endif
