// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


#ifndef DIRECTIONAL_RAW_TO_POLYVECTOR_H
#define DIRECTIONAL_RAW_TO_POLYVECTOR_H

#include <directional/CartesianField.h>
#include <directional/polyvector_to_raw.h>

namespace directional{

    //Conversion from raw to polyvector representation on a Cartesian field. This operator works directly on the intrinsic representations.
    //This computes the polynomial coefficients for the given roots.
    //Input:
    //  intField:       #TangentSpaces x 2N (x,y,x,y...) representation of the raw intrinsic field
    //  N:              Degree of the field.
    //  signSymmetry:   Whether the field has sign symmetry or not, only if N is even. Default: true
    //Output:
    //  pvField:        #TangentSpaces x N complex representation of the PolyVector.
    inline void raw_to_polyvector(const Eigen::MatrixXd& intField,
                                      const int N,
                                      Eigen::MatrixXcd& pvField,
                                      const bool signSymmetry=true){


        pvField=Eigen::MatrixXcd::Zero(intField.rows(),intField.cols()/2);
        Eigen::MatrixXcd actualRoots;
        int actualN;
        if ((signSymmetry)&&(N%2==0)){
            actualRoots.resize(intField.rows(), intField.cols()/4);
            actualN = N/2;
            for (int i=0;i<actualN;i++)
                for (int j=0;j<intField.rows();j++)
                    actualRoots(j,i)=std::complex<double>(intField(j,2*i), intField(j,2*i+1));

            actualRoots = actualRoots.array().square();
        }else {
            actualRoots.resize(intField.rows(), intField.cols()/2);
            actualN = N;
            for (int i=0;i<actualN;i++)
                for (int j=0;j<intField.rows();j++)
                    actualRoots(j,i)=std::complex<double>(intField(j,2*i), intField(j,2*i+1));
        }

        int jump = ((signSymmetry)&&(N%2==0) ? 2 : 1);
        Eigen::MatrixXcd actualPVField(pvField.rows(), 1);
        actualPVField.col(0)=-actualRoots.col(0);
        for (int i=1;i<actualN;i++)
            multiply_polynomials(actualPVField, actualRoots.col(i));

        for (int i=0;i<N;i+=jump)
            pvField.col(i)=actualPVField.col(i/jump);
    }


    inline void raw_to_polyvector(const CartesianField& rawField,
                                  CartesianField& pvField,
                                  const bool signSymmetry=true){

        Eigen::MatrixXcd pvFieldComplex;
        raw_to_polyvector(rawField.intField, rawField.N, pvFieldComplex, signSymmetry);
        pvField.set_intrinsic_field(pvFieldComplex);
    }

}

#endif //DIRECTIONAL_RAW_TO_POLYVECTOR_H
