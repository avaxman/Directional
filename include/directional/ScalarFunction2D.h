// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_SCALAR_FUNCTION_H
#define DIRECTIONAL_SCALAR_FUNCTION_H

#include <eigen/sparse>
#include <directional/TriMesh.h>
#include <directional/CartesianField.h>

namespace directional{

    //an interface to a single element in the class
    template<typename NumberType>
    class ScalarFunction2D{
    public:

        ScalarFunction2D(){}
        ~ScalarFunction2D(){}

        //Only good for triangle meshes right now
        NumberType virtual value(const int faceIndex,
                         const Eigen::VectorXd& baryCoords) = 0;

        //TODO: save the matrix for future use
        Eigen::SparseMatrix<NumberType> virtual gradient_matrix() = 0;

        void virtual gradient(directional::CartesianField& gradField) = 0;

        Eigen::SparseMatrix<NumberType> virtual mass_matrix() = 0;
        Eigen::SparseMatrix<NumberType> virtual inv_mass_matrix() = 0;
    };

}


#endif
