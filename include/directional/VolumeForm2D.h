// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_VOLUME_FUNCTION_H
#define DIRECTIONAL_VOLUME_FUNCTION_H

#include <eigen/sparse>

namespace directional{

    //an interface to a single element in the class
    template<typename NumberType>
    class VolumeForm2D{
    public:

        VolumeForm2D(){}
        ~VolumeForm2D(){}

        Eigen::SparseMatrix<NumberType> virtual mass_matrix() = 0;
        Eigen::SparseMatrix<NumberType> virtual inv_mass_matrix() = 0;
    };

}


#endif
