// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_READDMAT_H
#define DIRECTIONAL_READDMAT_H

#include <cmath>
#include <Eigen/Core>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>

namespace directional{

    //Reading DMAT files that contain matrices
    template<typename T>
    bool inline readDMAT(const std::string dmatFileName,
                         Eigen::PlainObjectBase<T>& dmat) {

        std::ifstream file(dmatFileName);
        if (!file.is_open())
            return false;

        int rows, cols;
        file >> cols >> rows;  //DMAT is column-major
        dmat.resize(rows, cols);

        for (int i = 0; i < cols; i++)
            for (int j = 0; j < rows; j++)
                file >> dmat(j, i);

        file.close();
        return true;
    }
}

#endif 
