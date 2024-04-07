//
// Created by Amir Vaxman on 07.04.24.
//

#ifndef DIRECTIONAL_READDMAT_H
#define DIRECTIONAL_READDMAT_H

#include <cmath>
#include <Eigen/Core>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>

namespace directional{

    //A wrapper around libigl readOBJ that uses the mesh class
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

#endif //DIRECTIONAL_READDMAT_H
