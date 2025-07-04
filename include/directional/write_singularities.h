// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_WRITE_SINGULARITIES_H
#define DIRECTIONAL_WRITE_SINGULARITIES_H

#include <cmath>
#include <Eigen/Core>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>


namespace directional
{
/*** Writes a list of singularities to a file.
Input:
fileName:               The file name to which the singularities should be saved.
N:                          The degree of the field
singLocalCycles:   The singular local cycles.
singIndices:            The integer index of the singularities, where the actual fractional index is singIndices/N.
Output:
 Whether or not the file was written successfully
 ***/
bool inline write_singularities(const std::string &fileName,
                                    const int N,
                                    const Eigen::VectorXi &singLocalCycles,
                                    const Eigen::VectorXi &singIndices)
{
    std::ofstream f(fileName, std::ios::trunc);
    f << N << " " << singIndices.size() << std::endl;
    
    for (int i=0;i<singIndices.rows();i++)
        f << singLocalCycles(i) << " " << singIndices(i) << std::endl;
    
    f.close();
    return !f.fail();
}
}

#endif
