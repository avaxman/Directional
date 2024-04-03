// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_READ_SINGULARITIES_H
#define DIRECTIONAL_READ_SINGULARITIES_H
#include <cmath>
#include <Eigen/Core>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include <directional/CartesianField.h>


namespace directional
{

    // Reads a list of element singularities from a file. The identity of the singularities (face or vertex) depend on the file type given.
    // Input:
    //   fileName: The to be loaded file.
    // Output:
    //   N:             The degree of the singularities (so true fractional index is singIndices/N).
    //   singElements:  The singular elements.
    //   singIndices:   The numerator of the index of the singularities.
    // Return:
    //   Whether or not the file was written successfully
    bool inline read_singularities(const std::string &fileName,
                                       int& N,
                                       Eigen::VectorXi& singElements,
                                       Eigen::VectorXi& singIndices)
    {
        try
        {
            std::ifstream f(fileName);
            int numSings;
            f >> N;
            f >> numSings;

            singElements = Eigen::VectorXi::Zero(numSings);
            singIndices = Eigen::VectorXi::Zero(numSings);

            for (int i=0;i<numSings;i++)
                f >> singElements.coeffRef(i)>> singIndices.coeffRef(i);

            f.close();
            return f.fail();
        }
        catch (std::exception e)
        {
            return false;
        }
    }


    //This version reads directly into a field object.
    bool inline read_singularities(const std::string &fileName,
                                       directional::CartesianField& field)
    {
        try
        {
            std::ifstream f(fileName);
            int numSings,N;
            f >> N;
            assert(N==field.N && "Read singularities should be of the same degree as the field");
            f >> numSings;

            Eigen::VectorXi singElements = Eigen::VectorXi::Zero(numSings);
            Eigen::VectorXi singIndices = Eigen::VectorXi::Zero(numSings);

            for (int i=0;i<numSings;i++)
                f >> singElements.coeffRef(i)>> singIndices.coeffRef(i);

            f.close();
            field.set_singularities(singElements, singIndices);
            return f.fail();
        }
        catch (std::exception e)
        {
            return false;
        }
    }
}

#endif
