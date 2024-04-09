// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_READ_RAW_FIELD_H
#define DIRECTIONAL_READ_RAW_FIELD_H
#include <cmath>
#include <Eigen/Core>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include <directional/TangentBundle.h>
#include <directional/CartesianField.h>

namespace directional
{

    // Reads a raw *extrinsic* cartesian field from a file and initializes a Cartesian file object, including projecting to the intrinsic tangent spaces
    // Inputs:
    //   fileName: The to be loaded file.
    //    tb: the underlying tangent bundle to the read field
    // Outputs:
    //   N: The degree of the field
    //   field: the read RAW_FIELD field.
    // Return:
    //   Whether or not the file was read successfully
    bool IGL_INLINE read_raw_field(const std::string &fileName,
                                   const directional::TangentBundle& tb,
                                   int& N,
                                   directional::CartesianField& field)
    {
        try
        {
            std::ifstream f(fileName);
            if (!f.is_open()) {
                return false;
            }
            int numT;
            f>>N;
            f>>numT;
            Eigen::MatrixXd extField;
            extField.conservativeResize(numT, 3*N);

            //Can we do better than element-wise reading?
            for (int i=0;i<extField.rows();i++)
                for (int j=0;j<extField.cols();j++)
                    f>>extField(i,j);

            f.close();
            assert(tb.sources.rows()==extField.rows());
            assert(tb.hasEmbedding() && "This tangent bundle doesn't admit an extrinsic embedding");
            field.init(tb, fieldTypeEnum::RAW_FIELD, N);
            field.set_extrinsic_field(extField);
            return f.good();
        }
        catch (std::exception e)
        {
            return false;
        }
    }
}

#endif
