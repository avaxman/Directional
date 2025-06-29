// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2024 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_WRITE_OBJ_H
#define DIRECTIONAL_WRITE_OBJ_H


#include <Eigen/Core>
#include <string>
#include <vector>
#include <directional/TriMesh.h>

namespace directional
{

bool writeOBJ(const std::string& fileName,
              const directional::TriMesh& mesh,
              const Eigen::MatrixXd& TC,
              const Eigen::MatrixXi& FTC,
              const std::string& mtlFileName = "",
              const std::string& textureName = "")
{
    std::ofstream out(fileName);
    if (!out.is_open())
        throw std::runtime_error("Failed to open file: " + fileName);

    // Write reference to MTL file if provided
    if (!mtlFileName.empty())
        out << "mtllib " << mtlFileName << "\n";

    // Use material if name is provided
    if (!textureName.empty())
        out << "usemtl " << textureName << "\n";

    // Write vertices
    for (int i = 0; i < mesh.V.rows(); ++i)
        out << "v " << mesh.V(i, 0) << " " << mesh.V(i, 1) << " " << mesh.V(i, 2) << "\n";

    // Write texture coordinates
    for (int i = 0; i < TC.rows(); ++i)
        out << "vt " << TC(i, 0) << " " << TC(i, 1) << "\n";

    // Write faces with texture indices
    for (int i = 0; i < mesh.F.rows(); ++i) {
        out << "f";
        for (int j = 0; j < 3; ++j) {
            int v_idx = mesh.F(i, j) + 1;    // OBJ uses 1-based indexing
            int vt_idx = FTC(i, j) + 1;
            out << " " << v_idx << "/" << vt_idx;
        }
        out << "\n";
    }

    out.close();
    return true;
}
}

#endif /* writeOBJ_h */
