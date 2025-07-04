// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_READ_OFF_H
#define DIRECTIONAL_READ_OFF_H
#include <cmath>
#include <Eigen/Core>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>
#include <cassert>
#include <directional/TriMesh.h>

namespace directional
{
//Reading an OFF file into a TriMesh class
bool inline readOFF(const std::string off_file_name,
                    directional::TriMesh& mesh){
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    std::ifstream file(off_file_name);
    if (!file.is_open()) {
        std::cerr << "Failed to open OFF file" << std::endl;
        return false;
    }
    
    std::string line;
    int numVertices, numFaces;
    file >> line >> numVertices >> numFaces;
    std::getline(file, line); // Skip remaining characters on first line
    
    V.resize(numVertices, 3);
    for (int i = 0; i < numVertices; ++i)
        file >> V(i,0) >> V(i,1) >> V(i,2);
    
    F.resize(numFaces,3);
    for (int i = 0; i < numFaces; ++i) {
        int numVerticesPerFace;
        file >> numVerticesPerFace;
        assert(numVerticesPerFace==3 && "readOFF() is only intended for triangle meshes.");
        file >> F(i,0) >> F(i,1) >> F(i,2);
    }
    mesh.set_mesh(V,F);
    return true;
}
}

#endif
