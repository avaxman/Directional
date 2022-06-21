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
#include <igl/readOFF.h>
#include <directional/TriMesh.h>


namespace directional
{
  bool IGL_INLINE readOFF(const std::string obj_file_name, directional::TriMesh& mesh){
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOFF(obj_file_name,V,F);
    mesh.set_mesh(V,F);
  }
}

#endif
