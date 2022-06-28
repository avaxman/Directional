// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2022 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_WRITE_OBJ_H
#define DIRECTIONAL_WRITE_OBJ_H


#include <Eigen/Core>
#include <string>
#include <vector>
#include <igl/igl_inline.h>
#include <igl/writeOBJ.h>
#include <directional/TriMesh.h>

namespace directional
{

  bool writeOBJ(const std::string& fileName,
                const directional::TriMesh& mesh,
                const Eigen::MatrixXd& TC,
                const Eigen::MatrixXd& FTC)
  {
    Eigen::MatrixXd emptyMat;
    return igl::writeOBJ(fileName, mesh.V, mesh.F, emptyMat, emptyMat, TC, FTC);
  }

}

#endif /* writeOBJ_h */
