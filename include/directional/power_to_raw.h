// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POWER_TO_RAW_H
#define DIRECTIONAL_POWER_TO_RAW_H

#include <igl/local_basis.h>
#include <directional/rotation_to_representative.h>
#include <directional/representative_to_raw.h>
#include <directional/power_to_representative.h>
#include <directional/TriMesh.h>
#include <directional/FaceField.h>

namespace directional
{
  // Converts the power complex representation to raw representation.
  // Input:
  //  B1, B2, B3: bases for each face from igl::local_base().
  //  powerField: #F x 1 Representation of the field as complex numbers
  //  N: the degree of the field.
  // normalize: whether to produce a normalized result (length = 1)
  // Output:
  //  rawField: #F by 3*N matrix with all N explicit vectors of each directional in the order X,Y,Z,X,Y,Z, ...
  IGL_INLINE void power_to_raw(const directional::TriMesh& mesh,
                               const Eigen::MatrixXcd& powerField,
                               int N,
                               directional::FaceField& field,
                               bool normalize=false)
  {
    Eigen::MatrixXd representative,rawField;
    power_to_representative(mesh.Bx, mesh.By, powerField, N, representative);
    if (normalize)
      representative.rowwise().normalize();
    representative_to_raw(mesh.faceNormals, representative, N, rawField);
    field.set_field(rawField,mesh);
  }
  
}

#endif
