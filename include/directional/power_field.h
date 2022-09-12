// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_POWER_FIELD_H
#define DIRECTIONAL_POWER_FIELD_H

#include <iostream>
#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <directional/polyvector_field.h>
#include <directional/TangentBundle.h>
#include <directional/CartesianField.h>


namespace directional
{
  // Computes a power field on the entire mesh from given values at the prescribed indices.
  // If no constraints are given the lowest-eigenvalue (of smoothness energy) field will be returned.
  // Inputs:
  //  mesh: a triangle mesh object.
  //  constFaces: the faces on which the polyvector is prescribed. If a face is repeated and the alignment is hard then all but the first vector in the face will be ignored.
  //  constVectors: #F by 3 in representative form of the N-RoSy's on the tangent spaces.
  //  alignWeights: #constFaces x 1 soft weights for alignment (negative values = fixed faces).
  //  N: The degree of the field.
  // Outputs:
  //  powerField: a cartesian power field object.
  IGL_INLINE void power_field(const TangentBundle& tb,
                              const Eigen::VectorXi& constSpaces,
                              const Eigen::MatrixXd& constVectors,
                              const Eigen::VectorXd& alignWeights,
                              const int N,
                              directional::CartesianField& field)
  {
    //TODO: have the field be initialized here
    polyvector_field(tb,constSpaces,constVectors,1.0, -1.0, alignWeights, N,field);
    field.fieldType = POWER_FIELD;
    //getting rid of the redundant zeros, in case they were allocated.
    field.intField.conservativeResize(field.intField.rows(),2);
    field.extField.conservativeResize(field.extField.rows(),3);
    //powerField=-pvField.col(0);  //powerfield is represented positively
  }
}


#endif
