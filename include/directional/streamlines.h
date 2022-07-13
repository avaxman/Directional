// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2016 Francisca Gil Ureta <gilureta@cs.nyu.edu>, 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_STREAMLINES_H
#define DIRECTIONAL_STREAMLINES_H

#include <igl/igl_inline.h>
#include <directional/definitions.h>

#include <Eigen/Core>
#include <vector>

namespace directional
{
  struct StreamlineData
  {
    directional::FaceField field;
    //      (N degrees stacked horizontally for each triangle)
    //Eigen::MatrixXi match_ab;   //  #E by N matrix, describing for each edge the matching a->b, where a
    //      and b are the faces adjacent to the edge (i.e. vector #i of
    //      the vector set in a is matched to vector #mab[i] in b)
    // Eigen::MatrixXi match_ba;   //  #E by N matrix, describing the inverse relation to match_ab
    int nsample;                //  #S, number of sample points
    Eigen::VectorXi samples;    //all original faces
  };
  
  struct StreamlineState
  {
    Eigen::MatrixXd start_point;        //  #N*S by 3 starting points of segment (stacked vertically for each degree)
    Eigen::MatrixXd end_point;          //  #N*S by 3 endpoints points of segment (stacked vertically for each degree)
    Eigen::MatrixXi current_face;       //  #S by N face indices (stacked horizontally for each degree)
    Eigen::MatrixXi current_direction;  //  #S by N field direction indices (stacked horizontally for each degree)
    int numSteps;                       // number of steps taken so far
    
    Eigen::MatrixXd P1, P2;             //entire set of streamline segments
    Eigen::VectorXi origFace, origVector; //original vectors from faces
    Eigen::VectorXi timeSignature;        //time (in steps) of the current segment
    
  };
  
  
  // Given a mesh and a field the function computes the /data/ necessary for tracing the field'
  // streamlines, and creates the initial /state/ for the tracing.
  // Input:
  //   V                #V by 3 list of mesh vertex coordinates
  //   F                #F by 3 list of mesh faces
  //   rawField         #F by 3n list of the 3D coordinates of the per-face vectors (n-degrees stacked horizontally for each triangle)
  //   seedLocations    indices into F of the seeds for streaming. Can be Eigen::VectorXi() for automatic generation.
  //   ringDistance     Samples are generated automatically in case seedLocations.size()=0 according to exclusion of two seeds < ringDistance faces apart.
  // Output:
  //   data          struct containing topology information of the mesh and field
  //   state         struct containing the state of the tracing
  IGL_INLINE void streamlines_init(const directional::FaceField& field,
                                   const Eigen::VectorXi& seedLocations,
                                   const int ringDistance,
                                   StreamlineData &data,
                                   StreamlineState &state);


  
  // The function computes the next state for each point in the sample
  //   V             #V by 3 list of mesh vertex coordinates
  //   F             #F by 3 list of mesh faces
  //   data          struct containing topology information
  //   state         struct containing the state of the tracing
  IGL_INLINE void streamlines_next(const StreamlineData & data,
                                   StreamlineState & state
                                   );
}

#include "streamlines.cpp"

#endif
