// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2016 Francisca Gil Ureta <gilureta@cs.nyu.edu>, 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_STREAMLINES_H
#define DIRECTIONAL_STREAMLINES_H

#include <Eigen/Core>
#include <vector>
#include <igl/igl_inline.h>
#include <directional/CartesianField.h>
#include <directional/IntrinsicFaceTangentBundle.h>


namespace directional
{
#define SL_VERTEX 0
#define SL_EDGE 1
#define SL_FACE 2


  struct StreamlineData
  {
    directional::CartesianField field;
    directional::IntrinsicFaceTangentBundle stb;  //the subdivision tangent bundle which is used for the streamlining
    const TriMesh* slMesh;
    Eigen::MatrixXd slField;  //raw extrinsic field extracted by sampling
    //      (N degrees stacked horizontally for each triangle)
    //Eigen::MatrixXi match_ab;   //  #E by N matrix, describing for each edge the matching a->b, where a
    //      and b are the faces adjacent to the edge (i.e. vector #i of
    //      the vector set in a is matched to vector #mab[i] in b)
    // Eigen::MatrixXi match_ba;   //  #E by N matrix, describing the inverse relation to match_ab
    Eigen::VectorXi sampleFaces;    //all original faces
    Eigen::MatrixXd samplePoints;  //3d point that must lie on the respective faces
  };
  
  struct StreamlineState
  {

    //current element details (index, type, direction)
    Eigen::VectorXi currElements;
    Eigen::VectorXi currElementTypes;
    Eigen::VectorXi currDirectionIndex;  //TODO: what to do if this is not a face!
    Eigen::MatrixXd currStartPoints;
    Eigen::VectorXd currTimes;
    Eigen::VectorXi currSegmentIndex;  //index within the traced segments

    //next element details
    Eigen::VectorXi nextElements;
    Eigen::VectorXi nextElementTypes;
    Eigen::MatrixXd nextStartPoints;
    Eigen::VectorXd nextTimes;
    Eigen::VectorXi nextDirectionIndex;
    Eigen::Matrix<bool,Eigen::Dynamic, 1> segmentAlive;

    //current time (live segments are such that beginTimes <= currTime < endTime
    double currTime;

    std::vector<Eigen::RowVector3d> segStart, segEnd, segNormal;            //traced segments features
    std::vector<int> segOrigFace, segOrigVector;                   //original vectors and faces
    std::vector<double> segTimeSignatures;  //the time of the beginning of the segment

  };
  
  
  // Given a mesh and a field the function computes the /data/ necessary for tracing the field'
  // streamlines, and creates the initial /state/ for the tracing.
  // Input:
  //   field            Cartesian field to be traced.
  //   seedLocations    indices into F of the seeds for streaming. Can be Eigen::VectorXi() for automatic generation.
  //   ringDistance     Samples are generated automatically in case seedLocations.size()=0 according to exclusion of two seeds < ringDistance faces apart.
  // Output:
  //   data          struct containing topology information of the mesh and field
  //   state         struct containing the state of the tracing
  IGL_INLINE void streamlines_init(const directional::CartesianField& field,
                                   const Eigen::VectorXi& seedLocations,
                                   const double distRatio,
                                   StreamlineData &data,
                                   StreamlineState &state);


  // The function computes the next state for each point in the sample
  //   data          struct containing topology information
  //   state         struct containing the state of the tracing
  IGL_INLINE void streamlines_next(const StreamlineData & data,
                                   StreamlineState & state,
                                   const double dTime);
}

#include "streamlines.cpp"

#endif
