// This file is part of libdirectional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_SINGULARITY_SPHERES_H
#define DIRECTIONAL_SINGULARITY_SPHERES_H

#include <cmath>
#include <igl/igl_inline.h>
#include <igl/barycenter.h>
#include <igl/per_face_normals.h>
#include <directional/representative_to_raw.h>
#include <directional/point_spheres.h>
#include <Eigen/Core>
#include <igl/avg_edge_length.h>


namespace directional
{
  // Returns a list of faces, vertices and color values that can be used to draw singularities for non-zero index values.
  // Input:
  //  V:              #V X 3 vertex coordinates.
  //  F:              #F X 3 mesh triangles.
  //  indices:        #V x 1 index (/N) per vertex (must be 0<index<N-1)
  //  positiveColors: N x 3 colos per positive index
  // negativeColors:  N x 3 colos per negative index
  // colorPerVertex in the output mesh
  // extendMesh   if to extend the singV,singT,singC, or to overwrite them
  // Output:
  //  singV:          The vertices of the singularity spheres.
  //  singF:          The faces of the singularity spheres.
  //  singC:         The colors of the singularity spheres.
  void IGL_INLINE singularity_spheres(const Eigen::MatrixXd& V,
                                      const Eigen::MatrixXi& F,
                                      const Eigen::VectorXi& singPositions,
                                      const Eigen::VectorXi& singIndices,
                                      const Eigen::MatrixXd& positiveColors,
                                      const Eigen::MatrixXd& negativeColors,
                                      const bool colorPerVertex,
                                      const bool extendMesh,
                                      Eigen::MatrixXd& singV,
                                      Eigen::MatrixXi& singF,
                                      Eigen::MatrixXd& singC)
  
  {

    Eigen::MatrixXd points(singPositions.size(), 3);
    Eigen::MatrixXd colors(singIndices.size(), 3);
    for (int i = 0; i < singIndices.rows(); i++)
    {
      points.row(i) = V.row(singPositions(i));
      if (singIndices(i) > 0)
        colors.row(i) = positiveColors.row((singIndices(i)-1 > positiveColors.rows()-1 ? positiveColors.rows()-1  : singIndices(i)-1) );
      else
        colors.row(i) = negativeColors.row((-singIndices(i)-1 > positiveColors.rows()-1 ? positiveColors.rows()-1  : -singIndices(i)-1));
      
    }
    double l = igl::avg_edge_length(V, F);
    directional::point_spheres(points, l/5.0, colors, 8, colorPerVertex, extendMesh, singV, singF, singC);
  
  }
  

  //version that provides all vertex indices instead of only singularities
  void IGL_INLINE singularity_spheres(const Eigen::MatrixXd& V,
                                      const Eigen::MatrixXi& F,
                                      const Eigen::VectorXi& fullIndices,
                                      const Eigen::MatrixXd& positiveColors,
                                      const Eigen::MatrixXd& negativeColors,
                                      const bool colorPerVertex,
                                      const bool extendMesh,
                                      Eigen::MatrixXd& singV,
                                      Eigen::MatrixXi& singF,
                                      Eigen::MatrixXd& singC)
  
  {
    
    std::vector<int> singPositionsList;
    std::vector<int> singIndicesList;
    for (int i=0;i<V.rows();i++)
      if (fullIndices(i)!=0){
        singPositionsList.push_back(i);
        singIndicesList.push_back(fullIndices(i));
      }
    
    Eigen::VectorXi singPositions(singPositionsList.size());
    Eigen::VectorXi singIndices(singIndicesList.size());
    for (int i=0;i<singPositionsList.size();i++){
      singPositions(i)=singPositionsList[i];
      singIndices(i)=singIndicesList[i];
    }
    
    singularity_spheres(V,singPositions,singIndices,positiveColors,negativeColors,colorPerVertex,extendMesh,singV, singF, singC);
  }
  
}

#endif
