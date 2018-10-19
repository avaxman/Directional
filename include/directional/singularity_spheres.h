// This file is part of libdirectional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_SINGULARITY_SPHERES_H
#define DIRECTIONAL_SINGULARITY_SPHERES_H

#include <cmath>
#include <Eigen/Core>
#include <igl/igl_inline.h>
#include <igl/barycenter.h>
#include <igl/per_face_normals.h>
#include <igl/avg_edge_length.h>
#include <igl/jet.h>
#include <directional/representative_to_raw.h>
#include <directional/point_spheres.h>

namespace directional
{
  
  //returns the default libdirectional colors
  //supports up to N=6
  Eigen::MatrixXd IGL_INLINE defaultSingularityColors(int N){
    assert ((N>=1) && (N<=6));
    Eigen::MatrixXd fullColors;
    Eigen::VectorXd NList(2*N);
    for (int i=0;i<N;i++){
      NList(i)=-N+i;
      NList(N+i)=i+1;
    }
    igl::jet(-NList,true,fullColors);
    return fullColors;
  }
  
  // Returns a list of faces, vertices and color values that can be used to draw singularities for non-zero index values.
  // Input:
  //  V:              #V X 3 vertex coordinates.
  //  F:              #F X 3 mesh triangles.
  //  indices:        #V x 1 index (/N) per vertex (must be 0<index<N-1)
  //  singularityColors: 2*N x 3 colos per positive index in order [-N,..-1, 1, N]
  // Output:
  //  singV:          The vertices of the singularity spheres.
  //  singF:          The faces of the singularity spheres.
  //  singC:         The colors of the singularity spheres.
  void IGL_INLINE singularity_spheres(const Eigen::MatrixXd& V,
                                      const Eigen::MatrixXi& F,
                                      const Eigen::VectorXi& singPositions,
                                      const Eigen::VectorXi& singIndices,
                                      const Eigen::MatrixXd& singularityColors,
                                      Eigen::MatrixXd& singV,
                                      Eigen::MatrixXi& singF,
                                      Eigen::MatrixXd& singC)
  
  {

    Eigen::MatrixXd points(singPositions.size(), 3);
    Eigen::MatrixXd colors(singIndices.size(), 3);
    Eigen::MatrixXd positiveColors=singularityColors.block(singularityColors.rows()/2,0,singularityColors.rows()/2,3);
    Eigen::MatrixXd negativeColors=singularityColors.block(0,0,singularityColors.rows()/2,3);
    for (int i = 0; i < singIndices.rows(); i++)
    {
      points.row(i) = V.row(singPositions(i));
      if (singIndices(i) > 0)
        colors.row(i) = positiveColors.row((singIndices(i)-1 > positiveColors.rows()-1 ? positiveColors.rows()-1  : singIndices(i)-1) );
      else if (singIndices(i)<0)
        colors.row(i) = negativeColors.row((negativeColors.rows()+singIndices(i) > 0 ? negativeColors.rows()+singIndices(i) : 0));
      else
        colors.row(i).setZero(); //this shouldn't have been input
      
    }
    double l = igl::avg_edge_length(V, F);
    directional::point_spheres(points, l/5.0, colors, 8, singV, singF, singC);
  
  }
  

  //version that provides all vertex indices instead of only singularities
  void IGL_INLINE singularity_spheres(const Eigen::MatrixXd& V,
                                      const Eigen::MatrixXi& F,
                                      const Eigen::VectorXi& fullIndices,
                                      const Eigen::MatrixXd& singularityColors,
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
    
    singularity_spheres(V,F, singPositions,singIndices,singularityColors,singV, singF, singC);
  }
  
}

#endif
