// This file is part of Directional, a library for directional field processing.
//
// Copyright (C) 2016 Francisca Gil Ureta <gilureta@cs.nyu.edu>, 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#include <Eigen/Geometry>
#include <igl/edge_topology.h>
#include <igl/sort_vectors_ccw.h>
#include <igl/per_face_normals.h>
#include <igl/segment_segment_intersect.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/barycenter.h>
#include <igl/slice.h>
#include <directional/principal_matching.h>
#include <directional/streamlines.h>


IGL_INLINE void directional::streamlines_init(
                                      const Eigen::MatrixXd V,
                                      const Eigen::MatrixXi F,
                                      const Eigen::MatrixXd &temp_field,
                                      StreamlineData &data,
                                      StreamlineState &state,
                                      double percentage
                                      ){
  using namespace Eigen;
  using namespace std;
  
  igl::edge_topology(V, F, data.EV, data.FE, data.EF);
  igl::triangle_triangle_adjacency(F, data.TT);
  
  // prepare vector field
  // --------------------------
  int degree = temp_field.cols()/3;
  data.degree = degree;
  
  Eigen::MatrixXd FN;
  Eigen::VectorXi order;
  Eigen::RowVectorXd sorted;
  
  igl::per_face_normals(V, F, FN);
  data.field.setZero(F.rows(), degree * 3);
  for (unsigned i = 0; i < F.rows(); ++i){
    const Eigen::RowVectorXd &n = FN.row(i);
    Eigen::RowVectorXd temp(1, degree * 3);
    temp = temp_field.row(i);
    igl::sort_vectors_ccw(temp, n, order, sorted);
    
    // project vectors to tangent plane
    for (int j = 0; j < degree; ++j)
    {
      Eigen::RowVector3d pd = sorted.segment(j * 3, 3);
      pd = (pd - (n.dot(pd)) * n).normalized();
      data.field.block(i, j * 3, 1, 3) = pd;
    }
  }
  Eigen::VectorXd effort;
  directional::principal_matching(V, F, data.EV, data.EF, data.FE, data.field, data.matching, data.effort);
  
  // create seeds for tracing
  // --------------------------
  Eigen::VectorXi samples;
  int nsamples;
  
  nsamples = percentage * F.rows();
  Eigen::VectorXd r;
  r.setRandom(nsamples, 1);
  r = (1 + r.array()) / 2.;
  samples = (r.array() * F.rows()).cast<int>();
  data.nsample = nsamples;
  
  Eigen::MatrixXd BC, BC_sample;
  igl::barycenter(V, F, BC);
  igl::slice(BC, samples, 1, BC_sample);
  
  // initialize state for tracing vector field
  
  state.start_point = BC_sample.replicate(degree,1);
  state.end_point = state.start_point;
  
  state.current_face = samples.replicate(1, degree);
  
  state.current_direction.setZero(nsamples, degree);
  for (int i = 0; i < nsamples; ++i)
    for (int j = 0; j < degree; ++j)
      state.current_direction(i, j) = j;
  
}

IGL_INLINE void directional::streamlines_next(
                                      const Eigen::MatrixXd V,
                                      const Eigen::MatrixXi F,
                                      const StreamlineData & data,
                                      StreamlineState & state
                                      ){
  using namespace Eigen;
  using namespace std;
  
  int degree = data.degree;
  int nsample = data.nsample;
  
  state.start_point = state.end_point;
  
  for (int i = 0; i < degree; ++i)
  {
    for (int j = 0; j < nsample; ++j)
    {
      int f0 = state.current_face(j,i);
      if (f0 == -1) // reach boundary
        continue;
      int m0 = state.current_direction(j, i);
      
      // the starting point of the vector
      const Eigen::RowVector3d &p = state.start_point.row(j + nsample * i);
      // the direction where we are trying to go
      const Eigen::RowVector3d &r = data.field.block(f0, 3 * m0, 1, 3);
      
      
      // new state,
      int f1, m1;
      
      for (int k = 0; k < 3; ++k)
      {
        f1 = data.TT(f0, k);
        
        // edge vertices
        const Eigen::RowVector3d &q = V.row(F(f0, k));
        const Eigen::RowVector3d &qs = V.row(F(f0, (k + 1) % 3));
        // edge direction
        Eigen::RowVector3d s = qs - q;
        
        double u;
        double t;
        if (igl::segments_intersect(p, r, q, s, t, u))
        {
          // point on next face
          state.end_point.row(j + nsample * i) = p + t * r;
          state.current_face(j,i) = f1;
          
          // matching direction on next face
          int e1 = data.FE(f0, k);
          if (data.EF(e1, 0) == f0)
            m1 = (data.matching(e1)+m0)%data.degree; //m1 = data.match_ab(e1, m0);
          else
            m1 = (-data.matching(e1)+m0+data.degree)%data.degree;  //data.match_ba(e1, m0);
          
          state.current_direction(j, i) = m1;
          break;
        }
        
      }
    }
  }
}
