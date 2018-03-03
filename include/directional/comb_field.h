// Copyright (C) 2017a Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef COMB_FIELD_H
#define COMB_FIELD_H
#include <igl/igl_inline.h>
#include <igl/gaussian_curvature.h>
#include <igl/local_basis.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/edge_topology.h>
#include <directional/tree.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>

#include <Eigen/Core>
#include <vector>
#include <cmath>


namespace directional
{
  // Takes a field in raw form and computes both the principal effort and the consequent principal matching on every edge.
  // Note: equals to rotation_angle*N in N-RoSy fields
  // Important: if the Raw field in not CCW ordered (e.., resulting from all functions X_t_raw() in libdirectional), the result is unpredictable.
  // Input:
  //  V:      #V x 3 vertex coordinates
  //  F:      #F x 3 face vertex indices
  //  EV:     #E x 2 edges to vertices indices
  //  EF:     #E x 2 edges to faces indices
  //  raw:    The directional field, assumed to be ordered CCW, and in xyzxyzxyz...xyz (3*N cols) form. The degree is inferred by the size.
  // Output:
  // matching: #E matching function, where vector k in EF(i,0) matches to vector (i+matching(i))%N in EF(i,1). In case of boundary, there is a -1.
  //  effort: #E principal matching efforts.
  //TODO: also return matching
  IGL_INLINE void comb_field(const Eigen::MatrixXd& V,
                             const Eigen::MatrixXi& F,
                             const Eigen::MatrixXi& EV,
                                     const Eigen::MatrixXi& EF,
                                     const Eigen::MatrixXi& FE,
                                     const Eigen::MatrixXd& rawField,
                             Eigen::MatrixXd combedField,
                             Eigen::VectorXi prinIndices,
                                     Eigen::VectorXi& combedMatching,
                                     Eigen::VectorXd& combedEffort)
  {
    using namespace Eigen;
    VectorXi matching, effort;
    principal_matching(V,F,EV, EF,FE,rawField,matching,effort);
    //flood-filling through the matching to comb field
    combedField.conservativeResize(rawField.size());
    
    //dual tree to find combing routes
    VectorXi tE, fathers;
    tree(EF,tE,fathers);
    

    
  }
  
  //representative version
  IGL_INLINE void comb_field(const Eigen::MatrixXd& V,
                                     const Eigen::MatrixXi& F,
                                     const Eigen::MatrixXi& EV,
                                     const Eigen::MatrixXi& EF,
                                     const Eigen::MatrixXi& FE,
                                     const Eigen::MatrixXd& representativeField,
                                     const int N,
                             Eigen::MatrixXd combedField,
                             Eigen::VectorXi prinIndices,
                             Eigen::VectorXi& combedMatching,
                             Eigen::VectorXd& combedEffort)
  {
    Eigen::MatrixXd raw;
    representative_to_raw(V, F, representativeField, N, raw);
    comb_field(V, F, EV, EF, FE, raw, combedField, prinIndices, matching, combedMatching, combedEffort);
  }
  
  
  
  // (V, F) only raw version
  IGL_INLINE void comb_field(const Eigen::MatrixXd& V,
                                     const Eigen::MatrixXi& F,
                                     const Eigen::MatrixXd& raw,
                                     Eigen::MatrixXd combedField,
                                     Eigen::VectorXi prinIndices,
                                     Eigen::VectorXi& combedMatching,
                                     Eigen::VectorXd& combedEffort)
  {
    Eigen::MatrixXi EV, FE, EF;
    igl::edge_topology(V, F, EV, FE, EF);
    comb_field(V, F, EV, EF, FE, raw, combedField, prinIndices, matching, combedMatching, combedEffort);
  }
  
  // (V,F) only representative version
  IGL_INLINE void comb_field(const Eigen::MatrixXd& V,
                                     const Eigen::MatrixXi& F,
                                     const Eigen::MatrixXd& representativeField,
                                     const int N,
                                     Eigen::MatrixXd combedField,
                                     Eigen::VectorXi prinIndices,
                                     Eigen::VectorXi& combedMatching,
                                     Eigen::VectorXd& combedEffort)
  {
    Eigen::MatrixXi EV, FE, EF;
    igl::edge_topology(V, F, EV, FE, EF);
    comb_field(V, F, EV, EF, FE, representativeField, N combedField, prinIndices, matching, combedMatching, combedEffort);
  }
}




#endif


