// Copyright (C) 2017 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef GET_INDICES_H
#define GET_INDICES_H
#include <igl/igl_inline.h>
#include <igl/edge_topology.h>
#include <igl/parallel_transport_angles.h>
#include <igl/per_face_normals.h>
#include <igl/parallel_transport_angles.h>
#include <directional/cycle_curvature.h>

#include <Eigen/Core>
#include <vector>
#include <cmath>


namespace directional
{
	// Computes cycle-based indices from edge-based efforts.
    // Note: input is effort (sum of rotation angles), and not individual rotation angles
	// Input:
	//   basisCycleMat:     #basisCycles x #E the oriented basis cycles around which the indices are measured
	//   effort:            #basisCyles x 1 the effort required to match between two faces, equal to N*adjustment angles for N-RoSy fields.
	//   cycleCurvature:    #basisCycles x 1 the natural curvature (angle defect) for each basis cycle. e.g., from directional::cycle_curvature. Note: this is not in general just the cycle holonomy.
	//   N:                 The degree of the field
	// Output:
	//   singularities:     #basisCycles x 1 the index around each cycle times N (to get an integer).
	IGL_INLINE void get_indices(const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycleMat,
                                const Eigen::VectorXd& effort,
                                const Eigen::VectorXd& cycleCurvature,
                                const int N,
                                Eigen::VectorXi& indices)
	{

		Eigen::VectorXd dIndices = ((basisCycleMat * effort + N*cycleCurvature).array() / (2.0*igl::PI));  //this should already be an integer up to numerical precision
		indices.conservativeResize(dIndices.size());
		for (int i = 0; i < dIndices.size(); i++)
			indices(i) = (int)round(dIndices(i));
	}

	// Version without any precomputed curvature
	IGL_INLINE void get_indices(const Eigen::MatrixXd& V,
                                const Eigen::MatrixXi& F,
                                const Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycleMat,
                                const Eigen::VectorXd& effort,
                                int N,
                                Eigen::VectorXi& indices)
	{
		Eigen::VectorXd cycleCurvature;
		cycle_curvature(V, F, basisCycleMat, cycleCurvature);
		directional::get_indices(basisCycleMat, effort, cycleCurvature, N, indices);
	}

}




#endif


