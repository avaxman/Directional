// Copyright (C) 2016 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef POLY_VECTOR_H
#define POLY_VECTOR_H

#include <Eigen/Geometry>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/Polynomials>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/local_basis.h>
#include <iostream>

namespace directional
{
	// Method to precalculate the solvers for a poly_vector field. Must be recalculated whenever 
	// soft_ids changes or the mesh changes.
	// Inputs:
	//  V: #V by 3 vertex coordinates.
	//  F: #F by 3 face vertex indices.
	//  TT: #F by 3 Triangle-triangle adjecencies.
	//  B1, B2: #F by 3 matrices representing the local base of each face.
	//  soft_id: The face ids of the soft constraints described in soft_value.
	//  N: The degree of the field.
	// Outputs:
	//  solvers: A vecetor of pre-computed solver pointers. Must be recomputed if soft_id changes.
	//  energy: The energy matrices for the given problem, must be recomputed along with the solvers.
	IGL_INLINE void poly_vector_prepare_solvers(
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& TT, 
		const Eigen::MatrixXd& B1,
		const Eigen::MatrixXd& B2,
		const Eigen::VectorXi& soft_id,
		const int N,
		std::vector<Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>*>& solvers,
		std::vector<Eigen::SparseMatrix<std::complex<double>>>& energy)
	{
		using namespace std;
		using namespace Eigen;

		for (int n = 0; n < N; n++)
		{
			if (solvers.size() <= n)
				solvers.push_back(new Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>());
			if (energy.size() <= n)
				energy.push_back(Eigen::SparseMatrix<std::complex<double>>());
			// Build the sparse matrix, with an energy term for each edge
			std::vector< Triplet<std::complex<double> > > t;

			int count = 0;
			for (int f = 0; f < F.rows(); ++f)
			{
				for (int ei = 0; ei < F.cols(); ++ei)
				{
					// Look up the opposite face
					int g = TT(f, ei);

					// If it is a boundary edge, it does not contribute to the energy
					if (g == -1) continue;

					// Avoid to count every edge twice
					if (f > g) continue;

					// Compute the complex representation of the common edge
					Vector3d e = (V.row(F(f, (ei + 1) % 3)) - V.row(F(f, ei)));
					Vector2d vef = Vector2d(e.dot(B1.row(f)), e.dot(B2.row(f))).normalized();
					std::complex<double> ef(vef(0), vef(1));
					Vector2d veg = Vector2d(e.dot(B1.row(g)), e.dot(B2.row(g))).normalized();
					std::complex<double> eg(veg(0), veg(1));

					// Add the term conj(f)^n*ui - conj(g)^n*uj to the energy matrix
					t.push_back(Triplet<std::complex<double> >(count, f, std::pow(std::conj(ef), N-n)));
					t.push_back(Triplet<std::complex<double> >(count, g, -1.*std::pow(std::conj(eg), N-n)));

					++count;
				}
			}
			double lambda = 1000000000000;
			for (int r = 0; r < soft_id.size(); ++r)
			{
				int f = soft_id(r);
				t.push_back(Triplet<std::complex<double> >(count, f, lambda));
				++count;
			}

			// Prepare the solver
			energy[n].resize(count, F.rows());
			energy[n].setFromTriplets(t.begin(), t.end());
			solvers[n]->compute(energy[n].adjoint()*energy[n]);
		}
	}


	// Version of the poly_vector method that allows reusing the solvers. Note that solvers need to
	// be recalculated whenever the soft_ids change.
	// Returns a polyvector field that attempts to make each vector as parallel as possible to the
	// example directionals given in soft_id and soft_values. Also known as "Globally Optimal" and 
	// "As Parallel As Possible". The field will be returned in the form of a complex polynomial 
	// and can be transformed into a raw vector field using the poly_to_raw function.
	// If no constraints are given the zero field will be returned.
	// Inputs:
	//  B1, B2: #F by 3 matrices representing the local base of each face.
	//  soft_id: The face ids of the soft constraints described in soft_value.
	//  soft_value: The directionals on the faces indicated by soft_id around which the field is
	//              generated. Should be in the form X1,Y1,Z1,X2,Y2,Z2,Xn,Yn,Zn.
	//  solvers: An aray of pre-computed solver pointers from poly_vector_prepare_solvers. 
	//           Must be recomputed if soft_id changes.
	//  energy: The energy matrices for the given problem, must be recomputed along with the solvers.
	//  N: The degree of the field.
	// Outputs:
	//  poly: Representation of the field as a complex polynomial
	IGL_INLINE void poly_vector(
		const Eigen::MatrixXd& B1,
		const Eigen::MatrixXd& B2,
		const Eigen::VectorXi& soft_id,
		const Eigen::MatrixXd& soft_value,
		const std::vector<Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>*>& solvers,
		const std::vector<Eigen::SparseMatrix<std::complex<double>>>& energy,
		const int N, 
		Eigen::MatrixXcd& poly)
	{
		using namespace std;
		using namespace Eigen;
		if (soft_id.size() == 0)
		{
			poly = MatrixXcd::Constant(B1.rows(), N, std::complex<double>());
			return;
		}

		// Build the sparse matrix, with an energy term for each edge
		std::vector< Triplet<std::complex<double> > > tb;

		int count = energy[0].rows() - soft_id.size();

		// Convert the constraints into the complex polynomial coefficients and add them as soft constraints
		double lambda = 1000000000000;
		for (int r = 0; r < soft_id.size(); ++r)
		{
			int f = soft_id(r);
			Eigen::RowVectorXcd roots(N);
			for (int i = 0; i < N; i++)
			{
				Vector3d v = soft_value.block<1, 3>(r, i * 3);
				std::complex<double> c(v.dot(B1.row(f)), v.dot(B2.row(f)));
				roots(i) = c;
			}
			Eigen::VectorXcd poly;
			roots_to_monicPolynomial(roots, poly);
			for (int n = 0; n < N; n++)
			{
				tb.push_back(Triplet<std::complex<double> >(count, n, poly(n) * std::complex<double>(lambda, 0)));
			}
			++count;
		}

		poly.resize(B1.rows(), N);
		typedef SparseMatrix<std::complex<double>> SparseMatrixXcd;
		SparseMatrixXcd b(count, N);
		b.setFromTriplets(tb.begin(), tb.end());

		for (int n = 0; n < N; n++)
		{
			// Solve the linear system
			assert(solvers[n]->info() == Success);
			poly.col(n) = solvers[n]->solve(energy[n].adjoint()*MatrixXcd(b).col(n));
			assert(solvers[n]->info() == Success);
		}
	}

	// Returns a polyvector field that attempts to make each vector as parallel as possible to the
	// example directionals given in soft_id and soft_values. Also known as "Globally Optimal" and 
	// "As Parallel As Possible". The field will be returned in the form of a complex polynomial 
	// and can be transformed into a raw vector field using the complex_to_polynomial function.
	// If no constraints are given the zero field will be returned.
	// Inputs:
	//  V: #V by 3 vertex coordinates.
	//  F: #F by 3 face vertex indices.
	//  TT: #F by 3 Triangle-triangle adjecencies.
	//  B1, B2: #F by 3 matrices representing the local base of each face.
	//  soft_id: The face ids of the soft constraints described in soft_value.
	//  soft_value: The directionals on the faces indicated by soft_id around which the field is
	//              generated. Should be in the form X,Y,Z.
	//  N: The degree of the field.
	// Outputs:
	//  poly: Representation of the field as a complex polynomial
	IGL_INLINE void poly_vector
	(
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& TT,
		const Eigen::MatrixXd& B1,
		const Eigen::MatrixXd& B2,
		const Eigen::VectorXi& soft_id,
		const Eigen::MatrixXd& soft_value,
		const int N,
		Eigen::MatrixXcd& poly
	)
	{
		std::vector<Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>*> solvers;
		std::vector<Eigen::SparseMatrix<std::complex<double>>> As;
		poly_vector_prepare_solvers(V, F, TT, B1, B2, soft_id, N, solvers, As);
		poly_vector(B1, B2, soft_id, soft_value, solvers, As, N, poly);
		for (std::vector< Eigen::SimplicialLDLT<Eigen::SparseMatrix<std::complex<double>>>* >::iterator it = solvers.begin(); it != solvers.end(); ++it)
		{
			delete (*it);
		}
	}


	// Returns a polyvector field that attempts to make each vector as parallel as possible to the
	// example directionals given in soft_id and soft_values. Also known as "Globally Optimal" and 
	// "As Parallel As Possible". The field will be returned in the form of a complex polynomial 
	// and can be transformed into a raw vector field using the complex_to_polynomial function.
	// If no constraints are given the zero field will be returned.
	// Inputs:
	//  V: #V X 3 vertex coordinates.
	//  F: #F by 3 face vertex indices.
	//  soft_id: The face ids of the soft constraints described in soft_value.
	//  soft_value: The directionals on the faces indicated by soft_id around which the field is
	//              generated. Should be in the form X,Y,Z.
	//  N: The degree of the field.
	// Outputs:
	//  poly: Representation of the field as a complex polynomial
	IGL_INLINE void poly_vector
	(
		const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::VectorXi& soft_id,
		const Eigen::MatrixXd& soft_value,
		const int N,
		Eigen::MatrixXcd& poly
	)
	{
		Eigen::MatrixXi TT;
		igl::triangle_triangle_adjacency(F, TT);
		Eigen::MatrixXd B1, B2, x;
		igl::local_basis(V, F, B1, B2, x);
		poly_vector(V, F, TT, B1, B2, soft_id, soft_value, N, poly);
	}
}
#endif