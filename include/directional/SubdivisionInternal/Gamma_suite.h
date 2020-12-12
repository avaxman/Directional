// This file is part of Directional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef DIRECTIONAL_GAMMA_SUITE_H
#define DIRECTIONAL_GAMMA_SUITE_H

#include <queue>
#include <vector>
#include <cmath>
#include <Eigen/Core>
#include <igl/igl_inline.h>
#include <igl/local_basis.h>
#include <igl/edge_topology.h>
#include <directional/FEM_masses.h>
#include <igl/per_face_normals.h>
#include <directional/FEM_suite.h>
#include <igl/cat.h>


namespace directional
{
	/**
	 * Gamma3 based curl.
	 */
    inline void Gamma3_Curl(
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& EI,
		int faceCount,
		Eigen::SparseMatrix<double>& C)
	{
		using namespace Eigen;
		std::vector<Triplet<double>> trips;
		trips.reserve(faceCount * 3);
		for(int e = 0; e < EF.rows(); e++)
		{
			// No curl on boundary
			if (EF(e, 0) == -1 || EF(e, 1) == -1) continue;

			// Right gamma minus left gamma, relative to global edge orientation
			trips.emplace_back(e, 3 * EF(e, 1) + EI(e, 1), 1.0);
			trips.emplace_back(e, 3 * EF(e, 0) + EI(e, 0), -1.0);
		}
		C = SparseMatrix<double>(EF.rows(), 3 * faceCount);
		C.setFromTriplets(trips.begin(), trips.end());
	}

	/**
	 * Averaging integrated edge field to faces by simply summing contributions
	 */
	inline void Edge_To_Face_Average(const Eigen::MatrixXi& sFE, int edgeCount, Eigen::SparseMatrix<double>& A_EToF)
	{
		A_EToF = Eigen::SparseMatrix<double>(sFE.rows(), edgeCount);
		std::vector<Eigen::Triplet<double>> trips;
		trips.reserve(3 * sFE.rows());
		for(int f = 0; f < sFE.rows(); f++)
		{
			trips.emplace_back(f, sFE(f, 0), 1);
			trips.emplace_back(f, sFE(f, 1), 1);
			trips.emplace_back(f, sFE(f, 2), 1);
		}
		A_EToF.setFromTriplets(trips.begin(), trips.end());
	}

	/**
	 * Projects a PCVF into the gamma3 space by taking the innerproduct of the globally directed edges per face with the
	 * vector field value per face. Projections are stored at the index of the corner opposite the edge within the face.
	 * Input:
	 * V - |V| x 3 matrix of vertex locations
	 * F - |F| x 3 matrix of vertex indices per face
	 * EV - |E| x 2 matrix of edge to vertex mapping. Indices in the first column refer to start vertices, column 1 end vertices.
	 * sFE - |F| x 6 matrix of face to edge mappings, including orientation. For face f, the edge opposite the ith corner is stored
	 *			at sFE(f,i). Its orientation with respect to the face normal is located at sFE(f,i+3), where 0 denotes CCW and 1 CW.
	 * EF - |E| x 2 matrix of edge to face mapping. Elements in column 0 refere to face to the left of the edge (in global orientation),
	 *			for column 1 the faces to the right. For missing faces due to boundary, a -1 is used.
	 * Output:
	 * P - 3|F| x 3|F| sparse matrix mapping PCVFs to gamma3 elements. PCVF is assumed to be in [x1;y1;z1;x2;y2;z2;...etc] order in Matlab notation (column vec of xyz values)
	 *					Note: if the PCVF is truly tangent to the local face, then the sum of the gamma values per face is zero.
	 *
	 */
	IGL_INLINE void Gamma3_projector(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EF,
		Eigen::SparseMatrix<double>& P)
	{
		// Projects FEM vector field into Gamma3 space.
		P = Eigen::SparseMatrix<double>(3 * F.rows(), 3 * F.rows());

		std::vector<Eigen::Triplet<double>> trips;
		trips.reserve(9 * F.rows());
		for (int f = 0; f < F.rows(); f++) {
			// Iterate over corners
			for (int c = 0; c < 3; c++) {
				//Get edge opposite corner in global orientation
				const int eI = sFE(f, c);
				const int eS = sFE(f, c + 3);
				Eigen::RowVector3d edge = V.row(EV(eI, 1 - eS)) - V.row(EV(eI, eS));
				for (int i = 0; i < 3; i++) {
					trips.emplace_back(3 * f + c, 3 * f + i, edge[i]);
				}
			}
		}
		P.setFromTriplets(trips.begin(), trips.end());
	}

	/**
	 * Projects a PCVF into the gamma2 space by taking the innerproduct of the globally directed edges per face with the
	 * vector field value per face, omitting the third edge per face. Projections are stored at the index of the corner opposite
	 * the edge within the face.
	 * Input:
	 * V - |V| x 3 matrix of vertex locations
	 * F - |F| x 3 matrix of vertex indices per face
	 * EV - |E| x 2 matrix of edge to vertex mapping. Indices in the first column refer to start vertices, column 1 end vertices.
	 * sFE - |F| x 6 matrix of face to edge mappings, including orientation. For face f, the edge opposite the ith corner is stored
	 *			at sFE(f,i). Its orientation with respect to the face normal is located at sFE(f,i+3), where 0 denotes CCW and 1 CW.
	 * EF - |E| x 2 matrix of edge to face mapping. Elements in column 0 refere to face to the left of the edge (in global orientation),
	 *			for column 1 the faces to the right. For missing faces due to boundary, a -1 is used.
	 * Output:
	 * P - 2|F| x 3|F| sparse matrix mapping PCVFs to gamma3 elements. PCVF is assumed to be in [x1;y1;z1;x2;y2;z2;...etc] order in Matlab notation (column vec of xyz values)
	 *
	 */
	IGL_INLINE void Gamma2_projector(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EF,
		Eigen::SparseMatrix<double>& P)
	{
		// Projects FEM vector field into Gamma3 space.
		P = Eigen::SparseMatrix<double>(2 * F.rows(), 3 * F.rows());

		std::vector<Eigen::Triplet<double>> trips;
		trips.reserve(9 * F.rows());
		for (int f = 0; f < F.rows(); f++) {
			// Iterate over corners
			for (int c = 0; c < 2; c++) {
				//Get edge opposite corner in global orientation
				const int eI = sFE(f, c);
				const int eS = sFE(f, c + 3);
				Eigen::RowVector3d edge = V.row(EV(eI, 1)) - V.row(EV(eI, 0));
				for (int i = 0; i < 3; i++) {
					trips.emplace_back(2 * f + c, 3 * f + i, edge[i]);
				}
			}
		}
		P.setFromTriplets(trips.begin(), trips.end());
	}

	/**
	 * Reprojects a gamma2 field to a PCVF.
	 * Input:
	 * V - |V| x 3 matrix of vertex locations
	 * F - |F| x 3 matrix of vertex indices per face
	 * EV - |E| x 2 matrix of edge to vertex mapping. Indices in the first column refer to start vertices, column 1 end vertices.
	 * sFE - |F| x 6 matrix of face to edge mappings, including orientation. For face f, the edge opposite the ith corner is stored
	 *			at sFE(f,i). Its orientation with respect to the face normal is located at sFE(f,i+3), where 0 denotes CCW and 1 CW.
	 * EF - |E| x 2 matrix of edge to face mapping. Elements in column 0 refere to face to the left of the edge (in global orientation),
	 *			for column 1 the faces to the right. For missing faces due to boundary, a -1 is used.
	 * Output:
	 * P - 3|F| x 2|F| sparse matrix mapping PCVFs to gamma3 elements. PCVF is assumed to be in [x1;y1;z1;x2;y2;z2;...etc] order in Matlab notation (column vec of xyz values)
	 *
	 */
	IGL_INLINE void Gamma2_reprojector(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EF,
		Eigen::SparseMatrix<double>& RP)
	{
		using namespace Eigen;
		RP = Eigen::SparseMatrix<double>(3 * F.rows(), 2 * F.rows());
		std::vector<Triplet<double>> trips;
		trips.reserve(9 * F.rows());
		Eigen::VectorXd DA;
		Eigen::MatrixXd N;
		igl::doublearea(V, F, DA);
		igl::per_face_normals(V, F, N);
		for (int f = 0; f < F.rows(); f++)
		{
			RowVector3d localN = N.row(f);
			const int e0ind = sFE(f, 0), e1ind = sFE(f, 1);
            // Compute gamma signs
			const double s0 = 1.0 - 2.0 * sFE(f, 3), s1 = 1.0 - 2.0 * sFE(f, 4);
            // Get the edges
			RowVector3d e0 = V.row(EV(e0ind, 1)) - V.row(EV(e0ind, 0));
			RowVector3d e1 = V.row(EV(e1ind, 1)) - V.row(EV(e1ind, 0));
            // Compute basis elements
			RowVector3d vals0 = localN.cross(e0) / DA(f);
			RowVector3d vals1 = -localN.cross(e1) / DA(f);
			const double signProd = s0 * s1;
			for (int i = 0; i < 3; i++)
			{
				trips.emplace_back(3 * f + i, 2 * f, signProd * vals1(i));
				trips.emplace_back(3 * f + i, 2 * f + 1, signProd * vals0(i));
			}
		}
		RP.setFromTriplets(trips.begin(), trips.end());
	}

	IGL_INLINE void Gamma2_To_Gamma3(
		const Eigen::MatrixXi& sFE,
		Eigen::SparseMatrix<double>& G2_To_G3)
	{
		G2_To_G3 = Eigen::SparseMatrix<double>(3 * sFE.rows(), 2 * sFE.rows());
		using namespace Eigen;
		using namespace std;
		vector<Triplet<double>> trips;
		trips.reserve(4 * sFE.rows());
		for (int f = 0; f < sFE.rows(); f++)
		{
			RowVector3d signs = RowVector3d::Constant(1.0) - 2.0 * sFE.block(f, 3, 1, 3).cast<double>();
			trips.emplace_back(3 * f, 2 * f, 1.);
			trips.emplace_back(3 * f + 1, 2 * f + 1, 1.);
			// Reconstruct last halfedge form from zero sum constraint.
			trips.emplace_back(3 * f + 2, 2 * f, -signs(0)*signs(2));
			trips.emplace_back(3 * f + 2, 2 * f + 1, -signs(1)*signs(2));
		}
		G2_To_G3.setFromTriplets(trips.begin(), trips.end());
	}
	IGL_INLINE void Gamma3_To_Gamma2(
		const Eigen::MatrixXi& sFE,
		Eigen::SparseMatrix<double>& G3_To_G2)
	{
		G3_To_G2 = Eigen::SparseMatrix<double>(2 * sFE.rows(), 3 * sFE.rows());
		using namespace Eigen;
		using namespace std;
		vector<Triplet<double>> trips;
		trips.reserve(2 * sFE.rows());
		for (int f = 0; f < sFE.rows(); f++)
		{
			trips.emplace_back(2 * f, 3 * f, 1.);
			trips.emplace_back(2 * f + 1, 3 * f + 1, 1.);
			// Drop third halfedge form element per face
		}
		G3_To_G2.setFromTriplets(trips.begin(), trips.end());
	}

	/**
	 * Averages the gamma2 form to one form representation by taking the average of adjacent halfedge forms to the edge.
	 */
	IGL_INLINE void Gamma3_To_Oneform(
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& EI,
		int faceCount,
		Eigen::SparseMatrix<double>& G3_To_OneForm)
	{
		using namespace std;
		using namespace Eigen;
		G3_To_OneForm = Eigen::SparseMatrix<double>(EF.rows(), faceCount * 3);
		vector<Triplet<double>> trips;
		trips.reserve(3 * faceCount);
		for (int e = 0; e < EF.rows(); e++)
		{
			const int isBoundary = EF(e, 0) == -1 || EF(e, 1) == -1;
			if(EF(e,0) != -1)
			{
				trips.emplace_back(e, 3 * EF(e, 0) + EI(e, 0), isBoundary ? 1.0 : 0.5);
			}
			if(EF(e,1) != -1)
			{
				trips.emplace_back(e, 3 * EF(e, 1) + EI(e, 1), isBoundary ? 1.0 : 0.5);
			}
		}
		G3_To_OneForm.setFromTriplets(trips.begin(), trips.end());
	}
	/**
	 * Averages the gamma2 form to one form representation by taking the average of adjacent halfedge forms to the edge.
	 */
	IGL_INLINE void Gamma2_To_Oneform(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& EI,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EF,
		Eigen::SparseMatrix<double>& G2_To_OneForm)
	{
		using namespace std;
		using namespace Eigen;
		SparseMatrix<double> G3ToOneForm, G2ToG3;
		Gamma3_To_Oneform(EF, EI, sFE.rows(), G3ToOneForm);
		Gamma2_To_Gamma3(sFE, G2ToG3);
		G2_To_OneForm = G3ToOneForm * G2ToG3;

		//G2_To_OneForm = Eigen::SparseMatrix<double>(EV.rows(), sFE.rows() * 2);
		//vector<Triplet<double>> trips;
		//trips.reserve(2 * EV.rows());
		//for (int f = 0; f < sFE.rows(); f++)
		//{
		//	// Signs
		//	const double s0 = 1. - 2. * sFE(f, 3);
		//	const double s1 = 1. - 2. * sFE(f, 4);
		//	const double s2 = 1. - 2. * sFE(f, 5);
		//	trips.emplace_back(sFE(f, 0), 2 * f, 0.5);
		//	trips.emplace_back(sFE(f, 1), 2 * f + 1, 0.5);
		//	// "Reconstruct" third halfedge form and assign half to one form.
		//	trips.emplace_back(sFE(f, 2), 2 * f, -s0 * s2 * 0.5);
		//	trips.emplace_back(sFE(f, 2), 2 * f + 1, -s1 * s2 * 0.5);
		//}
		//G2_To_OneForm.setFromTriplets(trips.begin(), trips.end());
	}
	
	inline void Gamma3_To_Decomp(const Eigen::MatrixXi& EF, const Eigen::MatrixXi& EI, int faceCount, Eigen::SparseMatrix<double>& G3ToDC)
	{
		using SMat = Eigen::SparseMatrix<double>;
		SMat G3ToOneForm, C;
		Gamma3_To_Oneform(EF, EI , faceCount, G3ToOneForm);
		Gamma3_Curl(EF, EI, faceCount, C);
		assert(C.rows() == G3ToOneForm.rows());
		igl::cat(1, G3ToOneForm, Eigen::SparseMatrix<double>(0.5 * C), G3ToDC);
	}

    /**
     * \brief Projects the average oneform-half curl decomposition back to a gamma3 field.
     * \param EF The edge-to-face connectivity
     * \param EI The edge indexing in its adjacent faces
     * \param faceCount The number of faces
     * \param DCToG3 Output decomposition-to-gamma3 matrix operator
     */
    inline void Decomp_To_Gamma3(const Eigen::MatrixXi& EF, const Eigen::MatrixXi& EI, int faceCount, Eigen::SparseMatrix<double>& DCToG3)
	{
		using SMat = Eigen::SparseMatrix<double>;
		SMat OneFormToG3(3 * faceCount,EF.rows()), CToG3(3*faceCount, EF.rows());

		std::vector<Eigen::Triplet<double>> tripsCBack, tripsOneFormBack;
		tripsCBack.reserve(2 * EF.rows());
		tripsOneFormBack.reserve(2 * EF.rows());
		for(int e = 0; e < EF.rows(); e++)
		{
			const bool isBoundary = EF(e, 0) == -1 || EF(e, 1) == -1;
			if(isBoundary)
			{
				const int s = EF(e, 0) == -1 ? 1 : 0;
				tripsOneFormBack.emplace_back(3 * EF(e, s) + EI(e, s), e, 1.);
			}
			else
			{
				tripsOneFormBack.emplace_back(3 * EF(e, 0) + EI(e, 0), e, 1.);
				tripsOneFormBack.emplace_back(3 * EF(e, 1) + EI(e, 1), e, 1.);
				tripsCBack.emplace_back(3 * EF(e, 0) + EI(e, 0), e, -1.);
				tripsCBack.emplace_back(3 * EF(e, 1) + EI(e, 1), e, 1.);
			}
		}
		OneFormToG3.setFromTriplets(tripsOneFormBack.begin(), tripsOneFormBack.end());
		CToG3.setFromTriplets(tripsCBack.begin(), tripsCBack.end());
		igl::cat(2, OneFormToG3, CToG3, DCToG3);
	}


	// Creating non-conforming mid-edge mesh, where the faces are between the midedges of each original face. This is generally only for visualization
	  // Input:
	  //  VMesh:      #V x 3 conforming mesh vertices
	  //  FMesh:      #F x 3 conforming mesh faces, defines vertices per face, where F(f,i) is the vertex in face f at the ith corner.
	  //  EV:         #E x 2 edges to vertices indices, defines the global orientation of the edge e, which points from EV(e,0) to EV(e,1)
	  //  EF:         #E x 2 edges to faces indices, where EF(e,0) is the face left of the edge (in global orientation) and EF(e,1) the fact to the right.
	  //  sFE:         #F x 6 faces to edges indices, where edge sFE(f,0) is opposite corner 0 in face F (so has vertices F(f,1) and F(f,2)). sFE(f,3) tells whether this edge
	  //                is oriented CCW with respect to the face normal (sFE(f,3) == 0) or CW (sFE(f,3) == 1). The orientation is defined by the entry of the edge in EV (points from vertex in column 0 to vertex in column 1)
	  // Output:
	  //  Gv:    #3f x V Conforming gradient matrix, returning vector of xyzxyz per face gradient vectors
	  //  Ge:    #3f x V Non-conforming gradient of the same style, but for mid-edge functions
	  //  J:    #3F x 3F rotation operator [Nx] per face
	  //  C:  Curl operator which is basically (JGe)^T * Mchi
	  //  C:  Divergence operator which is basically Gv^T * Mchi

	IGL_INLINE void Gamma2_suite(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& EI,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EF,
		Eigen::SparseMatrix<double>& Gv,
		Eigen::SparseMatrix<double>& Ge,
		Eigen::SparseMatrix<double>& J,
		Eigen::SparseMatrix<double>& C,
		Eigen::SparseMatrix<double>& D,
		Eigen::SparseMatrix<double>& A_G2_To_E,
		Eigen::SparseMatrix<double>& G2ToDC,
		Eigen::SparseMatrix<double>& DCToG2)
	{
		using namespace std;
		using namespace Eigen;

		// Construct the projector
		SparseMatrix<double> P2;
		Gamma2_projector(V, F, EV, sFE, EF, P2);
		// Construct the reprojector
		SparseMatrix<double> RP2;
		Gamma2_reprojector(V, F, EV, sFE, EF, RP2);

		// Acquire the fem suite
		FEM_suite(V, F, EV, sFE, EF, Gv, Ge, J, C, D);

		// Build decomposition matrices
		SparseMatrix<double> G3To1Form, G3ToDC, DCToG3;
		directional::Gamma3_To_Decomp(EF, EI, F.rows(), G3ToDC);
		directional::Decomp_To_Gamma3(EF, EI, F.rows(), DCToG3);
		directional::Gamma3_To_Oneform(EF, EI, sFE.rows(), G3To1Form);

		// Augment the operators with the proper projection/reprojection
		Gv = P2 * Gv;
		Ge = P2 * Ge;
		J = P2 * J * RP2;
		D = D * RP2;

		//Averager to oneforms
		SparseMatrix<double> G2_To_G3, G3_To_G2;

		// Conversions between Gamma2/Gamma3 spaces.
		Gamma2_To_Gamma3(sFE, G2_To_G3);
		Gamma3_To_Gamma2(sFE, G3_To_G2);

		// Make operators work with Gamma2 space.
		C =  C *  G2_To_G3;
		A_G2_To_E = G3To1Form * G2_To_G3;
		G2ToDC = G3ToDC * G2_To_G3;
		DCToG2 = G3_To_G2 * DCToG3;
	}

	// Creating non-conforming mid-edge mesh, where the faces are between the midedges of each original face. This is generally only for visualization
	// Input:
	//  VMesh:      #V x 3 conforming mesh vertices
	//  FMesh:      #F x 3 conforming mesh faces, defines vertices per face, where F(f,i) is the vertex in face f at the ith corner.
	//  EV:         #E x 2 edges to vertices indices, defines the global orientation of the edge e, which points from EV(e,0) to EV(e,1)
	//  EF:         #E x 2 edges to faces indices, where EF(e,0) is the face left of the edge (in global orientation) and EF(e,1) the fact to the right.
	//  sFE:         #F x 6 faces to edges indices, where edge sFE(f,0) is opposite corner 0 in face F (so has vertices F(f,1) and F(f,2)). sFE(f,3) tells whether this edge
	//                is oriented CCW with respect to the face normal (sFE(f,3) == 0) or CW (sFE(f,3) == 1). The orientation is defined by the entry of the edge in EV (points from vertex in column 0 to vertex in column 1)
	// Output:
	//  Gv:    #3f x V Conforming gradient matrix, returning vector of xyzxyz per face gradient vectors
	//  Ge:    #3f x V Non-conforming gradient of the same style, but for mid-edge functions
	//  J:    #3F x 3F rotation operator [Nx] per face
	//  C:  Curl operator which is basically (JGe)^T * Mchi
	//  C:  Divergence operator which is basically Gv^T * Mchi

	IGL_INLINE void Gamma3_suite(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& EI,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EF,
		Eigen::SparseMatrix<double>& Gv,
		Eigen::SparseMatrix<double>& Ge,
		Eigen::SparseMatrix<double>& J,
		Eigen::SparseMatrix<double>& C,
		Eigen::SparseMatrix<double>& D,
		Eigen::SparseMatrix<double>& A_G3_To_E,
		Eigen::SparseMatrix<double>& G3ToDC,
		Eigen::SparseMatrix<double>& DCToG3)
	{
        
		using namespace Eigen;
		using namespace std;
		SparseMatrix<double> G2ToG3, G3ToG2;
		Gamma2_To_Gamma3(sFE, G2ToG3);
		Gamma3_To_Gamma2(sFE, G3ToG2);
		Gamma2_suite(V, F, EV, EI, sFE, EF, Gv, Ge, J, C, D, A_G3_To_E, G3ToDC, DCToG3);
		Gv = G2ToG3 * Gv;
		Ge = G2ToG3 * Ge;
		J = G2ToG3 * J * G3ToG2;
		C = C * G3ToG2;
		D = D * G3ToG2;
		A_G3_To_E = A_G3_To_E * G3ToG2;
		G3ToDC = G3ToDC * G3ToG2;
		DCToG3 = G2ToG3 * DCToG3;
	}

	inline void DEC_d0(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EF,
		Eigen::SparseMatrix<double>& D0)
	{
		D0 = Eigen::SparseMatrix<double>(EV.rows(), V.rows());
		std::vector<Eigen::Triplet<double>> trips;
		trips.reserve(2 * EV.rows());
		for (int e = 0; e < EV.rows(); e++) {
			trips.emplace_back(e, EV(e, 0), -1);
			trips.emplace_back(e, EV(e, 1), 1);
		}
		D0.setFromTriplets(trips.begin(), trips.end());
	}
	inline void DEC_d1(const Eigen::MatrixXd& V,
		const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EF,
		Eigen::SparseMatrix<double>& D1)
	{
		D1 = Eigen::SparseMatrix<double>(sFE.rows(), EV.rows());
		std::vector<Eigen::Triplet<double>> trips;
		trips.reserve(3 * sFE.rows());
		for (int f = 0; f < sFE.rows(); f++) {
			for (int j = 0; j < 3; j++)
			{
				trips.emplace_back(f, sFE(f, j), 1.0 - 2.0 * sFE(f, j + 3));
			}
		}
		D1.setFromTriplets(trips.begin(), trips.end());
	}
}

#endif
