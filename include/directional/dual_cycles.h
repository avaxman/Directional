// Copyright (C) 2016 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DUAL_CYCLES_H
#define DUAL_CYCLES_H
#include <Eigen/Core>
#include <igl/boundary_loop.h>
#include <igl/is_border_vertex.h>
#include <igl/colon.h>
#include <igl/setdiff.h>
#include <igl/slice.h>
#include <igl/unique.h>
#include <igl/edge_topology.h>
#include <vector>
#include <unordered_map>
#include "tree.h"


namespace directional
{
	// create a matrix the encodes the sums over the basis dual cycles in the mesh
	// the basis cycle matrix first contains #V cycles for each vertex (including 
	// empty cycles for border vertices), than #borders border cycles and finally 
	// 2*genus generator cycles around all handles.
	//input:
	//  F: #F by 3 triangles.
	//  EV: #E by 2 matrix of edges (vertex indices)
	//  EF: #E by 2 matrix of oriented adjacent faces.
	//output:
	//  basisCycleMat: #C by #E basis cycles (summing over edges)
	IGL_INLINE void dual_cycles(const Eigen::MatrixXi& F,
		const Eigen::MatrixXi& EV,
		const Eigen::MatrixXi& EF,
		Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycleMat
	)
	{
		using namespace Eigen;
		std::vector<bool> isBorder;

		int numV = F.maxCoeff() + 1;
		int eulerCharacteristic = numV - EV.rows() + F.rows();
		int g = (2 - eulerCharacteristic);
		std::cout << "g:" << g << std::endl;

		std::vector<Triplet<double> > basisCycleTriplets(EV.rows() * 2);

		if (g == 0) {  //no homology cycles
			basisCycleMat.resize(numV, EV.rows());

			//contractible (1-ring) cycles:
			for (int i = 0; i < EV.rows(); i++) {
				basisCycleTriplets.push_back(Triplet<double>(EV(i, 0), i, -1.0));
				basisCycleTriplets.push_back(Triplet<double>(EV(i, 1), i, 1.0));
			}
			basisCycleMat.setFromTriplets(basisCycleTriplets.begin(), basisCycleTriplets.end());
			basisCycleMat.prune(0, 0);
			return;
		}

		std::vector<std::vector<int>> boundaryLoops;

		igl::boundary_loop(F, boundaryLoops);

		if (boundaryLoops.size() == 0)
		{
			isBorder.resize(numV, false);
			//contractible (1-ring) cycles + boundary loops:
			for (int i = 0; i < EV.rows(); i++) {
				basisCycleTriplets.push_back(Triplet<double>(EV(i, 0), i, -1.0));
				basisCycleTriplets.push_back(Triplet<double>(EV(i, 1), i, 1.0));
			}
		}
		else
		{
			isBorder = igl::is_border_vertex(MatrixXi(numV, 0), F);
			//contractible (1-ring) cycles/boundary cycle:
			for (int i = 0; i < EV.rows(); i++) {
				if (EF(i, 0) == -1 || EF(i, 1) == -1)
					continue;
				std::unordered_map<int, std::vector<Triplet<double>>> border;

				if (isBorder[EV(i, 0)]) {
					if (!border.count(EV(i, 0)))
						border[EV(i, 0)] = std::vector<Triplet<double>>();
					border[EV(i, 0)].push_back(Triplet<double>(numV, i, -1.0));
				}
				else
					basisCycleTriplets.push_back(Triplet<double>(EV(i, 0), i, -1.0));

				if (isBorder[EV(i, 1)])
					border[EV(i, 1)].push_back(Triplet<double>(numV, i, 1.0));
				else
					basisCycleTriplets.push_back(Triplet<double>(EV(i, 1), i, 1.0));

				for (size_t i = 0; i < boundaryLoops.size(); ++i)
					for (size_t j = 0; j < boundaryLoops[i].size(); ++j)
					{
						for (size_t k = 0; k < border[boundaryLoops[i][j]].size(); k++)
							basisCycleTriplets.push_back(Triplet<double>(numV + i, border[boundaryLoops[i][j]][k].col(), border[boundaryLoops[i][j]][k].value()));
					}
			}
		}

		if (g == boundaryLoops.size())
		{
			basisCycleMat.resize(numV + boundaryLoops.size(), EV.rows());
			basisCycleMat.setFromTriplets(basisCycleTriplets.begin(), basisCycleTriplets.end());
			basisCycleMat.prune(0, 0);
			return;
		}

		MatrixXi reducedEV(EV);
		for (int i = 1; i < reducedEV.rows(); i++)
			if(isBorder[reducedEV(i,0)] || isBorder[reducedEV(i, 1)])
				reducedEV(i,0) = -1;


		VectorXi primalTreeEdges, primalTreeFathers;
		VectorXi dualTreeEdges, dualTreeFathers;
		tree(reducedEV, primalTreeEdges, primalTreeFathers);
		//creating a set of dual edges that do not cross edges in the primal tree
		VectorXi fullIndices = VectorXi::LinSpaced(EV.rows(), 0, EV.rows() - 1);
		VectorXi reducedEFIndices, inFullIndices;
		MatrixXi reducedEF;
		igl::setdiff(fullIndices, primalTreeEdges, reducedEFIndices, inFullIndices);
		VectorXi Two = VectorXi::LinSpaced(2, 0, 1);

		igl::slice(EF, reducedEFIndices, Two, reducedEF);
		tree(reducedEF, dualTreeEdges, dualTreeFathers);
		//converting dualTreeEdges from reducedEF to EF
		for (int i = 0; i < dualTreeEdges.size(); i++)
			dualTreeEdges(i) = inFullIndices(dualTreeEdges(i));


		for (int i = 0; i < dualTreeFathers.size(); i++)
			if (dualTreeFathers(i) != -1 && dualTreeFathers(i) != -2)
				dualTreeFathers(i) = inFullIndices(dualTreeFathers(i));

		//building tree co-tree based homological cycles
		//finding dual edge which are not in the tree, and following their faces to the end
		VectorXi isinTree = VectorXi::Zero(EF.rows());
		for (int i = 0; i < dualTreeEdges.size(); i++) {
			isinTree(dualTreeEdges(i)) = 1;
		}
		for (int i = 0; i < primalTreeEdges.size(); i++) {
			isinTree(primalTreeEdges(i)) = 1;
		}

		//std::cout<<"#free edges: "<<EF.rows()-isinTree.sum()<<std::endl;

		int numCycle = boundaryLoops.size();
		for (int i = 0; i < isinTree.size(); i++) {
			if (isinTree(i))
				continue;


			//std::cout<<"New Cycle"<<std::endl;
			//otherwise, follow both end faces to the root and this is the dual cycle
			if (EF(i, 0) == -1 || EF(i, 1) == -1)
				continue;
			basisCycleTriplets.push_back(Triplet<double>(numCycle + numV, i, 1.0));
			Vector2i currLeaves; currLeaves << EF(i, 0), EF(i, 1);
			VectorXi visitedOnce = VectorXi::Zero(EF.rows());  //used to remove the tail from the LCA to the root
			std::vector<Triplet<double> > candidateTriplets;
			for (int i = 0; i < 2; i++) { //on leaves
				int currTreeEdge = -1;  //indexing within dualTreeEdges
				int currFace = currLeaves(i);
				currTreeEdge = dualTreeFathers(currFace);
				if (currTreeEdge == -2)
				{
					numCycle--;
					break;
				}

				while (currTreeEdge != -1) {
					//determining orientation of current edge vs. face
					double sign = ((EF(currTreeEdge, 0) == currFace) != (i == 0) ? 1.0 : -1.0);
					visitedOnce(currTreeEdge) = 1 - visitedOnce(currTreeEdge);
					candidateTriplets.push_back(Triplet<double>(numCycle + numV, currTreeEdge, sign));
					currFace = (EF(currTreeEdge, 0) == currFace ? EF(currTreeEdge, 1) : EF(currTreeEdge, 0));
					currTreeEdge = dualTreeFathers(currFace);
				};
			}
			numCycle++;

			//only putting in dual edges that are below the LCA
			for (size_t i = 0; i < candidateTriplets.size(); i++)
				if (visitedOnce(candidateTriplets[i].col()))
					basisCycleTriplets.push_back(candidateTriplets[i]);
		}

		basisCycleMat.resize(numV + numCycle, EV.rows());
		basisCycleMat.setFromTriplets(basisCycleTriplets.begin(), basisCycleTriplets.end());
		basisCycleMat.prune(0, 0);
	}


	// create a matrix the encodes the sums over the basis dual cycles in the mesh
	// the basis cycle matrix first contains #V cycles for each vertex (including 
	// empty cycles for border vertices), than #borders border cycles and finally 
	// 2*genus generator cycles around all handles.
	//input:
	//  F: #F by 3 triangles.
	//output:
	//  basisCycleMat: #C by #E basis cycles (summing over edges)
	IGL_INLINE void dual_cycles(const Eigen::MatrixXi& F,
		Eigen::SparseMatrix<double, Eigen::RowMajor>& basisCycleMat
	)
	{
		Eigen::MatrixXi EV, x, EF;
		igl::edge_topology(Eigen::MatrixXi(F.maxCoeff(), 0), F, EV, x, EF);
		directional::dual_cycles(F, EV, EF, basisCycleMat);
	}
}
#endif
