#ifndef DIRECTIONAL_VERTEX_SUBDIVISION_H
#define DIRECTIONAL_VERTEX_SUBDIVISION_H
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include "build_subdivision_operators.h"
#include "loop_coefficients.h"
namespace directional{
	namespace subdivision {
		template<typename CoefficientProvider>
		void Sv_triplet_provider(
			const int& vertexCount,
			const Eigen::MatrixXi& F0,
			const Eigen::MatrixXi& SFE0,
			const Eigen::MatrixXi& E0,
			const Eigen::MatrixXi& EI0,
			const Eigen::MatrixXi& EF0,
			const Eigen::MatrixXi& E0ToEk,
			const std::vector<int>& edges,
			const std::vector<int>& edgeOrients,
			const CoefficientProvider& coeffProvider,
			std::vector<Eigen::Triplet<double>>& output,
			int& maxRow
		)
		{
			maxRow = vertexCount + E0.rows();

			const bool isBoundary = edges.size() & 1;
			// Gets the ith vertex in the ring. 0 is the central vertex,
			// the others are numbered, staring at the boundary, moving in
			// CCW direction relative to the local normal.
			auto getVertex = [&edges, &edgeOrients, &E0](int vInd) {
				if (vInd == 0) return E0(edges[0], edgeOrients[0]);
				return E0(edges[2 * (vInd - 1)], 1 - edgeOrients[2 * (vInd - 1)]);
			};
			//Central vertex in loop
			const int vCentral = getVertex(0);

			const int valence = (edges.size() + isBoundary) / 2;
			const int vCount = vertexCount;

			// Retrieves the indices and coefficients required for even
			std::vector<int> inds; //Index zero is the central vertex of the ring
			std::vector<double> coeffs;
			coeffProvider(isBoundary, true, valence,0,  inds, coeffs);

			// Add the output
			for (int i = 0; i < inds.size(); i++) output.emplace_back(vCentral, getVertex(inds[i]), coeffs[i]);

			// Odd vertices 
			for (int v = 1; v <= valence; v++)
			{
				// Avoid duplicating stencils here
				if (edgeOrients[2 * (v - 1)] == 1) continue;

				// Clear output of previous coefficients
				inds.clear();
				coeffs.clear();

				// Get the odd stencil
				coeffProvider(isBoundary, false, valence, v, inds, coeffs);

				// The new row for the subdivided vertex is 
				const int target = vCount + edges[2 * (v - 1)];

				for (int i = 0; i < inds.size(); i++)
				{
					output.emplace_back(target, getVertex(inds[i]), coeffs[i]);
				}
			}
		}

	}
}
#endif