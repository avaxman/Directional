#ifndef DIRECTIONAL_FACE_SUBDIVISION_H
#define DIRECTIONAL_FACE_SUBDIVISION_H
#include <vector>
#include <Eigen/Eigen>
#include "build_subdivision_operators.h"
#include "loop_coefficients.h"
namespace directional{
	namespace subdivision {
	template<typename CoefficientProvider>
		void Sf_triplet_provider(
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
			maxRow = 4 * F0.rows();
			const bool isBoundary = edges.size() & 1;

			// Number of faces
			const int fCount = (edges.size() - isBoundary) / 2;
			// Valence
			const int valence = (edges.size() + isBoundary) / 2; //One extra boundary vert.

			auto getFace = [&EF0, &edges, &edgeOrients](int faceIndex)
			{
				// Edge index associated with face
				const int i = 2 * faceIndex + 1;
				return EF0(edges[i], edgeOrients[i]);
			};
			// Corner number opposite spoke that is used to represent face.
			auto getCorner = [&EI0, &edges, &edgeOrients](int faceIndex)
			{
				// Edge index associated with face
				const int i = 2 * faceIndex + 1;
				return EI0(edges[i], edgeOrients[i]);
			};

			//Face finder helper. Allows circular index lookup.
			std::vector<int> inds;
			std::vector<double> coeffs;

			for (int fI = 0; fI < fCount; fI++)
			{
				inds.clear();
				coeffs.clear();
				const int f = getFace(fI);
				const int c = getCorner(fI);
				// Odd stencil
				const int target = 4 * f + c;
				coeffProvider(isBoundary, false, valence, fI, inds, coeffs);
				for (int i = 0; i < inds.size(); i++)
				{
					output.emplace_back(target, getFace(inds[i]), coeffs[i]);
				}
				maxRow = std::max(maxRow, target);

				// Avoid duplicate application of even stencil
				if (c != 0) continue;

				// Apply even stencil
				// Count number of adjacent faces
				int neighborCount = 0;
				std::vector<int> neighs;
				neighs.push_back(f);
				for(int i = 0; i < 3; i++)
				{
					// Incident face
					const int twin = EF0(SFE0(f, i), 1 - SFE0(f, 3 + i));
					// Is the twin face present, then increment neighbour count
					if (twin != -1) {
						neighborCount++;
						neighs.push_back(twin);
					}
				}
				inds.clear();
				coeffs.clear();
				// Find the coefficients and apply
				coeffProvider(neighborCount != 3, true, neighborCount, 0, inds, coeffs); // TODO check args are correct here
				const int evenTarget = 4 * f + 3;
				maxRow = std::max(maxRow, evenTarget);
				for (int i = 0; i < inds.size(); i++)
				{
					output.emplace_back(evenTarget, neighs[inds[i]], coeffs[i]);
				}
			}
		}
    }
}
#endif