#ifndef DIRECTIONAL_ONEFORM_SUBDIVISION_H
#define DIRECTIONAL_ONEFORM_SUBDIVISION_H
#include <vector>
#include <Eigen/Eigen>
#include "build_subdivision_operators.h"
#include "loop_coefficients.h"
namespace directional{
	namespace subdivision {
        template<typename CoefficientProvider>
        void Se_triplet_provider(
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
            // Is the ring a boundary ring
            const bool isBoundaryRing = edges.size() & 1;

            std::vector<double> signs(edgeOrients.size(),1.0);
            for(int i = 0; i < signs.size(); i++) signs[i] -= 2.0 * edgeOrients[i];
            
            // Number of edges
            const int eCount = edges.size();
            // Vertex valence
            const int valence = (eCount + isBoundaryRing) / 2; //One extra boundary vert.

            //Local stencil description
            std::vector<int> inds;
            std::vector<double> coeffs;

            // All even elements
            bool isEven = true;
            for (int eI = 0; eI < eCount; eI += 2)
            {
                // Acquire stencil
                inds.clear();
                coeffs.clear();
                coeffProvider(isBoundaryRing, isEven, valence, eI, inds, coeffs);
                const int e = edges[eI];
                const double edgeSign = signs[eI];
                const int target = E0ToEk(e, edgeOrients[eI]);
                maxRow= std::max(maxRow,target+1);
                for (int i = 0; i < inds.size(); i++)
                {
                    output.emplace_back(target, edges[inds[i]], edgeSign * signs[inds[i]] * coeffs[i]);
                }
            }

            isEven = false;
            // All odd elements
            for (int eI = 1; eI < eCount; eI += 2)
            {
                // Avoid duplicates (may remove this by globally scaling by 0.5?)
                // Acquire stencil
                inds.clear();
                coeffs.clear();
                coeffProvider(isBoundaryRing, isEven, valence, eI, inds, coeffs);
                const int e = edges[eI];
                const double edgeSign = signs[eI];
                const int target = E0ToEk(e, 2 + edgeOrients[eI]);
				maxRow = std::max(maxRow, target + 1);
                for (int i = 0; i < inds.size(); i++)
                {
                    output.emplace_back(target, edges[inds[i]], edgeSign * signs[inds[i]] * coeffs[i]);
                }
            }
        }
    }
}
#endif