#ifndef DIRECTIONAL_ONEFORM_SUBDIVISION_H
#define DIRECTIONAL_ONEFORM_SUBDIVISION_H
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include "build_subdivision_operators.h"
#include "loop_coefficients.h"
namespace directional{
	namespace subdivision {
        template<typename CoefficientProvider>
        void Se_directional_triplet_provider(
            const int& vertexCount,
            const Eigen::MatrixXi& F0,
            const Eigen::MatrixXi& SFE0,
            const Eigen::MatrixXi& E0,
            const Eigen::MatrixXi& EI0,
            const Eigen::MatrixXi& EF0,
            const Eigen::MatrixXi& E0ToEk,
            const std::vector<int>& edges,
            const std::vector<int>& edgeOrients,
            const Eigen::MatrixXi& edgeLevels, // N/wraps x wraps * valence matrix containing the level for each edge that is part of the same branched cover of the field.
            const CoefficientProvider& coeffProvider,
            int wraps,
            int N,
            std::vector<Eigen::Triplet<double>>& output,
            int& maxRow
        )
        {
            const auto oldEdgeCount = E0.rows();
            const auto newEdgeCount = 2 * E0.rows() + 3 * SFE0.rows();

            // For the time being, we assume that there is no boundary

            // Number of edges
            const int eCount = edges.size();
            // Signs
            std::vector<double> signs(edgeOrients.size(), 1.0);
            for (int i = 0; i < signs.size(); i++) signs[i] -= 2.0 * edgeOrients[i];
            // Vertex valence
            const int valence = eCount / 2;

            // For each ''function'' around the vertex.
            for (int branch = 0; branch < edgeLevels.rows(); ++branch)
            {
                //Local stencil description
                std::vector<int> inds;
                std::vector<double> coeffs;

                // All elements
                for (int eI = 0; eI < eCount; ++eI)
                {
                    // Acquire stencil
                    inds.clear();
                    coeffs.clear();
                    auto isOdd = eI & 1;
                    // Isboundary, isEven, valence, location, output indices, output coefficients 
                    coeffProvider(false, !isOdd, valence, eI, inds, coeffs);
                    const int e = edges[eI];
                    const double edgeSign = signs[eI];
                    const int target = E0ToEk(e, 2 * isOdd + edgeOrients[eI]);
                    maxRow = std::max(maxRow, (int)(target + newEdgeCount * edgeLevels(branch, eI) + 1));
                    for (int i = 0; i < inds.size(); i++)
                    {
                        for (int w = 0; w < wraps; ++w)
                        {
                            output.emplace_back(
                                target + newEdgeCount * edgeLevels(branch, eI + w * eCount),
                                edges[inds[i]] + oldEdgeCount * edgeLevels(branch, inds[i] + w * eCount),
                                edgeSign * signs[inds[i]] * coeffs[i]
                            );
                        }

                    }
                }
            }
        }
    }
}
#endif