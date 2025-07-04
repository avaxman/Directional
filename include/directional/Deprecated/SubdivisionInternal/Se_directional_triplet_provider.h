#ifndef DIRECTIONAL_SE_DIRECTIONAL_TRIPLET_PROVIDER_H
#define DIRECTIONAL_SE_DIRECTIONAL_TRIPLET_PROVIDER_H
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
            const Eigen::MatrixXi& faceLevels, 
            const CoefficientProvider& coeffProvider,
            int wraps,
            int N,
            std::vector<Eigen::Triplet<double>>& output,
            int& maxRow
        )
        {
            // Old number of edges in mesh
            const auto oldEdgeCount = E0.rows();
            // Assuming no boundary: new number of edges in mesh
            const auto newEdgeCount = 2 * E0.rows() + 3 * SFE0.rows();

            // For the time being, we assume that there is no boundary

            // Number of edges
            const int eCount = edges.size();
            // Compute relative signs of edges
            std::vector<double> signs(edgeOrients.size(), 1.0);
            for (int i = 0; i < signs.size(); i++) signs[i] -= 2.0 * edgeOrients[i];
            // Maximum row value in the subdivision operator
            maxRow = N * newEdgeCount;

            assert(edgeLevels.cols() == edges.size() * wraps);

            // For each ''function'' around the vertex.
            for (int branch = 0; branch < edgeLevels.rows(); ++branch)
            {
                //Local stencil description
                std::vector<int> inds;
                std::vector<double> coeffs;

                // All elements
                for (int eI_branched = 0; eI_branched < edgeLevels.cols(); ++eI_branched)
                {
                    // Get edge index in non-wrapped ring
                    int eI = eI_branched % eCount;
                    // Acquire stencil
                    inds.clear(); // These are based on the branched function
                    coeffs.clear();
                    auto isOdd = eI & 1;
                    // Arguments: Isboundary, isEven, valence, location, output indices, output coefficients 
                    coeffProvider(false, !isOdd, edgeLevels.cols()/2, eI_branched, inds, coeffs); // Coefficients in branched space
                    const int e = edges[eI];
                    const double edgeSign = signs[eI];
                    const int target = E0ToEk(e, 2 * isOdd + edgeOrients[eI]);
                    for (int i = 0; i < inds.size(); i++)
                    {
                        output.emplace_back(
                            // Using facelevels for odd elements since it is assumed that the odd elements are assigned a matching of 0.
                            target + newEdgeCount * (isOdd ? faceLevels(branch, eI_branched / 2) : edgeLevels(branch, eI_branched)),
                            edges[inds[i] % eCount] + oldEdgeCount * edgeLevels(branch, inds[i]),
                            edgeSign * signs[inds[i] % eCount] * coeffs[i]
                        );
                    }
                }
            }
        }
    }
}
#endif