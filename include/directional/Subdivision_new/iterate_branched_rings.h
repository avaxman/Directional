#ifndef DIRECTIONAL_ITERATE_BRANCHED_RINGS_H
#define DIRECTIONAL_ITERATE_BRANCHED_RINGS_H
#include <Eigen/Eigen>
#include <directional/dir_assert.h>
#include "directional/DirectionalGamma_Suite.h"
#include <vector>

namespace directional
{

    /**
     * \brief 
     * \param edges The edge IDs of the one-ring
     * \param edgeSigns The signs relative to the local canonical orientation: all outer-ring edges in CCW direction around the vertex
     * and all spoke edges attached to the vertex pointing outwards
     * \param matching The matching for all edges
     * \param N The number of directionals per face
     * \param wraps Number of wrap arounds for each function aroudn this vertex.
     * \param faceLevels wraps x (N/wraps * valence) matrix with the level functions for the unique, face-based functions around the vertex
     * \param edgeLevels Same, but for edge-based functions
     */
    inline void unwrapField(const std::vector<int>& edges, const std::vector<int>& edgeSigns, const Eigen::VectorXi& matching, int N, int wraps, 
        Eigen::MatrixXi& faceLevels,
        Eigen::MatrixXi& edgeLevels)
    {
        // Only support non-boundary rings
        assert((edges.size() & 1) == 0);
        // The mesh valence
        const auto originalValence = edges.size() / 2;
        // Valence of the branched function
        const auto singularValence = wraps * originalValence;
        // Levels for face-based fields
        faceLevels = Eigen::MatrixXi::Zero(N/wraps, singularValence);
        faceLevels.col(0) = Eigen::VectorXi::LinSpaced(N/wraps, 0, N / wraps-1);

        // Make all matchings positive. The exact sign is immaterial to our subdivision, the amount of levels we jump is, i.e.
        // the amount of the circular shift of the directinal indices in the face.
        Eigen::VectorXi positiveMatching(matching.size());
        for(int i =0;i < matching.size(); ++i)
        {
            positiveMatching(i) = modulo(matching(i), N);
        }

        int meshF = 1;
        for(int i = 1; i < singularValence; ++i)
        {
            auto eI = (2 * meshF) % edges.size();
            auto match = positiveMatching(edges[eI]);
            auto levelJump = edgeSigns[eI] == 1 ? match : (N - match);
            faceLevels.col(i) = faceLevels.col(i - 1) + Eigen::VectorXi::Ones(faceLevels.rows(), 1) * levelJump;
            meshF = (meshF+1) % originalValence;
        }
        // The branch function should be complete under the matching
        assert(
            modulo(faceLevels(0, faceLevels.cols() - 1) + (edgeSigns[0] == 1 ? positiveMatching(edges[0]) : (N - positiveMatching(edges[0]))), N)
            != faceLevels(0, 0)
        );

        
        // Apply modulo N to get correct positive levels
        faceLevels = faceLevels.unaryExpr([&N](const int val) {return modulo(val , N); });

        // Levels for edge-based fields. Twice singular valence since we have 2 edges per face
        edgeLevels = Eigen::MatrixXi::Zero(N / wraps, 2*singularValence);

        // Do the first wrapping. Then, infer from difference in levels for faceLevels what the different levels should be for the function
        // that are edge-based.
        Eigen::RowVectorXd edgeLevelOffset(wraps * edges.size());
        for(auto j = 0; j < edgeLevelOffset.size() / 2; ++j)
        {
            // The start edge in the non-wrapped 1-ring.
            int eStart = (2 * j) % edges.size();
            int jump = edgeSigns[eStart] == 0 ? 0 : (N - positiveMatching(edges[eStart]));
            edgeLevelOffset( 2 * j) = jump;

            int jump2 = edgeSigns[eStart + 1] == 0 ? 0 : (N - positiveMatching(edges[eStart + 1]));
            edgeLevelOffset(2* j + 1) = jump2;
        }
        for(int j = 0 ; j < faceLevels.cols(); ++j)
        {
            for(int k =0; k < faceLevels.rows(); ++k)
            {
                edgeLevels(k, 2 * j) = faceLevels(k, j) + edgeLevelOffset(2 * j);
                edgeLevels(k, 2 * j+1) = faceLevels(k, j) + edgeLevelOffset(2 * j+1);
            }
        }
        // Make sure all edge levels are in [0, N)
        edgeLevels = edgeLevels.unaryExpr([&N](const int val) {return modulo(val,N); });
    }

    inline int gcd(int a, int b)
    {
        assert(a >= 0 && b >= 0);
        assert(a != 0 || b != 0);
        // Trivial solutions
        if (a == 0) return b;
        if (b == 0) return a;
        if (a == 1 || b == 1) return 1;

        // Apply GCD computations
        int bv = b;
        int av = a;
        int t = 0;
        while(bv != 0)
        {
            t = bv;
            bv = av % bv;
            av = t;
        }
        return av;
    };

    inline int lcm(int a, int b)
    {
        // GCD implementation without recursion

        return a * b / gcd(a, b);
    }

    /**
	 * \brief Iterates over all 1-rings in the mesh as specified by the input topology, calling the handler with each branched function around the 
	 * vertex as determined by the specified matching
	 * \tparam Handler Type of the visitor. Should be a function/function object that accepts as input
	 * - edges The edges of the rings, given as consecutive spoke and ring edges for the faces of the 1-ring
	 * in CCW order, 
	 *  - the edge signs to make the spoke point outward and the ring edges in CCW direction, given as 0=don't change, 1=negate, 
	 *  - the edge levels of the ring: gives the directional level (branch) for each edge per unique function around the vertex.
	 *  Given as a MatrixXi that has a row for each unique function and the columns contain the level for each of the given edges
	 *  with the given sign
	 *  - the number of wraps a function goes around the vertex
	 *  - the number of branches
	 * \param vCount The number of vertices
	 * \param E 
	 * \param EF 
	 * \param EI 
	 * \param SFE Face to edge connectivity in first 3 columns, and orientations relative to CCW in last 3 columns (per edge), 0 = CCW, 1 = CW.
	 * \param N Number of directionals
	 * \param matching 
	 * \param h A handler for the rings. Has input:
	 *  - edges The edges of the rings, given as consecutive spoke and ring edges for the faces of the 1-ring
	 * in CCW order, 
	 *  - the edge signs to make the spoke point outward and the ring edges in CCW direction, given as 0=don't change, 1=negate, 
	 *  - the edge levels of the ring: gives the directional level (branch) for each edge per unique function around the vertex.
	 *  Given as a MatrixXi that has a row for each unique function and the columns contain the level for each of the given edges
	 *  with the given sign
	 *  - the number of wraps a function goes around the vertex
	 *  - the number of branches
	 */
	template<typename Handler>
	void iterate_branched_rings(
		int vCount,
		const Eigen::MatrixXi& E,
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& EI,
		const Eigen::MatrixXi& SFE,
		int N,
        const Eigen::VectorXi& matching,
		Handler& h)
	{
        auto handleRing = [&matching, &h, &N,&E](const std::vector<int>& edges, const std::vector<int>& edgeSides)
        {
            int totalTransport = 0;

            // Compute index for ring
            for (int i = 0; i < edges.size(); i += 2)
            {
                totalTransport += edgeSides[i] == 0 ? N - matching(edges[i]) : matching(edges[i]);
            }
            Eigen::MatrixXi faceLevels, edgeLevels;
            totalTransport = modulo(totalTransport, N);
            // Determine associated number of branched functions
            const int branches = totalTransport == 0 ? N : gcd(totalTransport, N);
            // Determine number of wrap arounds
            const int wraps = N / branches; 
            // Get levels for edges and faces
            unwrapField(edges, edgeSides, matching, N, wraps, faceLevels, edgeLevels);
            // Apply handler
            h(edges, edgeSides, edgeLevels, faceLevels, wraps, N);
        };
        iterate_rings(vCount, E, EF, EI, SFE, handleRing);
	}
}
#endif