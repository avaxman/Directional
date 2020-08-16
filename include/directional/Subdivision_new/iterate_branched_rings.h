#ifndef DIRECTIONAL_ITERATE_BRANCHED_RINGS_H
#define DIRECTIONAL_ITERATE_BRANCHED_RINGS_H
#include <Eigen/Eigen>
#include <directional/dir_assert.h>
#include <vector>

namespace directional
{
    template<typename EigenVector>
    inline EigenVector circshift(const EigenVector& inputVector, int amount)
    {
        EigenVector output(inputVector.size());
        if(amount >  0)
        {
            output.tail(inputVector.size() - amount) = inputVector.head(inputVector.size() - amount);
            output.head(amount) = inputVector.tail(amount);
        }
        else
        {
            output.head(inputVector.size() - amount) = inputVector.tail(inputVector.size() - amount);
            output.tail(amount) = inputVector.head(amount);
        }
        return output;
    }

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

        auto originalValence = edges.size() / 2;
        const auto singularValence = wraps * originalValence;
        // Levels for face-based fields
        faceLevels = Eigen::MatrixXi::Zero(N/wraps, singularValence);
        faceLevels.col(0) = Eigen::VectorXi::LinSpaced(N/wraps, 0, N / wraps-1);

        // Make all matchings positive. The exact sign is immaterial to our subdivision, the amount of levels we jump is, i.e.
        // the amount of the circular shift of the directinal indices in the face.
        Eigen::VectorXi positiveMatching(matching);
        for(int i =0;i < positiveMatching.size(); ++i)
        {
            positiveMatching(i) = positiveMatching(i) < 0 ? N + positiveMatching(i) : positiveMatching(i);
        }

        int meshF = 1;
        for(int i = 1; i < singularValence; ++i)
        {
            auto e = (2 * meshF) % edges.size();
            auto match = positiveMatching(e);
            auto levelJump = edgeSigns[e] == 1 ? match : (N - match);
            faceLevels.col(i) = faceLevels.col(i - 1) + Eigen::VectorXi::Ones(faceLevels.rows(), 1) * levelJump;
            meshF = (meshF+1) % originalValence;
        }
        // Apply modulo N to get correct positive levels
        faceLevels = faceLevels.unaryExpr([&N](const int val) {return val % N; });

        // Levels for edge-based fields. Twice singular valence since we have 2 edges per face
        edgeLevels = Eigen::MatrixXi::Zero(N / wraps, 2*singularValence);

        // Do the first wrapping. Then, infer from difference in levels for faceLevels what the different levels should be for the function
        // that are edge-based.
        int lastJump = edgeSigns[0] == 0 ? N - positiveMatching[0] : positiveMatching[0];
        for(auto j = 0; j < originalValence; ++j)
        {
            int jump = edgeSigns[2 * j] == 0 ? 0 : (N - positiveMatching(2 * j));
            edgeLevels(0, 2 * j) = faceLevels(0, j) + jump;

            int jump2 = edgeSigns[2 * j + 1] == 0 ? 0 : (N - positiveMatching(2 * j + 1));
            edgeLevels(0, 2 * j + 1) = faceLevels(0, j) + jump2;
        }
        // Do the other wraps
        Eigen::RowVectorXi wrapOffset = Eigen::RowVectorXi::Ones(1, originalValence) * lastJump;
        for(int w = 1; w < wraps; ++w)
        {
            edgeLevels.block(0, w* originalValence, 1, originalValence) =
                edgeLevels.block(0, (w - 1)* originalValence, 1, originalValence) + wrapOffset;
        }
        // Now copy the unwrapped field by following the difference in face levels
        for(int i = 1; i < N/wraps; ++i)
        {
            edgeLevels.row(i) = edgeLevels.row(i - 1) + Eigen::RowVectorXi::Ones(edgeLevels.cols()) * (faceLevels(i, 0) - faceLevels(i-1, 0));
        }
        // Modulo the levels
        edgeLevels = edgeLevels.unaryExpr([&N](const int val) {return val % N; });
    }

    inline int lcm(int a, int b)
    {
        // GCD implementation without recursion
        auto gcd = [](int a, int b)
        {
            int vals[2] = { a,b };
            while(vals[0] != 0)
            {
                auto temp = vals[0];
                vals[0] = vals[1] % vals[0];
                vals[1] = temp;
            }
            return vals[1];
        };
        return a * b / gcd(a, b);
    }

    /**
	 * \brief 
	 * \tparam Handler 
	 * \param vCount 
	 * \param E 
	 * \param EF 
	 * \param EI 
	 * \param SFE 
	 * \param N 
	 * \param indices 
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
        // We assume here that the mesh has no boundary, though we do not explicitly check yet.

		// Ideally, find max valence.
		Eigen::VectorXi VE;
		VE.setConstant(vCount, -1);
		for(int i = 0; i < E.rows();i++)
		{
			if (VE(E(i, 0)) == -1) VE(E(i, 0)) = i;
			if (VE(E(i, 1)) == -1) VE(E(i, 1)) = i;
		}

		int edge = 0;
        // The current ''level'' of the field: one of the [0,N-1] values that describe the N vectors in the face.
		int side = 0;
		int face = -1;
		int corner = -1;

		auto toTwin = [&edge,&side,&face,&corner,&matching, &EF, &EI] ()
		{
			side = 1 - side;
			face = EF(edge, side);
			corner = EI(edge, side);
		};
		auto next = [&edge, &side, &face, &corner, &SFE]()
		{
			corner = (corner + 1) % 3;
			side = SFE(face, corner + 3);
			edge = SFE(face, corner);
		};

		// Handle regular vertex rings.
		for (int v = 0; v < VE.size(); v++)
		{
			// Already handled this vertex, so ignore
			if (VE(v) < 0) continue;

			// The target edge to start at
			const int eI = VE(v);

			// Mark handled
			VE(v) = -1;

			std::vector<int> edges;
			std::vector<int> edgeSides;
            Eigen::MatrixXi faceLevels,edgeLevels;
			// Start at the appropriate halfedge
			edge = eI; side = 0;

			// Make sure the halfedge we start with points away from the vertex in question
			if (E(edge,1-side) != v) toTwin();

            // Total separate branched fields = N / wraps = N^2 / lcm = N / indices(v) * gcd(N, indices(v))
            do
            {
                toTwin();
                edges.push_back(edge); edgeSides.push_back(side);
                next();
                edges.push_back(edge); edgeSides.push_back(side);
                next();
            } while (edge != eI); //Continue until we made a full loop

            int totalTransport = 0;

            // Compute index for ring
            for(int i = 0; i < edges.size(); i += 2)
            {
                totalTransport += edgeSides[i] == 0 ? N - matching(edges[i]) : matching(edges[i]);
            }
            while (totalTransport < 0) totalTransport += N;
            totalTransport = totalTransport % N;
            // Determine number of times a function wraps around the central vertex
            const int wraps = totalTransport == 0 ? 1 : lcm(N, totalTransport) / N;
            // Determine associated number of branched functions
            const int branches = N / wraps;
            // Get levels for edges and faces
            unwrapField(edges, edgeSides, matching, N, wraps, faceLevels, edgeLevels);
            // Apply handler
            // Invoke visitor once
            h(edges, edgeSides, edgeLevels, wraps, branches);
		}
	}
}
#endif