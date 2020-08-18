#ifndef DIRECTIONAL_ITERATE_RINGS_H
#define DIRECTIONAL_ITERATE_RINGS_H
#include <Eigen/Eigen>
#include <vector>

namespace directional
{
	/**
	 * \brief Retrieves a list of boundary edges. The list is given as a |BE| x 2 matrix,
	 * with |BE| the number of boundary edges. The first column contains all edge indices
	 * that are boundary. The second column contains 0 or 1 which designates whether the
	 * boundary is to the left of the oriented edge (0) or to the right (1).
	 * \param boundaryEdges Matrix of |BE| x 2 size, describing the boundary edges.
	 */
	inline void boundary_edges(const Eigen::MatrixXi& EF, Eigen::MatrixXi& boundaryEdges)
	{
		// Count number of boundary edges
		int leftCount = 0;
		int rightCount = 0;
		std::vector<int> leftBoundary, rightBoundary;
		for (int i = 0; i < EF.rows(); i++)
		{
			if (EF(i, 0) == -1) leftBoundary.push_back(i);
			if (EF(i, 1) == -1) rightBoundary.push_back(i);
			assert(EF(i, 0) != -1 || EF(i, 1) != -1, "Invalid edge found while looking for boundary edges. Edge has no faces");
		}
		boundaryEdges.resize(leftBoundary.size() + rightBoundary.size(), 2);

		int i = 0;
		for (; i < leftBoundary.size(); i++)
		{
			boundaryEdges(i, 0) = leftBoundary[i];
			boundaryEdges(i, 1) = 0; // Left side is empty.
		}
		const int offset = i;
		for (; i < rightBoundary.size() + offset; i++)
		{
			boundaryEdges(i, 0) = rightBoundary[i - offset];
			boundaryEdges(i, 1) = 1; // Right side is empty.
		}
	}

	template<typename Handler>
	void iterate_rings(
		int vCount,
		const Eigen::MatrixXi& E,
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& EI,
		const Eigen::MatrixXi& SFE,
		Handler& h)
	{
		// Retrieve the boundary edges along with their side
		Eigen::MatrixXi boundary;
		boundary_edges(EF, boundary);
		assert(boundary.col(1) == Eigen::VectorXi::Constant(boundary.rows(), 0), "Boundary edges are not in canonical position (with boundary at the left)");

		// Ideally, find max valence.
		Eigen::VectorXi VE;
		VE.setConstant(vCount, -1);
		for(int i = 0; i < E.rows();i++)
		{
			if (VE(E(i, 0)) == -1) VE(E(i, 0)) = i;
			if (VE(E(i, 1)) == -1) VE(E(i, 1)) = i;
		}

		const auto boundaryEdgeCount = boundary.rows();

		int edge = 0;
		int side = 0;
		int face = -1;
		int corner = -1;

		auto toTwin = [&edge,&side,&face,&corner, &EF, &EI] ()
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



		//DIR_ASSERT_M(boundaryEdgeCount == data.boundaryEdgeCount, "Found different number of boundary edges");

		//Construct ring per boundary vertex
		for (int eI = 0; eI < boundaryEdgeCount; eI++)
		{
			// Get the boundary edge
			const int e = boundary(eI, 0);
			std::vector<int> edges;
			std::vector<int> edgeSides;
			// Set edge oriented along outside of boundary
			edge = e; side = 0; //side = boundary(eI,1);
			// Mark handled
			VE(E(edge,1)) = -1;
			//DIR_ASSERT(face() == -1);
			//DIR_ASSERT(twinFace() != -1);
			do
			{
				toTwin();
				// Add spoke edge
				edges.push_back(edge);
				edgeSides.push_back(side);
				// Move to next edge
				next();
				//DIR_ASSERT(f == face());
				// Add ring edge
				edges.push_back(edge);
				edgeSides.push_back(side);
				next();
				//DIR_ASSERT(f == face());
			} while (EF(edge,1-side)!=-1); // Repeat as long as the twin is not outside the mesh (i.e. the next boundary edge)

			// Add the last edge, with the direction pointing away from the central vertex being 0.
			edges.push_back(edge);
			edgeSides.push_back(1 - side);

            // Invoke visitor
			h(edges, edgeSides);
		}

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
			// Start at the appropriate halfedge
			edge = eI; side = 0;

			// Make sure the halfedge we start with points away from the vertex in question
			if (E(edge,1-side) != v) toTwin();
			//DIR_ASSERT(face() != -1);

			//int prevFace = face();
			do
			{
				toTwin();
				//DIR_ASSERT(twinFace() == prevFace);
				//const int currF = face();
				edges.push_back(edge); edgeSides.push_back(side);
				next();
				//DIR_ASSERT(currF == face());
				edges.push_back(edge); edgeSides.push_back(side);
				next();
				//DIR_ASSERT(currF == face());
				//prevFace = face();
			} while (edge != eI); //Continue until we made a full loop

            // Invoke visitor
			h(edges, edgeSides);
		}
	}
}
#endif