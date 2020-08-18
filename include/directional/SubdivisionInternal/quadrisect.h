#ifndef DIRECTIONAL_QUADRISECT_H
#define DIRECTIONAL_QUADRISECT_H
#include <Eigen/Eigen>

namespace directional
{
	/**
	 * \brief Quadrisects the edge data in place.
	 * \param vCount The current vertex count
	 * \param E0ToEk Mapping from the original edge to 4 new edges, which are the newly created edges
	 * that are ''parallel'' in the subdivided edge flap.
	 */
	void quadrisect(
		const Eigen::MatrixXi& F0,
		const int& vCount,
		const Eigen::MatrixXi& E0,
		const Eigen::MatrixXi& SFE0,
		const Eigen::MatrixXi& EF0,
		const Eigen::MatrixXi& EI0,
		Eigen::MatrixXi& E0ToEk, 
		Eigen::MatrixXi& F1,
		Eigen::MatrixXi& E1,
		Eigen::MatrixXi& SFE1,
		Eigen::MatrixXi& EF1,
		Eigen::MatrixXi& EI1
	)
	{
		// New number of faces
		const int newFCount = 4 * F0.rows();
		const int newECount = 3 * F0.rows() + 2 * E0.rows();
		
		// Simply modulo 3 operator.
		const int mod3[6] = { 0,1,2,0,1,2 };

		// Reserve elements data
		F1.setConstant(newFCount, 3, -1);
		E1.setConstant(newECount, 2, -1);
		SFE1.setConstant(newFCount, 6, -1);
		E0ToEk.setConstant(E0.rows(), 4, -1);
		EF1.setConstant(newECount, 2, -1);
		EI1.setConstant(newECount, 2, -1);

		// Current edge in new data.
		int currE = 0;

		auto updateEdge = [&EF1,&EI1,&SFE1](int edge, int side, int face, int corner)
		{
			// Updates an edge in the data object with all related data
			EF1(edge, side) = face;
			EI1(edge, side) = corner; //Keep corner numbering the same
			SFE1(face, corner) = edge;
			SFE1(face, corner + 3) = side;
		};

		// Offset for new vertices
		for (int e = 0; e < E0.rows(); e++)
		{
			// Get IDs for new edge elements
			const int start = currE;
			const int end = currE + 1;
			const int leftOdd = EF0(e, 0) != -1 ? currE + 2 : -1;
			int rightOdd = leftOdd + 1;
			if (leftOdd == -1) rightOdd = currE + 2;
			else if (EF0(e, 1) == -1) rightOdd = -1;
			assert(EF0(e, 0) != -1 || EF0(e, 1) != -1, "Invalid edge detected, no faces connected");
			currE = std::max(rightOdd, leftOdd) + 1;

			// Update mapping
			E0ToEk.row(e) = Eigen::RowVector4i(start, end, leftOdd, rightOdd);

			// Left face not empty
			if (EF0(e, 0) != -1)
			{
				const int f = EF0(e, 0);
				const int corn = EI0(e, 0);
				const int offset = 4 * f;
				updateEdge(start, 0, offset + mod3[corn + 1], corn);
				updateEdge(end, 0, offset + mod3[corn + 2], corn);
				// Furthest face edge
				updateEdge(leftOdd, 0, offset + corn, corn);
				// Center face edge
				updateEdge(leftOdd, 1, offset + 3, corn);

				// Update vertex
				E1(leftOdd, 0) = vCount + SFE0(f, mod3[corn + 2]);
				E1(leftOdd, 1) = vCount + SFE0(f, mod3[corn + 1]);
			}
			if (EF0(e, 1) != -1)
			{
				const int f = EF0(e, 1);
				const int corn = EI0(e, 1);
				const int offset = 4 * f;
				updateEdge(start, 1, offset + mod3[corn + 2], corn);
				updateEdge(end, 1, offset + mod3[corn + 1], corn);
				//Furthest face
				updateEdge(rightOdd, 1, offset + corn, corn);
				// Center face
				updateEdge(rightOdd, 0, offset + 3, corn);

				// Update vertex
				E1(rightOdd, 0) = vCount + SFE0(f, mod3[corn + 1]);
				E1(rightOdd, 1) = vCount + SFE0(f, mod3[corn + 2]);
			}
			//Vertices for even edges
			E1(start, 0) = E0(e, 0);
			E1(start, 1) = vCount + e; // New odd vertex
			//Vertices for even edges
			E1(end, 0) = vCount + e; // New odd vertex
			E1(end, 1) = E0(e, 1);
		}
		assert(currE == newECount);

		// Reconstruct the F matrix for the quadrisected mesh
		for (int f = 0; f < SFE1.rows(); f++)
		{
			F1(f, 1) = E1(SFE1(f, 0), SFE1(f, 3));
			F1(f, 2) = E1(SFE1(f, 1), SFE1(f, 4));
			F1(f, 0) = E1(SFE1(f, 2), SFE1(f, 5));
		}
	}
}
#endif