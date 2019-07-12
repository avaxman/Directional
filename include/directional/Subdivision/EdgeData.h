#ifndef DIRECTIONAL_EDGEDATA_H
#define DIRECTIONAL_EDGEDATA_H
#include <Eigen/Eigen>
#include <vector>
#include <igl/sort.h>
#include <igl/sortrows.h>
#include <directional/Subdivision/RangeHelpers.h>
#include <directional/dir_assert.h>

namespace directional
{
	namespace intern
	{
		inline void flipEdgeDirection(Eigen::MatrixXi& E,
			Eigen::MatrixXi& EF,
			Eigen::MatrixXi& sFE,
			Eigen::MatrixXi& EI,
			int e)
		{
			std::swap(E(e, 0), E(e, 1));
			if (EF(e, 0) == -1)
			{
				sFE(EF(e, 1), EI(e, 1) + 3) = 1 - sFE(EF(e, 1), EI(e, 1) + 3);
			}
			else if (EF(e, 1) == -1)
			{
				sFE(EF(e, 0), EI(e, 0) + 3) = 1 - sFE(EF(e, 0), EI(e, 0) + 3);
			}
			else
			{
				std::swap(sFE(EF(e, 0), EI(e, 0) + 3), sFE(EF(e, 1), EI(e, 1) + 3));
			}
			std::swap(EF(e, 0), EF(e, 1));
			std::swap(EI(e, 0), EI(e, 1));
		}
	}
	/**
	 * Specific edge topology for subdivision operator constructors for SHC.
	 */
	void construct_edge_topology(const Eigen::MatrixXi& F,
		Eigen::MatrixXi& E,
		Eigen::MatrixXi& EF,
		Eigen::MatrixXi& sFE,
		Eigen::MatrixXi& EI)
	{
		assert(F.cols() == 3, "Only triangle meshes supported");
		sFE.resize(F.rows(), 6);
		Eigen::MatrixXi rawEdges(3 * F.rows(), 2);
		for (int f = 0; f < F.rows(); f++)
		{
			sFE.row(f) = Eigen::RowVectorXi::Constant(6, -1);
			const auto& row = F.row(f);
			rawEdges.row(3 * f) = Eigen::RowVector2i(row(1), row(2));
			rawEdges.row(3 * f + 1) = Eigen::RowVector2i(row(2), row(0));
			rawEdges.row(3 * f + 2) = Eigen::RowVector2i(row(0), row(1));
			// Sort rows ascending. Record swap as opposite orientation.
			for (int i = 0; i < 3; i++) {
				if (rawEdges(3 * f + i, 0) > rawEdges(3 * f + i, 1)) {
					std::swap(rawEdges(3 * f + i, 0), rawEdges(3 * f + i, 1));
					sFE(f, 3 + i) = 1;
				}
				else
				{
					sFE(f, 3 + i) = 0;
				}
			}
		}
		// Sort total rows by first vertex index. Same vertex edges should now be adjacent rows in the matrix
		Eigen::VectorXi edgeMap;
		igl::sortrows(rawEdges, true, E, edgeMap);
		EF = Eigen::MatrixXi::Constant(edgeMap.rows(), 2, -1);
		EI = Eigen::MatrixXi::Constant(edgeMap.rows(), 2, -1);
		// Deduplicate edges and construct other matrices.
		int currId = 0;
		for (int i = 0; i < edgeMap.rows(); i++)
		{
			const int f = edgeMap(i) / 3;
			const int corner = edgeMap(i) - 3 * f;
			// Edge is not same as previous edge.
			if (E.row(i) != E.row(currId))
			{
				currId++;
				// Assign the row to the appropriate row in E.
				E.row(currId) = E.row(i);
			}
			// Update connectivity data.
			sFE(f, corner) = currId;
			EF(currId, sFE(f, 3 + corner)) = f;
			EI(currId, sFE(f, 3 + corner)) = corner;
		}
		currId++;

		// Resize the edge matrices to the appropriate size
		E.conservativeResize(currId, 2);
		EF.conservativeResize(currId, 2);
		EI.conservativeResize(currId, 2);

		//Acquire the boundary edge count
		//boundaryEdgeCount = 0;
		//for (int e = 0; e < currId; e++)
		//{
		//	if (EF(e, 0) == -1 || EF(e, 1) == -1)boundaryEdgeCount++;
		//}

		// Fixup boundary
		for (int e = 0; e < currId; e++)
		{
			if (EF(e, 1) == -1)
			{
				intern::flipEdgeDirection(E, EF, sFE, EI, e);
			}
		}
	}

	bool is_edgetopology_consistent(const Eigen::MatrixXi& E, 
		const Eigen::MatrixXi& F, 
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& sFE, 
		const Eigen::MatrixXi& EI)
	{
		for (int e = 0; e < E.rows(); e++)
		{
			if (EF(e, 0) != -1)
			{
				const int f = EF(e, 0);
				const int c = EI(e, 0);
				// Check edge is at the right corner in the sFE matrix
				if (sFE(f, c) != e) return false;
				// Check edge has the appropriate orientation
				if (sFE(f, c + 3) != 0) return false;
				// Check that the vertices in F correspond to the vertices in E.
				if (F(f, (c + 1) % 3) != E(e, 0) || F(f, (c + 2) % 3) != E(e, 1))
					return false;
			}
			if (EF(e, 1) != -1)
			{
				const int f = EF(e, 1);
				const int c = EI(e, 1);
				if (sFE(f, c) != e) return false;
				if (sFE(f, c + 3) != 1) return false;
				if (F(f, (c + 1) % 3) != E(e, 1) || F(f, (c + 2) % 3) != E(e, 0))
					return false;
			}

		}
		for (int i = 0; i < sFE.rows(); i++)
		{
			for (int j = 0; j < 3; j++)
			{
				const int e = sFE(i,j);
				const int s = sFE(i,j+3);
				// Check all edges are set.
				if (e < 0 || e >= EF.rows()) return false;
				if (EF(e,s) != i ) return false;
			}
		}
		return true;
	}
}

struct EdgeData
{
	Eigen::MatrixXi sFE;
	Eigen::MatrixXi FE_O; // Orientations of edges within face. 0 = CCW, 1 = CW.
	Eigen::MatrixXi EF;
	Eigen::MatrixXi EI;
	Eigen::MatrixXi E;
	Eigen::MatrixXi F;

	int boundaryEdgeCount = 0;
	int vCount = 0;
	int eulerChar = 0;

	void getEdgeTopologyCompatibleFE(Eigen::MatrixXi& FE)
	{
		FE = Eigen::MatrixXi(sFE.rows(), 3);
		FE.col(0) = sFE.col(1);
		FE.col(1) = sFE.col(2);
		FE.col(2) = sFE.col(0);
	}

	bool isConsistent(bool verbose = false) const
	{
		bool fail = false;

		for (int e = 0; e < E.rows(); e++)
		{
			if (EF(e, 0) != -1)
			{
				const int f = EF(e, 0);
				const int c = EI(e, 0);
				if (sFE(f, c) != e) {
					if (verbose) {
						std::cout << "Edge and sFE not the same at e = " << e << ", f,c = (" << f << "," << c << ")" << std::endl;
						fail = true;
					}
					else  return false;
				}
				if (sFE(f, c + 3) != 0) {
					if (verbose) {
						std::cout << "Corrner of sFE not 0 at e = " << e << ", f,c = (" << f << "," << c << ")" << std::endl;
						fail = true;
					}
					else return false;
				}
				if (F(f, (c + 1) % 3) != E(e, 0) && F(f, (c + 2) % 3) != E(e, 1))
				{
					if (verbose) {
						std::cout << "Edge vertices do not concur at e = " << e << ", f,c = (" << f << "," << c << "), E:" << E.row(e) << ", sFE" << sFE.row(f) << std::endl;
						fail = true;
					}
					else return false;
				}
			}
			if (EF(e, 1) != -1)
			{
				const int f = EF(e, 1);
				const int c = EI(e, 1);
				if (sFE(f, c) != e) {
					if (verbose) {
						std::cout << "Edge and sFE not the same at e = " << e << ", f,c = (" << f << "," << c << ")" << std::endl;
						fail = true;
					}
					else return false;
				}
				if (sFE(f, c + 3) != 1) {
					if (verbose) {
						std::cout << "Corrner of sFE not 1 at e = " << e << ", f,c = (" << f << "," << c << ")" << std::endl;
						fail = true;
					}
					else return false;
				}
				if (F(f, (c + 1) % 3) != E(e, 1) && F(f, (c + 2) % 3) != E(e, 0))
				{
					if (verbose) {
						std::cout << "Edge vertices do not concur at e = " << e << ", f,c = (" << f << "," << c << "), E:" << E.row(e) << ", sFE" << sFE.row(f) << std::endl;
						fail = true;
					}
					else return false;
				}
			}

		}
		for(int i = 0; i < sFE.rows(); i++)
		{
			for(int j = 0; j < 3; j++)
			{
				const int e = sFE(i, j);
				const int s = sFE(i, j + 3);
				// Check all edges are set.
				if (e < 0 || e >= EF.rows()) {
					if (verbose) {
						std::cout << "Edge index out of bounds: e = " << e << ", at f  = " << i;
						fail = true;
					}
					else return false;
				}
				if (EF(e, s) != i) {
					if (verbose) {
						std::cout << "Edge flap with e,s = (" << e << ',' << s << ") incorrect face: expected " << i << ", got " << EF(e, s) << std::endl;
						fail = true;
					}
					else return false;
				}
				if (EI(e, s) != j) {
					if (verbose) {
						std::cout << "Edge flap with e,s = (" << e << ',' << s << ") incorrect corner: expected " << j << ", got " << EI(e, s) << std::endl;
						fail = true;
					}
					else return false;
				}
			}
		}
		const int vCount = vertexCount();
		if(eulerChar != vCount - edgeCount() + faceCount())
		{
			if(verbose)
			{
				std::cout << "Invalid vertex count: " << vCount << ", expected" << (eulerChar + edgeCount() - faceCount()) << std::endl;
				fail = true;
			}
			else return false;
		}
		return !fail;
	}

	/**
	 * \brief Create empty EdgeData object.
	 */
	EdgeData() {}

	/**
	 * \brief Construct an edge data object from a face matrix ( |F| x 3 matrix containing vertices in CCW order)
	 * \param F The face matrix
	 */
	EdgeData(const Eigen::MatrixXi& F) : F(F) { construct(); }

	/**
	 * \brief Sets the face matrix of the EdgeData object
	 * \param F The face matrix
	 */
	void setFaceMatrix(const Eigen::MatrixXi& F)
	{
		this->F = F;
	}

	/**
	 * \brief Reconstructs the EdgeData object from the given face matrix
	 * \param F The face matrix
	 */
	void construct(const Eigen::MatrixXi& F)
	{
		this->F = F;
		construct();
	}

	inline void flipEdgeDirection(int e)
	{
		std::swap(E(e, 0), E(e, 1));
		if(EF(e,0) == -1)
		{
			sFE(EF(e, 1), EI(e, 1) + 3) = 1 - sFE(EF(e, 1), EI(e, 1) + 3);
		}
		else if(EF(e,1) == -1)
		{
			sFE(EF(e, 0), EI(e, 0) + 3) = 1 - sFE(EF(e, 0), EI(e, 0) + 3);
		}
		else
		{
			std::swap(sFE(EF(e, 0), EI(e, 0) + 3), sFE(EF(e, 1), EI(e, 1) + 3));
		}
		std::swap(EF(e, 0), EF(e, 1));
		std::swap(EI(e, 0), EI(e, 1));
	}

	/**
	 * \brief Constructs the edge data, provided that an F matrix representing a triangle mesh was previously set.
	 */
	void construct()
	{
		assert(F.cols() == 3, "Only triangle meshes supported");
		sFE.resize(F.rows(), 6);
		Eigen::MatrixXi rawEdges(3 * F.rows(), 2);
		for (int f = 0; f < F.rows(); f++)
		{
			sFE.row(f) = Eigen::RowVectorXi::Constant(6, -1);
			const auto& row = F.row(f);
			rawEdges.row(3 * f) = Eigen::RowVector2i(row(1), row(2));
			rawEdges.row(3 * f + 1) = Eigen::RowVector2i(row(2), row(0));
			rawEdges.row(3 * f + 2) = Eigen::RowVector2i(row(0), row(1));
			// Sort rows ascending. Record swap as opposite orientation.
			for (int i = 0; i < 3; i++) {
				if (rawEdges(3 * f + i, 0) > rawEdges(3 * f + i, 1)) {
					std::swap(rawEdges(3 * f + i, 0), rawEdges(3 * f + i, 1));
					sFE(f, 3 + i) = 1;
				}
				else
				{
					sFE(f, 3 + i) = 0;
				}
			}
		}
		// Sort total rows by first vertex index. Same vertex edges should now be adjacent rows in the matrix
		Eigen::VectorXi edgeMap;
		igl::sortrows(rawEdges, true, E, edgeMap);
		EF = Eigen::MatrixXi::Constant(edgeMap.rows(), 2,-1);
		EI = Eigen::MatrixXi::Constant(edgeMap.rows(), 2, -1);
		// Deduplicate edges and construct other matrices.
		int currId = 0;
		for (int i = 0; i < edgeMap.rows(); i++)
		{
			const int f = edgeMap(i) / 3;
			const int corner = edgeMap(i) - 3 * f;
			// Edge is not same as previous edge.
			if (E.row(i) != E.row(currId))
			{
				currId++;
				// Assign the row to the appropriate row in E.
				E.row(currId) = E.row(i);
			}
			// Update connectivity data.
			sFE(f, corner) = currId;
			EF(currId, sFE(f, 3 + corner)) = f;
			EI(currId, sFE(f, 3 + corner)) = corner;
		}
		currId++;

		// Resize the edge matrices to the appropriate size
		E.conservativeResize(currId, 2);
		EF.conservativeResize(currId, 2);
		EI.conservativeResize(currId, 2);

		//Acquire the boundary edge count
		boundaryEdgeCount = 0;
		for (int e = 0; e < currId; e++)
		{
			if (EF(e, 0) == -1 || EF(e, 1) == -1)boundaryEdgeCount++;
			// Fixup boundary to have outside at left of edge
			if (EF(e, 1) == -1)
			{
				flipEdgeDirection(e);
			}
		}
		
		eulerChar = vertexCount() - edgeCount() + faceCount(); // Minus/Plus boundary loop
	}

	/**
	 * \brief Returns the edge count
	 * \return The number of edges in the mesh
	 */
	int edgeCount() const
	{
		return E.rows();
	}

	/**
	 * \brief Returns the face count
	 * \return The number of faces
	 */
	int faceCount() const
	{
		return sFE.rows();
	}

	int vertexCountFast() const
	{
		return eulerChar + edgeCount() - faceCount();
	}

	/**
	 * \brief Returns the vertex count, deduced from the edge matrix. Cache the result if
	 * performance is important.
	 * \return The number of vertices.
	 */
	int vertexCount() const
	{
		// Max of all edge indices + 1.
		return E.colwise().maxCoeff().maxCoeff() + 1;
	}

	void boundaryLogical(Eigen::VectorXi& target)
	{
		target.resize(edgeCount());
		for (int e = 0; e < edgeCount(); e++)
		{
			target(e) = (EF(e, 0) == -1 || EF(e, 1) == -1);
		}
	}

	static void evenFacesClasses(const Eigen::MatrixXi& faceConnections, Eigen::VectorXi& classes, Eigen::VectorXi& classCount)
	{
		const int fCount = faceConnections.rows();
		classes.resize(fCount);
		classCount.resize(4);
		int inds[4] = { 0,0, fCount - 1,fCount - 1 };
		for (int f = 0; f < fCount; f++)
		{
			int cnt = 0;
			for (int i = 0; i < 4; i++)
			{
				if (faceConnections(f, i) >= 0) cnt++;
			}
			// Update count
			classCount(cnt)++;
			//Place element in sorted order.
			if (cnt == 0)
			{
				classes(inds[1]) = classes(inds[0]);
				inds[1]++;
				classes(inds[0]) = f;
				inds[0]++;
			}
			else if (cnt == 1)
			{
				classes(inds[1]) = f;
				inds[1]++;
			}
			else if (cnt == 2)
			{
				classes(inds[2]) = f;
				inds[2]--;
			}
			else if (cnt == 3)
			{
				classes(inds[2]) = classes(inds[3]);
				inds[2]--;
				classes(inds[3]) = f;
				inds[3]--;
			}
		}
	}

	void triangleTriangleAdjacency(Eigen::MatrixXi& tt, Eigen::VectorXi& counts)
	{
		tt.setConstant(faceCount(), 3, 0);
		counts.resize(faceCount());
		for (int i = 0; i < sFE.rows(); i++)
		{
			int curr = 0;
			for (int j = 0; j < 3; j++)
			{
				const int twinFace = EF(sFE(i, j), 1 - sFE(i, j + 3));
				if (twinFace >= 0)
				{
					tt(i, curr) = twinFace;
					curr++;
				}
			}
			counts(i) = curr;
		}
	}
	void updateEdge(int edge, int side, int face, int corner)
	{
		// Updates an edge in the data object with all related data
		EF(edge, side) = face;
		EI(edge, side) = corner; //Keep corner numbering the same
		sFE(face, corner) = edge;
		sFE(face, corner + 3) = side;
	}
	/**
	 * \brief Rebuilds the F matrix from the current edge data (sFE and E)
			Remember that sFE stores the edge opposite the corner in the
			face. Taking the CCW first vert of every edge gives the
			vertices 1-2-0 per face.
	 */
	void rebuildF()
	{
		F.resize(sFE.rows(), 3);
		for (int f = 0; f < sFE.rows(); f++)
		{
			F(f, 1) = E(sFE(f, 0), sFE(f, 3));
			F(f, 2) = E(sFE(f, 1), sFE(f, 4));
			F(f, 0) = E(sFE(f, 2), sFE(f, 5));
		}
	}

	/**
	 * \brief Reserves space in the matrices of the object to be able to handle
	 * the given number of faces and edges.
	 * \param faceCount Face count
	 * \param edgeCount Edge count
	 */
	void reserveElements(int faceCount, int edgeCount)
	{
		E.setConstant(edgeCount, 2, -1);
		EF.setConstant(edgeCount, 2, -1);
		EI.setConstant(edgeCount, 2, -1);
		sFE.setConstant(faceCount, 6, -1);
	}


	int oppositeFace(int face, int corner)
	{
		int e = sFE(face, corner);
		int side = sFE(face, corner + 3);
		return EF(e, 2 - side);
	}

	/**
	 * \brief Retrieves a list of boundary edges. The list is given as a |BE| x 2 matrix,
	 * with |BE| the number of boundary edges. The first column contains all edge indices
	 * that are boundary. The second column contains 0 or 1 which designates whether the
	 * boundary is to the left of the oriented edge (0) or to the right (1).
	 * \param boundaryEdges Matrix of |BE| x 2 size, describing the boundary edges.
	 */
	void getBoundaryEdges(Eigen::MatrixXi& boundaryEdges)
	{
		// Count number of boundary edges
		int leftCount = 0;
		int rightCount = 0;
		std::vector<int> leftBoundary, rightBoundary;
		for (int i = 0; i < EF.rows(); i++)
		{
			if (EF(i, 0) == -1) leftBoundary.push_back(i);
			if (EF(i, 1) == -1) rightBoundary.push_back(i);
			assert(EF(i, 0) != -1 || EF(i, 1) != -1, "Invalid found while looking for boundary edges. Edge has no faces");
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
	void vertexToFirstEdge(Eigen::VectorXi& VE) const
	{
		VE = Eigen::VectorXi::Constant(vertexCount(), -1);
		for (int e = 0; e < E.rows(); e++)
		{
			if (VE(E(e, 0)) == -1) VE(E(e, 0)) = e;
			else if (VE(E(e, 1)) == -1) VE(E(e, 1)) = e;
		}
	}

	void vertexValence(Eigen::VectorXi& vertexValence) const
	{
		vertexValence.setConstant(vertexCount(), 0);
		for (int e = 0; e < E.rows(); e++)
		{
			vertexValence(E(e, 0)) += 1;
			vertexValence(E(e, 1)) += 1;
		}
	}

	void boundaryVerts(Eigen::VectorXi& BVs) const
	{
		BVs.setConstant(vertexCount(), 0);
		for(int f = 0; f < F.rows();f++)
		{
			for(int j =0; j < 3; j++)
			{
				BVs(F(f, j)) -= 1;
			}
		}
		for(int e = 0; e < E.rows(); e++)
		{
			BVs(E(e, 0)) += 1;
			BVs(E(e, 1)) += 1;
		}
	}

	/**
	 * \brief Quadrisects the edge data in place.
	 * \param vCount The current vertex count
	 * \param E0ToEk Mapping from the original edge to 4 new edges, which are the newly created edges
	 * that are ''parallel'' in the subdivided edge flap.
	 */
	void quadrisect(int vCount, Eigen::MatrixXi& E0ToEk, EdgeData& newData)
	{
		// New number of faces
		const int newFCount = 4 * F.rows();
		const int newECount = 3 * F.rows() + 2 * E.rows();

		const int mod3[6] = { 0,1,2,0,1,2 };

		// Quadrisected edge data with elements reserved.
		newData.reserveElements(newFCount, newECount);

		// Topological invariant
		newData.eulerChar = eulerChar;

		E0ToEk.setConstant(E.rows(), 4, -1);

		// Current edge in new data.
		int currE = 0;

		// Offset for new vertices
		for (int e = 0; e < E.rows(); e++)
		{
			// Get IDs for new edge elements
			const int start = currE;
			const int end = currE + 1;
			const int leftOdd = EF(e, 0) != -1 ? currE + 2 : -1;
			int rightOdd = leftOdd + 1;
			if (leftOdd == -1) rightOdd = currE + 2;
			else if (EF(e, 1) == -1) rightOdd = -1;
			DIR_ASSERT(EF(e, 0) != -1 || EF(e, 1) != -1, "Invalid edge detected, no faces connected");
			currE = std::max(rightOdd, leftOdd) + 1;

			// Update mapping
			E0ToEk.row(e) = Eigen::RowVector4i(start, end, leftOdd, rightOdd);
			
			// Left face not empty
			if (EF(e, 0) != -1)
			{
				const int f = EF(e, 0);
				const int corn = EI(e, 0);
				const int offset = 4 * f;
				newData.updateEdge(start, 0, offset + mod3[corn + 1], corn);
				newData.updateEdge(end, 0, offset + mod3[corn + 2], corn);
				// Furthest face edge
				newData.updateEdge(leftOdd, 0, offset + corn, corn);
				// Center face edge
				newData.updateEdge(leftOdd, 1, offset + 3, corn);

				// Update vertex
				newData.E(leftOdd, 0) = vCount + sFE(f, mod3[corn + 2]);
				newData.E(leftOdd, 1) = vCount + sFE(f, mod3[corn + 1]);
			}
			if (EF(e, 1) != -1)
			{
				const int f = EF(e, 1);
				const int corn = EI(e, 1);
				const int offset = 4 * f;
				newData.updateEdge(start, 1, offset + mod3[corn + 2], corn);
				newData.updateEdge(end, 1, offset + mod3[corn + 1], corn);
				//Furthest face
				newData.updateEdge(rightOdd, 1, offset + corn, corn);
				// Center face
				newData.updateEdge(rightOdd, 0, offset + 3, corn);

				// Update vertex
				newData.E(rightOdd, 0) = vCount + sFE(f, mod3[corn + 1]);
				newData.E(rightOdd, 1) = vCount + sFE(f, mod3[corn + 2]);
			}
			//Vertices for even edges
			newData.E(start, 0) = E(e, 0);
			newData.E(start, 1) = vCount + e; // New odd vertex
			//Vertices for even edges
			newData.E(end, 0) = vCount + e; // New odd vertex
			newData.E(end, 1) = E(e, 1);
		}
		DIR_ASSERT(currE == newECount);
		newData.rebuildF();
		newData.boundaryEdgeCount = boundaryEdgeCount * 2;
	}
};
#endif