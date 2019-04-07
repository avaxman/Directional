#ifndef DIRECTIONAL_DCELSUB_H
#define  DIRECTIONAL_DCELSUB_H
#include <Eigen/Eigen>
#include <cassert>
#include <directional/Subdivision/EdgeData.h>

struct DCEL
{
	//Class containing edge data
	EdgeData& data;

	int edge = 0;
	int side = 0;

	DCEL(EdgeData& data) : data(data){}

	inline int corner() const
	{
		return data.EI(edge, side);
	}

	inline int face() const
	{
		return data.EF(edge, side);
	}

	inline bool halfedgeIsOutside()
	{
		return data.EF(edge, side) == -1;
	}

	inline bool twinIsOutside()
	{
		return data.EF(edge, 1 - side) == -1;
	}

	static int mod3(int value)
	{
		static int mod3s[6] = { 0,1,2,0,1,2 };
		return mod3s[value];
	}

	/**
	 * \brief Moves to the next CCW half edge in the current face.
	 */
	void next()
	{
		int f = face();
		int corn = corner();
		//Moves to the next CCW edge in the face
		corn = mod3(corn+1);
		setFaceAndCorner(f, corn);
	}
	
	/**
	 * \brief Moves to the previous CCW half edge in the current face
	 */
	void prev()
	{
		int f = face();
		int corn = corner();
		corn = mod3(corn + 3);
		setFaceAndCorner(f, corn);
	}

	
	/**
	 * \brief Moves to the twin halfedge of the current halfedge. Moving over
	 * the boundary of a mesh is not allowed.
	 */
	void toTwin()
	{
		side = 1 - side;
	}

	int twinFace()
	{
		return data.EF(edge, 1 - side);
	}


	/**
	 * \brief The current edge at which the iterator is
	 * \return The current edge
	 */
	int currentEdge()
	{
		return edge;
	}

	//Returns the edge ID in the face, relative to the current
	//  edge.Offset should be in the[0, 3) integer interval.
	int edgeInFace(int offset  = 0)
	{
		assert(offset >= 0 && offset < 3, "Invalid offset");
		return data.sFE(face(), mod3(corner() + offset));
	}

	int edgeOrientationInFace(int offset = 0)
	{
		assert(offset >= 0 && offset < 3, "Invalid offset");
		return data.sFE(face(), mod3(corner() + offset)+3);
	}


	int startVInFace(int offset = 0)
	{
		assert(offset >= 0 && offset < 3, "Invalid offset");
		if(offset == 0)
			return startVertexId();
		int side = edgeOrientationInFace(offset);
		int e = edgeInFace(offset);
		if(side == 1)
			return data.E(e, 1);
		return data.E(e, 2);
	}
	int endVInFace(int offset = 0) {
		assert(offset >= 0 && offset < 3, "Invalid offset");
		if(offset == 0)
			return endVertexId();
		int side = edgeOrientationInFace(offset);
		int e = edgeInFace(offset);
		if(side == 1)
			return data.E(e, 2);
		return data.E(e, 1);
	}
	// Returns the side of the current edge
	int currentSide(){
		return side;
	}
	// Returns the CCW start of the current halfedge
	int  startVertexId()
	{
		return data.E(edge, side);
	}
	int  endVertexId()
	{
		return data.E(edge, 1 - side);
	}
	// Sets the current edge and side for the iterator
	void setEdge(int edgeId,int side)
	{
		edge = edgeId;
		this->side = side;
	}

	void setFaceAndCorner(int face, int corner)
	{
		edge = data.sFE(face, corner);
		side = data.sFE(face, corner + 3);
	}

	// Returns whether the twin halfedge of the current one is an
	// outer boubdary halfedge
	bool edgeIsBoundary()
	{
		return data.EF(edge, side) == -1 || data.EF(edge, 1 - side) == -1;
	}
	
	template<typename Handler>
	void iterateRings(Handler& h)
	{
		Eigen::MatrixXi boundary;
		data.getBoundaryEdges(boundary);
		assert(boundary.col(1) == Eigen::VectorXi::Constant(boundary.rows(), 0), "Boundary edges are not in canonical position (with boundary at the left)");

		// Ideally, find max valence.
		Eigen::VectorXi VE;
		data.vertexToFirstEdge(VE);

		const int boundaryEdgeCount = boundary.rows();

		//Construct ring per boundary vertex
		for(int eI = 0; eI < boundaryEdgeCount; eI++)
		{
			const int e = boundary(eI, 0);
			std::vector<int> edges;
			std::vector<int> edgeSides;
			// Set edge oriented along outside of boundary
			setEdge(e, 0);
			// Mark handled
			VE(endVertexId()) = -1;
			assert(face() == -1);
			assert(twinFace() != -1);
			do
			{
				toTwin();
				edges.push_back(currentEdge());
				edgeSides.push_back(currentSide());
				next();
				edges.push_back(currentEdge());
				edgeSides.push_back(currentSide());
				next();
			} while (!twinIsOutside());

			// Add the last edge, with the direction pointing away from the central vertex being 0.
			edges.push_back(currentEdge());
			edgeSides.push_back(1 -currentSide());

			h.handleBoundaryRing(edges, edgeSides);
		}

		// Handle regular vertex rings.
		for(int v = 0; v < VE.size(); v++)
		{
			if (VE(v) < 0) continue;
			const int eI = VE(v);
			std::vector<int> edges;
			std::vector<int> edgeSides;
			// Start at the appropriate halfedge
			setEdge(eI, 0);
			if (endVertexId() != v) toTwin();
			assert(face() != -1);
			do
			{
				toTwin();
				edges.push_back(currentEdge());
				edgeSides.push_back(currentSide());
				next();
				edges.push_back(currentEdge());
				edgeSides.push_back(currentSide());
				next();
			} while (currentEdge() != eI);

			h.handleRegularRing(edges, edgeSides);
		}
	}
};
#endif