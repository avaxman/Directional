#ifndef DIRECTIONAL_ITERATE_RINGS_H
#define DIRECTIONAL_ITERATE_RINGS_H
#include <directional/dir_assert.h>
#include <Eigen/Eigen>

/**
 * Iterates over all 1-rings in the mesh. Applies the operator() operations of the given handler on
 * vectors of edge IDs and corresponding sides with respect to the canonical orientation for the ring:
 * spoke edges pointing away from the central vertex and ring edges moving counterclockwise around the vertex.
 * The edgeIDs and sides are given as interleaved spoke and ring edges. It follows that boundary vertices can 
 * be dinstinguished by having an odd number of edges, whereas regular vertices have an even number.
 * Input:
 * - sFE: |F| x 6 matrix of face to edge connectivity, along with local orientations. sFE(f,c) with c=0,1,2 gives the edge opposite corner c. sFE(f,c+3) gives the orientation
 *		of the edge relative to the CCW orientation around the face normal: 0 if CCW, 1 if CW.
 * - EF: |E| x 2 matrix of edge to face connectivity, where EF(e,s) with s=0,1 gives the left (s=0) or right(s=1) face connected to edge e. If no face is connected (boundary),
 *		the value is -1.
 * - E: |E| x 2 matrix of edge to vertex connectivity, where edge e points from E(e,0) to E(e,1).
 * - EI: |E| x 2 matrix of edge to face connectivity, where EI(e,s), s=0,1 specifies the corner opposite the edge for the left(s=0) and right(s=1) faces. Again, on boundary,
 *		the value is -1.
 *	- boundaryEdges: |B| x 1 Vector with edge IDs of edges that form the boundary of the mesh.
 *	- vertexCount: Number of vertices in the mesh
 *	- hs: Handler objects. Should define void operator()(const std::vector<int>&, const std::vector<int>&) or something compatible to handle all rings.
 * Output:
 */
template<typename...Handler>
void iterate_rings(
	const Eigen::MatrixXi& sFE,
	const Eigen::MatrixXi& EF,
	const Eigen::MatrixXi& E,
	const Eigen::MatrixXi& EI,
	const Eigen::VectorXi& boundaryEdges,
	const int vertexCount,
	Handler&...hs)
{
    const int boundaryEdgeCount = boundaryEdges.rows();

    DIR_ASSERT_M(boundaryEdgeCount == data.boundaryEdgeCount, "Found different number of boundary edges");

	// Used in hack to make handlers work.
	using expand = int[];

	// Mark vertices when handled
	Eigen::VectorXi seenVerts;
	seenVerts.setConstant(vertexCount, 0);

	// Next function for iterating the ring.
	auto next = [&sFE, &EF](int& edge, int& side)
	{
		// The face we are in.
		const int f = EF(edge, side);
		// New corner
		const int newC = (EI(edge, side) + 1) % 3;
		// Get new (half)edge
		edge = sFE(f, newC);
		// Get orientation of the halfedge
		side = sFE(f, newC + 3);
	};

    //Construct ring per boundary vertex
    for(int eI = 0; eI < boundaryEdgeCount * 2; eI++)
    {
        int e = boundaryEdges(eI/2);
		int s = eI & 1;

		if (seenVerts(E(e, s)) == 1) continue;
		// Make sure that the edge is the start of the counterclockwise ring around the vertex. 
		// If not, another boundary edge will be
		if (EF(e, 1 - s) != -1) continue;

		// Mark seen
		seenVerts(E(e, s)) = 1;

		// Build the loop, composed of spoke and ring edges for the central vertex.
        std::vector<int> edges;
        std::vector<int> edgeSides;

		// Set halfedge oriented along outside of boundary
		if (EF(e, s) != -1) s = 1 - s; //To twin

        // Mark handled
        DIR_ASSERT(EF(e,s) == -1); // Current halfedge is outside
        DIR_ASSERT(EF(e,1-s) != -1); // Twin halfedge is not
        do
        {
			s = 1 - s; // To twin halfedge: switches face.
            const int f = EF(e,s);
            edges.push_back(e);
            edgeSides.push_back(s);
			// Move to next CCW halfedge in the face. This is a ring edge
            next(e,s);
            DIR_ASSERT(f == face());
            edges.push_back(e);
            edgeSides.push_back(s);
			// Move to next CCW halfedge in the face. This is a spoke edge
            next(e, s);
            DIR_ASSERT(f == face());
        } while (EF(e,1-s) != -1); //Loop while the twin halfedge is not outside.

        // Add the last edge, with the direction pointing away from the central vertex being 0.
        edges.push_back(e);
		edgeSides.push_back(1 - s);

		(void)expand {
			(h(edges, edgeSides), 0)...
		};
    }

    // Handle regular vertex rings.
    for(int eI = 0; eI < E.size() * 2; eI++)
    {
		int e = eI / 2;
		int s = eI & 1;
		if (seenVerts(E(e, s)) == 1) continue;

		// Mark seen
		seenVerts(E(e, s)) = 1;

		// Start 
		if (s == 0) s = 1 - s;
        
        std::vector<int> edges;
        std::vector<int> edgeSides;
        
        DIR_ASSERT(EF(e,s) != -1);
        
        int prevFace = EF(e,s);
		const int eInit = e;
        do
        {
			s = 1 - s; // To twin halfedge: switches face.
            DIR_ASSERT(EF(e,1-s) == prevFace);
            const int currF = EF(e,s);
            edges.push_back(e);
            edgeSides.push_back(s);
			// Move to the next CCW halfedge. This is a ring edge.
            next(e, s);
            DIR_ASSERT(currF == EF(e,s));
            edges.push_back(e);
            edgeSides.push_back(s);
			// Move to the next CCW halfedge. This is a spoke edge in opposite direction.
            next(e, s);
            DIR_ASSERT(currF == EF(e,s));
            prevFace = EF(e,s);
        } while (e != eInit); // Continue until we reach our starting edge.

		(void)expand {
			(h(edges, edgeSides), 0)...
		};
    }
}

#endif