#include <Eigen/Eigen>
#include <directional/Subdivision/EdgeData.h>
#include <directional/Subdivision/RangeHelpers.h>
struct VertFinder
{
	using IndexList = std::vector<int>;
	using CoeffList = std::vector<double>;
	const IndexList& edgeLoop;
	const IndexList& edgeOrient;
	EdgeData* data;
	VertFinder(const IndexList& edgeLoop, const IndexList& edgeOrient, EdgeData* data) :edgeLoop(edgeLoop),edgeOrient(edgeOrient),data(data){}

	/**
	 * Retrieves the given vertex. 0 is the central vertex, others are the other vertex of the loop (in CCW fashion wrt the normal).
	 * For boundary, the vertex 1 is on the boundary, as well as vertex at position ''valence'' of the centrall vertex.
	 */
	int vertex(int num)
	{
		if(num == 0) return data->E(edgeLoop[num], edgeOrient[num]);
		return data->E(edgeLoop[2 * (num-1)], 1 - edgeOrient[2 * (num - 1)]);
	}
	int boundaryValence() const
	{
		return (edgeLoop.size() + 1) / 2;
	}
	int regularValence() const
	{
		return edgeLoop.size() / 2;
	}
	int vertexEdge(int num) const
	{
		return edgeLoop[2 * (num-1)];
	}
	int vertexEdgeOrient(int num) const
	{
		return edgeOrient[2 * (num-1)];
	}
};

struct LoopCoefficientProvider
{
	using IndexList = std::vector<int>;
	using CoeffList = std::vector<double>;

	CoeffList ODD_COEFFS;

	

	LoopCoefficientProvider()
	{
		ODD_COEFFS = { 3. / 8., 3. / 8., 1. / 8., 1. / 8. };
	}
	static double factor(int valence)
	{
		if (valence == 3) return 3. / 16.;
		return 3. / (8. * valence);
	}
	void getEvenBoundaryStencil(int valence, IndexList& inds, CoeffList& coeffs)
	{
		// Boundary vert - central vert - boundary vert
		inds = { 1, 0, valence };
		coeffs = { 1.0 / 8.0, 3.0 / 4.0, 1.0 / 8.0 };
	}
	void getOddBoundaryStencil(int valence, int location, IndexList& inds, CoeffList& coeffs)
	{
		//On boundary 
		if(location == 1)
		{
			inds = { 0, 1 };
			coeffs = { 0.5, 0.5 };
		}
		// On other side boundary
		else if(location == valence)
		{
			inds = { 0, valence };
			coeffs = { 0.5, 0.5 };
		}
		// "Internal"
		else
		{
			inds = { 0, location, location - 1, location + 1 };
			coeffs = { 3.0 / 8., 3. / 8., 1. / 8.,1. / 8 };
		}
	}
	void getEvenRegularStencil(int valence, IndexList& inds, CoeffList& coeffs)
	{
		const double lFactor = factor(valence);
		coeffs = Helpers::Constant(valence + 1, lFactor);
		coeffs[0] = 1. - valence * lFactor; // Set the central vertex coefficient.

		// Set the target indices.
		Helpers::intRange(0, valence+1, inds);
	}
	void getOddRegularStencil(int valence, int position, IndexList& indices, CoeffList& coeffs)
	{
		// Indices: central vertex, ''target'' vertex, circularly CW vertex, circularly CCW vertex (wrt ''normal'').
		indices = { 0, position, position == 1 ? valence : position - 1, position == valence ? 1 : position + 1 };

		// Always same coefficients
		coeffs = ODD_COEFFS;
	}
};