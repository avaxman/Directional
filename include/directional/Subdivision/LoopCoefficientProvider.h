#include <Eigen/Eigen>

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
	static Eigen::VectorXd ODD_COEFFS;

	using IndexList = std::vector<int>;
	using CoeffList = std::vector<double>;

	LoopCoefficientProvider()
	{
		ODD_COEFFS = Eigen::VectorXd(4);
		ODD_COEFFS << 3. / 8., 3. / 8., 1. / 8., 1. / 8.;
	}
	static void intRange(int start, int endExcl, Eigen::VectorXi& target)
	{
		target.resize(endExcl - start);
		for (int i = start; i < endExcl; i++) target(i - start) = i;
	}

	static double factor(int valence)
	{
		if (valence == 3) return 3. / 16.;
		return 3. / (8. * valence);
	}
	void getEvenBoundaryStencil(int valence, Eigen::VectorXi& inds, Eigen::VectorXd& coeffs)
	{
		// Boundary vert - central vert - boundary vert
		inds = { 1, 0, valence };
		coeffs = { 1.0 / 8.0, 3.0 / 4.0, 1.0 / 8.0 };
	}
	void getOddBoundaryStencil(int valence, int location, Eigen::VectorXi& inds, Eigen::VectorXd& coeffs)
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
	void getEvenRegularStencil(int valence, Eigen::VectorXi& inds, Eigen::VectorXd& coeffs)
	{
		const double lFactor = factor(valence);
		coeffs = Eigen::VectorXd::Constant(valence + 1, lFactor);
		coeffs(0) = 1. - valence * lFactor; // Set the central vertex coefficient.

		// Set the target indices.
		intRange(0, valence+1, inds);
	}
	void getOddRegularStencil(int valence, int position, Eigen::VectorXi& indices, Eigen::VectorXd& coeffs)
	{
		// Indices: central vertex, ''target'' vertex, circularly CW vertex, circularly CCW vertex (wrt ''normal'').
		indices = { 0, position, position == 1 ? valence : position - 1, position == valence ? 1 : position + 1 };

		// Always same coefficients
		coeffs = ODD_COEFFS;
	}
};