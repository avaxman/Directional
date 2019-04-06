
struct FaceFinder
{
	CircularAccess edgeLoop;
	CircularAccess edgeOrient;
	EdgeData* data;
	FaceFinder(const Eigen::VectorXi& edgeLoop, const Eigen::VectorXi& edgeOrient, EdgeData* data) :
		edgeLoop(&edgeLoop),
		edgeOrient(&edgeOrient),
		data(data)
	{}
	int face(int ind)
	{
		return data->EF(edgeLoop.at(2 * ind + 1), edgeOrient.at(2 * ind + 1));
	}
	int corner(int ind)
	{
		return data->EI(edgeLoop.at(2 * ind + 1), edgeOrient.at(2 * ind + 1));
	}
};

template<typename CoefficientProvider>
class FaceFieldSubdivider
{
	SubdivisionBuilder builder;
	CoefficientProvider cp;
	EdgeData* ED;
public:
	FaceSubdivider(EdgeData& ED):ED(&ED), builder(ED.faceCount(), ED.faceCount()), cp({})
	{
		builder.makeId();
	}
	template<typename...T>
	void prepareNext(Subdivider<T...>& subdivider)
	{
		builder.cols = builder.rows;
		builder.rows = 4 * builder.rows;
	}

	void handleBoundaryRing(const Eigen::VectorXi& edges, const Eigen::VectorXi& edgeOrients)
	{
		// Number of faces
		const int fCount = (edges.size() - 1) / 2;
		// Valence
		const int valence = fCount + 1; //One extra boundary vert.

		//Face finder helper. Allows circular index lookup.
		FaceFinder ffind(edges, edgeOrients, ED);
		Eigen::VectorXi inds;
		Eigen::VectorXd coeffs;

		for (int fI : IntRange(0, fCount))
		{
			const int target = 4 * ffind.face(fI) + ffind.corner(fI);
			cp.getOddBoundaryStencil(valence, fI, inds, coeffs);
			for(int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(target, ffind.face(inds(i)), coeffs(i));
			}
		}
	}
	void handleRegularRing(const Eigen::VectorXi& edges, const Eigen::VectorXi& edgeOrients)
	{
		// Number of faces
		const int fCount = (edges.size()) / 2;
		// Valence
		const int valence = fCount + 1; //One extra boundary vert.

		//Face finder helper. Allows circular index lookup.
		FaceFinder ffind(edges, edgeOrients, ED);
		Eigen::VectorXi inds;
		Eigen::VectorXd coeffs;

		for (int fI : IntRange(0, fCount))
		{
			const int target = 4 * ffind.face(fI) + ffind.corner(fI);
			cp.getOddRegularStencil(valence, fI, inds, coeffs);
			for (int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(target, ffind.face(inds(i)), coeffs(i));
			}
		}
	}
	void finalize()
	{
		//Do all even stencils here

		// Get triangle triangle adjacency.
		Eigen::MatrixXi tt;
		Eigen::VectorXi counts;
		ED->triangleTriangleAdjacency(tt, counts);

		// Dummy vector
		Eigen::VectorXi inds;

		// Pre-acquire stencils
		Eigen::VectorXd coeffs[5];
		for(int i = 0; i <= 3; i++) cp.getEvenBoundaryStencil(i, inds, coeffs[i]);
		cp.getEvenRegularStencil(4, 0, inds, coeffs[4]);


		// Per face
		for(int f = 0; f < counts.size(); f++)
		{
			const int target = 4 * f + 3;
			const int count = counts(f);
			// Coefficient of self
			builder.addCoeff(target, f, coeffs[count](0));
			// Coefficients for others
			for(int i = 0; i < count; i++)
			{
				builder.addCoeff(target, tt(f,i), coeffs[count](i+1));
			}
		}


		builder.finalize();
	}
};
