#ifndef DIRECTIONAL_FACEFIELDSUBDIVIDER_H
#define DIRECTIONAL_FACEFIELDSUBDIVIDER_H
#include <directional/Subdivision/SubdivisionBuilder.h>
struct FaceFinder
{
	template<typename T>
	using Vec = std::vector<T, std::allocator<T>>;
	CircularLookup<int,Vec> edgeLoop;
	CircularLookup<int, Vec> edgeOrient;
	EdgeData* data;
	FaceFinder(const std::vector<int>& edgeLoop, const std::vector<int>& edgeOrient, EdgeData* data) :
		edgeLoop(&edgeLoop),
		edgeOrient(&edgeOrient),
		data(data)
	{}
	// Lookup values via ring edge. Corner gives the corner of the face adjacent to the central vertex.
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
	LeveledSparseConstructor builder;
	CoefficientProvider cp;
	EdgeData* ED;
public:
	FaceFieldSubdivider():cp({}){}

	Eigen::SparseMatrix<double> getMatrix() const
	{
		return builder.matrix;
	}

	void setup(EdgeData& ED)
	{
		this->ED = &ED;
		builder = LeveledSparseConstructor(ED.faceCount(), ED.faceCount());
		builder.makeId();
	}
	template<typename...T>
	void prepareNext(SubdivisionBuilder<T...>& subdivider)
	{
		builder.cols = builder.rows;
		builder.rows = 4 * builder.rows;
	}

	void handleBoundaryRing(const std::vector<int>& edges, const std::vector<int>& edgeOrients)
	{
		// Number of faces
		const int fCount = (edges.size() - 1) / 2;
		// Valence
		const int valence = fCount + 1; //One extra boundary vert.

		//Face finder helper. Allows circular index lookup.
		FaceFinder ffind(edges, edgeOrients, ED);
		std::vector<int> inds;
		std::vector<double> coeffs;

		for (int fI = 0; fI <fCount; fI++)
		{
			const int target = 4 * ffind.face(fI) + ffind.corner(fI);
			cp.getOddBoundaryStencil(valence, fI, inds, coeffs);
			for(int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(target, ffind.face(inds[i]), coeffs[i]);
			}
		}
	}
	void handleRegularRing(const std::vector<int> &edges, const std::vector<int>& edgeOrients)
	{
		// Number of faces
		const int fCount = (edges.size()) / 2;
		// Valence
		const int valence = fCount; //One extra boundary vert.

		//Face finder helper. Allows circular index lookup.
		FaceFinder ffind(edges, edgeOrients, ED);
		std::vector<int> inds;
		std::vector<double> coeffs;

		for (int fI = 0; fI < fCount; fI++)
		{
			const int target = 4 * ffind.face(fI) + ffind.corner(fI);
			cp.getOddRegularStencil(valence, fI, inds, coeffs);
			for (int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(target, ffind.face(inds[i]), coeffs[i]);
			}
		}
	}
	void finalize()
	{
		//Do all even stencils here
		SimpleTimer s;
		s.start();
		// Get triangle triangle adjacency.
		Eigen::MatrixXi tt;
		// Number of adjacent triangles per triangle.
		Eigen::VectorXi counts;
		ED->triangleTriangleAdjacency(tt, counts);

		// Dummy vector
		std::vector<int> inds;

		// Pre-acquire stencils
		std::vector<double> coeffs[4];
		// Boundary stencils
		for(int i = 0; i <= 2; i++) cp.getEvenBoundaryStencil(i, inds, coeffs[i]);
		// Regular stencil
		cp.getEvenRegularStencil(3, 0, inds, coeffs[3]);


		// Per face
		for(int f = 0; f < counts.size(); f++)
		{
			const int target = 4 * f + 3;
			const int count = counts(f);
			// Coefficient of self
			builder.addCoeff(target, f, coeffs[count][0]);
			// Coefficients for others
			for(int i = 0; i < count; i++)
			{
				builder.addCoeff(target, tt(f,i), coeffs[count][i+1]);
			}
		}
		s.stop();
		std::cout << "Facefield finalize coefficients: " << s.elapsed() << " ms " << std::endl;


		builder.finalize();
	}
};
#endif