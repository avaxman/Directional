#include "EdgeData.h"
#include "SubdivisionBuilder.h"
template<typename CoefficientProvider>
class EdgeFieldSubdivider
{
	using IndexList = std::vector<int>;
	using CoefficientList = std::vector<double>;
	LeveledSparseConstructor builder;
	CoefficientProvider cp;
	EdgeData* ED;
	// Maps edges to next level edges. |E| x 4 elements, with |E| number of original edges.
	// The values per row e are: 
	// - index for new even element with same start vertex as e
	// - index for new even lement with same end as e
	// - index for odd element in left face, ''parallel'' to e
	// - index for odd element in right face,  ''parallel'' to e
	// If the left or right face does not exist (boundary), the index is -1.
	Eigen::MatrixXi* E0ToEk = nullptr;
public:
	EdgeFieldSubdivider() :cp({}){}

	Eigen::SparseMatrix<double> getMatrix() const
	{
		return builder.matrix;
	}

	void setup(EdgeData& ED)
	{
		builder = LeveledSparseConstructor(ED.edgeCount(), ED.edgeCount());
		builder.makeId();
		this->ED = &ED;
	}
	template<typename...T>
	void prepareNext(SubdivisionBuilder<T...>& subdivider)
	{
		const int previousECount = builder.rows;
		builder.cols = builder.rows;
		builder.rows = ED->faceCount() * 3 + 2 * ED->edgeCount();
		E0ToEk = &subdivider.E0ToEk;
	}

	void handleBoundaryRing(const std::vector<int>& edges, const std::vector<int>& edgeOrients)
	{
		// Number of faces
		const int eCount = edges.size();
		// Valence
		const int valence = (eCount + 1) / 2; //One extra boundary vert.

		//Face finder helper. Allows circular index lookup.
		IndexList inds;
		CoefficientList coeffs;

		// All even elements
		for (int eI = 0; eI < eCount; eI += 2)
		{
			// Acquire stencil
			cp.getEvenBoundaryStencil(valence, eI, inds, coeffs);
			const int e = edges[eI];
			for (int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff((*E0ToEk)(eI, edgeOrients[eI]), edges[inds[i]], coeffs[i]);
			}
		}
		// All odd elements
		for (int eI = 1; eI < eCount; eI += 2)
		{
			// Avoid duplicates (may remove this by globally scaling by 0.5?)
			if (edgeOrients[eI] == 1) continue;
			// Acquire stencil
			cp.getOddBoundaryStencil(valence, eI, inds, coeffs);
			const int e = edges[eI];
			for (int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff((*E0ToEk)(eI, 2 + edgeOrients[eI]), edges[inds[i]], coeffs[i]);
			}
		}
	}
	void handleRegularRing(const std::vector<int>& edges, const std::vector<int>& edgeOrients)
	{
		// Number of faces
		const int eCount = edges.size();
		// Valence
		const int valence = eCount / 2; //One extra boundary vert.

		//Face finder helper. Allows circular index lookup.
		IndexList inds;
		CoefficientList coeffs;

		// All even elements
		for (int eI = 0; eI < eCount; eI += 2)
		{
			// Avoid duplicates (may remove this by globally scaling by 0.5?)
			if (edgeOrients[eI] == 1) continue;
			// Acquire stencil
			cp.getEvenRegularStencil(valence, eI, inds, coeffs);
			const int e = edges[eI];
			const int target = (*E0ToEk)(eI, edgeOrients[eI]);
			for (int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(target, edges[inds[i]], coeffs[i]);
			}
		}
		// All odd elements
		for (int eI = 1; eI < eCount; eI += 2)
		{
			// Avoid duplicates (may remove this by globally scaling by 0.5?)
			if (edgeOrients[eI] == 1) continue;
			// Acquire stencil
			cp.getOddRegularStencil(valence, eI, inds, coeffs);
			const int e = edges[eI];
			const int target = (*E0ToEk)(eI, 2 + edgeOrients[eI]);
			for (int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(target, edges[inds[i]], coeffs[i]);
			}
		}
	}
	void finalize()
	{
		builder.finalize();
	}
};