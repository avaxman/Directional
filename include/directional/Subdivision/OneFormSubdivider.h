#include <ostream>
template<typename T>
std::ostream& operator<<(std::ostream& str, const std::vector<T>& vec)
{
	for (int i = 0; i < vec.size(); i++) str << vec[i] << ",";
	return str;
}

template<typename CoefficientProvider>
class OneFormSubdivider
{
	LeveledSparseConstructor builder;
	CoefficientProvider cp;
	EdgeData* ED;
	Eigen::MatrixXi* E0ToEk;
	Eigen::VectorXi seenEs;
public:
	OneFormSubdivider():cp({}){}

	Eigen::SparseMatrix<double> getMatrix() const
	{
		return builder.matrix;
	}
	void setup(EdgeData& ED)
	{
		this->ED = &ED;
		builder = LeveledSparseConstructor(ED.edgeCount(), ED.edgeCount());
		builder.makeId();
	}
	template<typename...T>
	void prepareNext(SubdivisionBuilder<T...>& subdivider)
	{
		builder.cols = builder.rows;
		builder.rows = 2 * subdivider.ED.edgeCount() + 3 * subdivider.ED.faceCount();
		E0ToEk = &subdivider.E0ToEk;
		seenEs.setConstant(builder.rows, 0);
	}

	void handleBoundaryRing(const std::vector<int>& edges, const std::vector<int>& edgeOrients)
	{
		// Number of faces
		const int eCount = edges.size();
		// Valence
		const int valence = (eCount+1)/2; //One extra boundary vert.

		std::vector<double> signs(edgeOrients.size());
		signs = Helpers::Constant(eCount, 1.0);
		for(int i = 0; i < signs.size(); i++) signs[i] -= 2.0 * edgeOrients[i];

		//Face finder helper. Allows circular index lookup.
		std::vector<int> inds;  
		std::vector<double> coeffs;

		// All even elements
		for (int eI = 0; eI < eCount; eI +=2)
		{
			inds.clear();
			coeffs.clear();
			// Acquire stencil
			cp.getEvenBoundaryStencil(valence, eI, inds, coeffs);
			const int e = edges[eI];
			const double edgeSign = signs[eI];
			const int target = (*E0ToEk)(e, edgeOrients[eI]);
			DIR_ASSERT(target >= 0);
			DIR_ASSERT(coeffs.size() > 0);
			seenEs(target) += 1;
			if(target == 0)
			{
				std::cout << "Target 0: deg =" << valence << ", pos " << eI << std::endl;

				std::cout << inds << std::endl;
				std::cout << coeffs << std::endl;
			}
			for(int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(target, edges[inds[i]], edgeSign * signs[inds[i]] * coeffs[i]);
			}
		}
		// All odd elements
		for (int eI = 1; eI < eCount; eI += 2)
		{
			inds.clear();
			coeffs.clear();
			// Acquire stencil
			cp.getOddBoundaryStencil(valence, eI, inds, coeffs);
			const int e = edges[eI];
			const double edgeSign = signs[eI];
			const int target = (*E0ToEk)(e, 2 + edgeOrients[eI]);
			DIR_ASSERT(target >= 0);
			DIR_ASSERT(coeffs.size() > 0);
			seenEs(target) += 1;
			for (int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(target, edges[inds[i]], edgeSign * signs[inds[i]] * coeffs[i]);
			}
		}
	}
	void handleRegularRing(const std::vector<int>& edges, const std::vector<int>& edgeOrients)
	{
		// Number of faces
		const int eCount = edges.size();
		// Valence
		const int valence = eCount / 2;

		std::vector<double> signs = Helpers::Constant(eCount, 1.0);
		for (int i = 0; i < signs.size(); i++) signs[i] -= 2.0 * edgeOrients[i];

		//Indices within loop and corresponding coefficients.
		std::vector<int> inds;
		std::vector<double> coeffs;

		// All even elements
		for (int eI = 0 ; eI < eCount; eI += 2)
		{
			// Acquire stencil
			cp.getEvenRegularStencil(valence, eI, inds, coeffs);
			DIR_ASSERT(inds.size() == coeffs.size());
			const int edgeSign = signs[eI];
			const int target = (*E0ToEk)(edges[eI], edgeOrients[eI]);
			seenEs(target) += 1;
			for (int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(target, edges[inds[i]], edgeSign * signs[inds[i]] * coeffs[i]);
			}
		}
		// All odd elements
		for (int eI = 1; eI < eCount; eI += 2)
		{
			// Acquire stencil
			cp.getOddRegularStencil(valence, eI, inds, coeffs);
			DIR_ASSERT(inds.size() == coeffs.size());
			const int target = (*E0ToEk)(edges[eI], 2 + edgeOrients[eI]);
			seenEs(target)+= 1;
			const int edgeSign = signs[eI];
			for (int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(target, edges[inds[i]], edgeSign * signs[inds[i]] * coeffs[i]);
			}
		}
	}
	void finalize()
	{
		builder.finalize();
		for(int i = 0; i < seenEs.rows(); i++)
		{
			DIR_ASSERT(seenEs(i) == 1);
		}
		std::cout << "Fin oneformsub" << std::endl;
	}
};