
template<typename CoefficientProvider>
class OneFormSubdivider
{
	SubdivisionBuilder builder;
	CoefficientProvider cp;
	EdgeData* ED;
public:
	OneFormSubdivider(EdgeData& ED) :ED(&ED), builder(ED.faceCount(), ED.faceCount()), cp({})
	{
		builder.makeId();
	}
	template<typename...T>
	void prepareNext(Subdivider<T...>& subdivider)
	{

	}

	void handleBoundaryRing(const Eigen::VectorXi& edges, const Eigen::VectorXi& edgeOrients)
	{
		// Number of faces
		const int eCount = edges.size();
		// Valence
		const int valence = (eCount+1)/2; //One extra boundary vert.

		Eigen::VectorXd signs(edgeOrients.size());
		signs = Eigen::VectorXd::Constant(eCount, 1.0) - 2.0 * edgeOrients;

		//Face finder helper. Allows circular index lookup.
		Eigen::VectorXi inds;
		Eigen::VectorXd coeffs;

		// All even elements
		for (int eI : IntRange(0, eCount, 2))
		{
			// Avoid duplicates (may remove this by globally scaling by 0.5?)
			if (edgeOrients(eI) == 1) continue;
			// Acquire stencil
			cp.getEvenBoundaryStencil(valence, eI, inds, coeffs);
			const int e = edges(eI);
			const int edgeSign = signs(eI);
			for(int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(4 * e, edges(inds(i)), sign * signs(inds(i)) * coeffs(i));
			}
		}
		// All odd elements
		for (int eI : IntRange(1, eCount, 2))
		{
			// Avoid duplicates (may remove this by globally scaling by 0.5?)
			if (edgeOrients(eI) == 1) continue;
			// Acquire stencil
			cp.getOddBoundaryStencil(valence, eI, inds, coeffs);
			const int e = edges(eI);
			const int edgeSign = signs(eI);
			for (int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(4 * e, edges(inds(i)), edgeSign * signs(inds(i)) * coeffs(i));
			}
		}
	}
	void handleRegularRing(const Eigen::VectorXi& edges, const Eigen::VectorXi& edgeOrients)
	{
		// Number of faces
		const int eCount = edges.size();
		// Valence
		const int valence = eCount / 2; //One extra boundary vert.

		Eigen::VectorXd signs(edgeOrients.size());
		signs = Eigen::VectorXd::Constant(eCount, 1.0) - 2.0 * edgeOrients;

		//Face finder helper. Allows circular index lookup.
		Eigen::VectorXi inds;
		Eigen::VectorXd coeffs;

		// All even elements
		for (int eI : IntRange(0, eCount, 2))
		{
			// Avoid duplicates (may remove this by globally scaling by 0.5?)
			if (edgeOrients(eI) == 1) continue;
			// Acquire stencil
			cp.getEvenBRegularStencil(valence, eI, inds, coeffs);
			const int e = edges(eI);
			const int edgeSign = signs(eI);
			for (int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(4 * e, edges(inds(i)), edgeSign * signs(inds(i)) * coeffs(i));
			}
		}
		// All odd elements
		for (int eI : IntRange(1, eCount, 2))
		{
			// Avoid duplicates (may remove this by globally scaling by 0.5?)
			if (edgeOrients(eI) == 1) continue;
			// Acquire stencil
			cp.getOddRegularStencil(valence, eI, inds, coeffs);
			const int e = edges(eI);
			const int edgeSign = signs(eI);
			for (int i = 0; i < inds.size(); i++)
			{
				builder.addCoeff(4 * e, edges(inds(i)), edgeSign * signs(inds(i)) * coeffs(i));
			}
		}
	}
	void finalize()
	{
		builder.finalize();
	}
};