#ifndef DIRECTIONAL_HCCOEFFICIENTPROVIDER_H
#define DIRECTIONAL_HCCOEFFICIENTPROVIDER_H
#include <directional/Subdivision/SubdivisionBuilder.h>

struct HalfCurlCoefficientProvider
{
	template<typename T>
	using Vec = std::vector<T, std::allocator<T>>;
	using IndexList = std::vector<int>;
	using CoefficientList = std::vector<double>;
	void getEvenBoundaryStencil(int valence, int location, IndexList& inds, CoefficientList& coeffs)
	{
		const double z = 3. / 32.;
		const int maxEdgeIndex = (valence - 1) * 2;

		if(location == 0 || location == maxEdgeIndex)
		{
			coeffs = { 1. / 4. };
			inds = { location };
		}
		else if(location == 2 || location == maxEdgeIndex - 2)
		{
			if(valence == 3)
			{
				Helpers::assign(coeffs, { -1, 1, 8, 1, -1 }, 1. / 32.);
				Helpers::assignRange(0, 5, inds);
			}
			else
			{
				Helpers::assign(coeffs, { 
					z - 5. / 32.,
					z - 3. / 32., 
					7. / 32., // Target
					1. / 8. - z,
					3. / 32. - z,
					1. / 32.,
					1. / 32. });
				if (location == 2) Helpers::assignRange(0, 7, inds);
				else
				{
					Helpers:: assignDecreasingRange(maxEdgeIndex - 6, maxEdgeIndex + 1, inds);
				}
			}
		}
		else
		{
			Helpers::assign(coeffs, { 1.,1.,-1.,0, 6., 0, -1., 1., 1. }, 1. / 32.);
			// Range always is within allowed elements
			Helpers::assignRange(location - 4, location + 5, inds);
		}
		DIR_ASSERT(coeffs.size() == inds.size());
	}
	void getOddBoundaryStencil(int valence, int location, IndexList& inds, CoefficientList& coeffs)
	{
		const double z = 3. / 32.;
		const int maxEdgeIndex = valence * 2 - 2;
		if(location == 1 || location == maxEdgeIndex-1)
		{
			if(valence == 2)
			{
				inds = { 1 };
				coeffs = { 1. / 4. };
			}
			else
			{
				coeffs = { 3. / 32. - z, 9. / 32. - z, 0, z - 3. / 32., z - 1. / 32. };
				if (location == 1) Helpers::assignRange(0, 5, inds);
				else
				{
					Helpers::assignDecreasingRange(maxEdgeIndex - 4, maxEdgeIndex + 1, inds);
				}
			}
		}
		else
		{
			Helpers::assign(coeffs, { 1.,2.,1. }, 1. / 16.);
			inds = { location - 3,location, location + 3 };
		}
	}
	void getEvenRegularStencil(int valence, int location, IndexList& inds, CoefficientList& coeffs)
	{
		IndexList localInds;
		Helpers::assignRange(0, 2 * valence, localInds);
		CircularLookup<int,Vec> ca(&localInds);
		switch(valence)
		{
		case 3:
			Helpers::assign(coeffs, { -1.,1.,11.,1.,-1., 1. }, 1. / 48.);
			ca.range(location - 2, location + 4, inds);
			break;
		case 4:
			Helpers::assign(coeffs, { -2, 1, 2, 1 ,10, 1, 2, 1 }, 1. / 64.);
			ca.range(location - 4, location + 4, inds);
			break;
		case 5:
		{
			const double fact = 1. / (8. * (2. * std::sqrt(5) + 10.));
			Helpers::assign(coeffs, {
				1. / 32. - fact,
				1. / 32. - fact,
				-1. / 32.,
				fact,
				3. / 16. + fact * 2.,
				fact,
				-1. / 32.,
				1. / 32. - fact,
				1. / 32. - fact
				});
			ca.range(location - 4, location + 5, inds);
		}
			break;
		case 6:
			Helpers::assign(coeffs, { 1.,1., 5., 1. }, 1. / 32.);
			inds = {
				ca[location+6],
				ca[location-3],
				location,
				ca[location+3]
				};
			break;
		default:
			Helpers::assign(coeffs, { 1.,1.,-1, 0, 6, 0, -1, 1, 1 }, 1. / 32.);
			ca.range(location - 4, location + 5, inds);
			break;
		}
		DIR_ASSERT(coeffs.size() == inds.size());
	}
	void getOddRegularStencil(int valence, int location, IndexList& inds, CoefficientList& coeffs)
	{
		IndexList localInds;
		Helpers::assignRange (0, 2 * valence, localInds);
		CircularLookup<int,Vec> ca(&localInds);
		if (valence == 3)
		{
			coeffs = { 1. / 8.,1. / 8. };
			inds = { location, ca[location - 3] };
		}
		else
		{
			coeffs = { 1. / 16.,1. / 8., 1. / 16. };
			inds = { ca[location - 3],location,ca[location + 3] };
		}
	}
};
#endif