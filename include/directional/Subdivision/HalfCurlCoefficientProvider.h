#include "SubdivisionBuilder.h"

struct HalfCurlCoefficientProvider
{
	static void assign(CoefficientList& v, const CoefficientList& data, double factor = 1.0)
	{
		v = CoefficientList(data.size());
		for (int i = 0; i < data.size(); i++) v[i] = data[i] * factor;
	}

	static void assign(IndexList& v, const IndexList& data)
	{
		v = IndexList(data.size());
		for (int i = 0; i < data.size(); i++) v[i] = data[i];
	}
	static void assignRange(int start, int endExcl, IndexList& target)
	{
		for (int i = 0; i < endExcl - start; i++)target[i] = start + i;
	}
	static void assignDecreasingRange(int start, int endExcl, IndexList& target)
	{
		for (int v = endExcl - 1, i = 0; v >= start; v--, i++)target[i] = v;
	}


	void getEvenBoundaryStencil(int valence, int location, IndexList& inds, CoefficientList& coeffs)
	{
		const double z = 3. / 32.;
		const int maxEdgeIndex = valence * 2 - 1;

		if(location == 0 || location == maxEdgeIndex)
		{
			assign(coeffs, { 1. / 4. });
			assign(inds, { location });
		}
		else if(location == 2 || location == maxEdgeIndex - 2)
		{
			if(valence == 3)
			{
				assign(coeffs, { -1, 1, 8, 1, -1 }, 1. / 32.);
				int_range(0, 5, inds);
			}
			else
			{
				assign(coeffs, { z - 5. / 32., z - 3. / 32., 7. / 32., 1. / 8. - z, 3. / 32. - z,1. / 32.,1. / 32 });
				if (location == 2) int_range(0, 7, inds);
				else
				{
					assignDecreasingRange(maxEdgeIndex - 6, maxEdgeIndex + 1, inds);
				}
			}
		}
		else
		{
			assign(coeffs, { 1.,1.,-1.,0, 6., 0, -1., 1., 1. }, 1. / 32.);
			// Range always is within allowed elements
			int_range(location - 4, location + 5, inds);
		}
	}
	void getOddBoundaryStencil(int valence, int location, IndexList& inds, CoefficientList& coeffs)
	{
		IndexList localInds;
		int_range(0, 2 * valence-1, localInds);
		
		CircularAccess<int,Vec> ca(&localInds);

		const double z = 3. / 32.;
		const int maxEdgeIndex = valence * 2 - 2;
		if(location == 1 || location == maxEdgeIndex-1)
		{
			if(valence == 2)
			{
				assign(inds, { 1 });
				assign(coeffs, { 1. / 4. });
			}
			else
			{
				assign(coeffs, { 3. / 32. - z, 9. / 32. - z, 0, z - 3. / 32., z - 1. / 32. });
				if (location == 1) int_range(0, 5, inds);
				else
				{
					assignDecreasingRange(maxEdgeIndex - 4, maxEdgeIndex + 1, inds);
				}
			}
		}
		else
		{
			assign(coeffs, { 1.,2.,1. }, 1. / 16.);
			assign(inds, { ca[location - 3],location, ca[location + 3] });
		}
	}
	void getEvenRegularStencil(int valence, int location, IndexList& inds, CoefficientList& coeffs)
	{
		IndexList localInds;
		int_range(0, 2 * valence, localInds);
		CircularAccess<int,Vec> ca(&localInds);
		switch(valence)
		{
		case 3:
			assign(coeffs, { -1.,1.,11.,1.,-1. }, 1. / 48.);
			ca.range(location - 2, location + 3, inds);
			break;
		case 4:
			assign(coeffs, { -2, 1, 2, 1 ,10, 1, 2, 1 }, 1. / 64.);
			ca.range(location - 4, location + 4, inds);
			break;
		case 5:
		{
			const double fact = 1. / (8. * (2. * std::sqrt(5) + 10));
			assign(coeffs, {
				1. / 32 - fact,
				1. / 32 - fact,
				-1. / 32,
				fact,
				3. / 16. + fact * 2.,
				fact,
				-1. / 32,
				1. / 32 - fact,
				1. / 32 - fact
				});
			ca.range(location - 4, location + 5, inds);
		}
			break;
		case 6:
			assign(coeffs, { 1.,1., 5., 1. }, 1. / 32.);
			assign(inds, {
				ca[location-6],
				ca[location-3],
				location,
				ca[location+3]
				});
			break;
		default:
			assign(coeffs, { 1.,1.,-1, 0, 6, 0, -1, 1, 1 }, 1. / 32.);
			ca.range(location - 4, location + 5, inds);
			break;
		}
	}
	void getOddRegularStencil(int valence, int location, IndexList& inds, CoefficientList& coeffs)
	{
		IndexList localInds;
		int_range(0, 2 * valence, localInds);
		CircularAccess<int,Vec> ca(&localInds);
		if (valence == 3)
		{
			coeffs = { 1. / 8.,1. / 8. };
			inds = { location, ca[location - 3] };
		}
		else
		{
			coeffs = { 1. / 16.,1. / 8., 1. / 16. };
			inds = { location - 3,location,location + 3 };
		}
	}
};