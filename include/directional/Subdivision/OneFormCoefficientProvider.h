#ifndef DIRECTIONAL_ONEFORMCOEFFICIENTPROVIDER_H
#define DIRECTIONAL_ONEFORMCOEFFICIENTPROVIDER_H
#include <directional/Subdivision/RangeHelpers.h>
struct OneFormCoefficientProvider
{
	using CoeffList = std::vector<double>;
	using IndList = std::vector<int>;
	template<typename T>
	using Vec = std::vector<T, std::allocator<T>>;


	void getEvenBoundaryStencil(int valence, int location, IndList& inds, CoeffList& coeffs)
	{
		const int maxEdgeIndex = 2 * valence - 2;
		if(location == 0 || location == maxEdgeIndex)
		{
			coeffs = { 3. / 8., -1. / 8. };
			inds = { location == 0 ? 0 : maxEdgeIndex, location == 0 ? maxEdgeIndex : 0 };
		}
		else if(location == 2 || location == maxEdgeIndex - 2)
		{
			if(valence == 3)
			{
				Helpers::assign(coeffs, {
					5./3.,
					5./3.,
					26./3.,
					-5./3.,
					5./3.
					},1./32.);
				if (location == 2) inds = { 0,1,2,3,4 };
				else inds = { maxEdgeIndex, maxEdgeIndex-1, maxEdgeIndex-2, maxEdgeIndex-3, maxEdgeIndex-4 };
			}
			else if(valence == 4)
			{
				Helpers::assign(coeffs, {
					5./3.,
					5./3.,
					10,
					-1./3.,
					3.,
					-4./3.,
					-8./3.
					}, 1./32.);
				if (location == 2) Helpers::assignRange(0, 7, inds);
				else inds = { maxEdgeIndex, maxEdgeIndex - 1, maxEdgeIndex - 2, maxEdgeIndex - 3, maxEdgeIndex - 4, maxEdgeIndex - 5, maxEdgeIndex - 6 };
			}
			else
			{
				Helpers::assign(coeffs, {
					5./3.,
					5./3.,
					10.,
					-1./3.,
					3.,
					-4./3.,
					4./3.,
					-4.
					},1./32.);
				if (location == 2) {
					Helpers::assignRange(location - 2, location + 6, inds);
					inds[7] = maxEdgeIndex;
				}
				else 
				{
					Helpers::assignDecreasingRange(maxEdgeIndex - 6, maxEdgeIndex+1, inds);
					inds.push_back(0);
				}
			}
		}
		else
		{
			coeffs = {
				1./32,
				1./32,
				4./32,
				1./32,
				10./32.,
				-1./32.,
				4./32.,
				-1./32.,
				1./32.,
				-4/32., // Boundary coefficient
				-4./32.// Boundary coefficient
			};
			Helpers::assignRange(location - 4, location + 7, inds);
			inds[9] = 0;
			inds[10] = maxEdgeIndex;
		}
	}
	void getOddBoundaryStencil(int valence, int location, IndList& inds, CoeffList& coeffs)
	{
		const int maxEdgeIndex = valence * 2 - 2;
		if (location == 1 || location == maxEdgeIndex - 1)
		{
			if(valence == 2)
			{
				coeffs = { -1. / 4., 1. / 4, 1. / 4. };
				if (location == 1) inds = { 0,1,2 };
				else inds = { maxEdgeIndex, maxEdgeIndex - 1, maxEdgeIndex - 2 };
			}
			else
			{
				coeffs = { -5. / 32., 7. / 32, 6./32, 1./32, 3./32.};
				if (location == 1) inds = { 0, 1, 2, 3, 4 };
				else inds = { maxEdgeIndex, maxEdgeIndex - 1, maxEdgeIndex - 2, maxEdgeIndex - 3, maxEdgeIndex - 4 };
			}
		}
		else
		{
			coeffs = {
				-3./32.,
				1./32,
				-3./32.,
				6./32.,
				3./32.,
				1./32.,
				3./32.
				};
			Helpers::assignRange(location - 3, location + 4, inds);
		}
	}
	void getEvenRegularStencil(int valence, int location, IndList& inds, CoeffList& coeffs)
	{
		// Loop subdivision factor
		const double alfa = LoopCoefficientProvider::factor(valence);
		// Halfbox spline subdivision factor
		const double beta = HalfboxSplineCoefficientProvider::factor(valence);
		// Indices
		IndList localInds;
		Helpers::assignRange(0, 2 * valence, localInds);
		// 
		CircularLookup<int,Vec> ca(&localInds);
		switch(valence)
		{
		case 0:
		case 1:
		case 2:
			throw std::runtime_error("Invalid regular stencil valence");
		case 3:
			coeffs = { 1. / 8. - alfa + beta / 8.,
				beta / 8.,
				3. / 8. - alfa - .25 * beta,
				-beta / 8.,
				1. / 8. - alfa + beta / 8. };
			ca.range(location - 2, location + 3, inds);
			break;
		case 4:
			coeffs = { 1. / 8. - alfa,
				beta / 8.,
				3. / 8. - alfa - .25 * beta, // Target edge
				-beta / 8.,
				1. / 8. - alfa,
				-beta / 8.,
				beta / 4. - alfa };
			ca.range(location - 2, location + 4, inds);
			break;
		default:
			// Higher valence
			break;
		}

	}
	void getOddRegularStencil(int valence, int location, IndList& inds, CoeffList& coeffs)
	{
		IndList localInds;
		Helpers::assignRange(0, valence * 2, localInds);
		CircularLookup<int,Vec> ca(&localInds);

		if(valence == 3)
		{
			coeffs = 
				{
				1./32., 
				-3./32.,
				6./32.,
				3./32.,
				1./32
				};
			ca.range(location - 2, location + 3, inds);
		}
		else
		{
			coeffs =
				{
				-3./32.,
				1. / 32.,
				-3. / 32.,
				6. / 32.,
				3. / 32.,
				1. / 32,
				3. /32.
				};
			ca.range(location - 3, location + 4, inds);
		}
	}
};
#endif