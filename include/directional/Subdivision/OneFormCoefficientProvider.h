#ifndef DIRECTIONAL_ONEFORMCOEFFICIENTPROVIDER_H
#define DIRECTIONAL_ONEFORMCOEFFICIENTPROVIDER_H
#include <directional/Subdivision/RangeHelpers.h>
#include <directional/Subdivision/LoopCoefficientProvider.h>
#include <directional/Subdivision/HalfBoxSplineCoefficientProvider.h>
struct OneFormCoefficientProvider
{
	using CoeffList = std::vector<double>;
	using IndList = std::vector<int>;
	template<typename T>
	using Vec = std::vector<T, std::allocator<T>>;

	// Parameter for boundary stencils
	const double z = 3. / 32.;


	void getEvenBoundaryStencil(int valence, int location, IndList& inds, CoeffList& coeffs)
	{
		const int maxEdgeIndex = 2 * valence - 2;
		// Should fully handle valence = 2 case
		if(location == 0 || location == maxEdgeIndex)
		{
			coeffs = { 3. / 8., -1. / 8. };
			inds = { 0, maxEdgeIndex };
			if (location != 0) std::swap(inds[0], inds[1]);
		}
		else if(location == 2 || location == maxEdgeIndex - 2)
		{
			const int isMirrored = location != 2;
			if (valence == 2) throw std::runtime_error("Invalid valence on boundary");
			if(valence == 3)
			{
				coeffs = {
					2. * z - 3 / 16., 
					2 * z - 3 / 16., 
					3 / 4. - 4. * z, 
					3 / 16. - 2. * z, 
					2. * z - 3 / 16.
				};
				if (!isMirrored) {
					Helpers::intRange(0, 5, inds);
				}
				else {
					Helpers::assignDecreasingRange(maxEdgeIndex - 4, maxEdgeIndex + 1, inds);
					// Flip coefficients for ring edges due to mirroring of stencil
					for (int i = 1; i < coeffs.size(); i += 2) coeffs[i] *= -1;
				}
			}
			else if(valence == 4)
			{
				coeffs = {
					z - 1 / 16.,
					z - 1 / 16., 
					17 / 32. - 2 * z, 
					3 / 32. - z,           
					z,    
					-1 / 32.,    
					-3 / 32.
				};
				if (!isMirrored) {
					Helpers::assignRange(0, 7, inds);
				}
				else {
					// Flip coefficients for ring edges due to mirroring of stencil
					for (int i = 1; i < coeffs.size(); i += 2) coeffs[i] *= -1;
					Helpers::assignDecreasingRange(maxEdgeIndex - 6, maxEdgeIndex + 1, inds);
				}
			}
			else
			{
				coeffs = {
					 z - 1 / 16., 
					z - 1 / 16., 
					17 / 32. - 2 * z, 
					3 / 32. - z, 
					z, 
					-1 / 32., 
					1 / 32., 
					-1 / 8. // On boundary
				};
				if (!isMirrored) {
					Helpers::assignRange(location - 2, location + 6, inds);
					inds[7] = maxEdgeIndex;
				}
				else 
				{
					Helpers::assignDecreasingRange(maxEdgeIndex - 6, maxEdgeIndex+1, inds);
					for (int i = 1; i < coeffs.size(); i += 2) coeffs[i] *= -1;
					// Last element should not be compensated.
					coeffs.back() *= -1;
					inds.push_back(0);
				}
			}
		}
		else
		{
			Helpers::assign(coeffs, {
				1,
				1,
				4,
				1,
				10, // Target
				-1,
				4,
				-1,
				1,
				-4, // Boundary element
				-4 // Boundary element
				}, 1. / 32.);
			Helpers::assignRange(location - 4, location + 7, inds);
			// Set proper indices for boundary elements
			inds[9] = 0;
			inds[10] = maxEdgeIndex;
		}
		DIR_ASSERT(coeffs.size() == inds.size());
	}
	void getOddBoundaryStencil(int valence, int location, IndList& inds, CoeffList& coeffs)
	{
		const int maxEdgeIndex = valence * 2 - 2;

		const bool mirror = location > maxEdgeIndex / 2;
		// Location relative to closest boundary element in the ring
		const int genLoc = mirror ? maxEdgeIndex - location: location;

		if (location == 1 || location == maxEdgeIndex - 1)
		{
			if(valence == 2)
			{
				coeffs = { -1. / 4., 1. / 4., 1. / 4. };
				inds = { 0,1,2 };
			}
			else
			{
				coeffs = {
					z - 1. / 4.,
					z + 1. / 8., // Target
					3. / 8. - 2. * z,
					1./8. - z,
					z
				};
				if (location == 1) {
					Helpers::intRange(0, 5, inds);
				}
				else {
					Helpers::assignDecreasingRange(maxEdgeIndex - 4, maxEdgeIndex + 1, inds);
					// Compensate for mirror: negate all coefficients for target edge orientation,
					// then negate ring coefficients again for being in the opposite direction.
					for (int i = 0; i < coeffs.size(); i += 2) coeffs[i] *= -1;
				}
			}
		}
		else
		{
			Helpers::assign(coeffs,
				{
				-3.,
				1.,
				-3.,
				6., // Target
				3.,
				1.,
				3.
				},1./32.); // Common factor
			Helpers::assignRange(location - 3, location + 4, inds);

			// TODO is this correct?
			//if(mirror) for (int i = 0; i < coeffs.size(); i += 2) coeffs[i] *= -1;
		}
		DIR_ASSERT(coeffs.size() == inds.size());
	}
	void getEvenRegularStencil(int valence, int location, IndList& inds, CoeffList& coeffs)
	{
		// Loop subdivision factor
		const double alfa = LoopCoefficientProvider::factor(valence);
		// Halfbox spline subdivision factor
		const double beta = HalfboxSplineCoefficientProvider::factor(valence);

		const double b8 = beta / 8.;

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
			coeffs = { 
				1. / 8. - alfa + b8,
				b8,
				3. / 8. - alfa - 2.0 * b8,
				-b8,
				1. / 8. - alfa + b8 };
			ca.range(location - 2, location + 3, inds);
			DIR_ASSERT(coeffs.size() == 5);
			break;
		case 4:
			coeffs = { 
				b8,
				1. / 8. - alfa,
				b8,
				3. / 8. - alfa - b8 * 2.0, // Target edge
				-b8,
				1. / 8. - alfa,
				-b8,
				2.0 * b8 - alfa 
			};
			ca.range(location - 3, location + 5, inds);
			DIR_ASSERT(coeffs.size() == 8);
			break;
		default:
			coeffs = {};
			inds = {};
			// Higher valence
			coeffs.reserve(9 + valence);
			inds.reserve(9 + valence);
			Helpers::append(coeffs, {
				b8,
				b8,
				1. / 8.,
				b8,
				3. / 8. - 2. * b8, // Target edge
				-b8,
				1. / 8.,
				-b8,
				b8
				});
			std::vector<int> baseInds;
			ca.range(location - 4, location + 5, baseInds);
			Helpers::append(inds, baseInds);
			for(int i = 0; i < valence; i++)
			{
				coeffs.push_back(-alfa);
				inds.push_back(2 * i);
			}
			break;
		}

		DIR_ASSERT(inds.size() == coeffs.size());
	}
	void getOddRegularStencil(int valence, int location, IndList& inds, CoeffList& coeffs)
	{
		IndList localInds;
		Helpers::assignRange(0, valence * 2, localInds);
		CircularLookup<int,Vec> ca(&localInds);

		if(valence == 3)
		{
			//This is actually a special case of the higher valence one, with th
			//other coefficients ''folded'' to the same element, which gives a coefficient of 3/32 + -3/32 = 0.
			Helpers::assign(coeffs,
				{
				1., 
				-3.,
				6.,
				3.,
				1.
				}, 1/32.);
			ca.range(location - 2, location + 3, inds);
			DIR_ASSERT(inds.size() == 5);
		}
		else
		{
			Helpers::assign(coeffs,
				{
				-3.,
				1.,
				-3.,
				6., //Center
				3.,
				1.,
				3.
				}, 1 / 32.);
			ca.range(location - 3, location + 4, inds);
			DIR_ASSERT(inds.size() == 7);
		}
		DIR_ASSERT(inds.size() == coeffs.size());
	}
};
#endif