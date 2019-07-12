#ifndef DIRECTIONAL_HALFBOXSPLINECOEFFICIENTPROVIDER_H
#define DIRECTIONAL_HALFBOXSPLINECOEFFICIENTPROVIDER_H
#include <igl/PI.h>
#include <directional/Subdivision/RangeHelpers.h>

struct HalfboxSplineCoefficientProvider
{
	// Boundary parameter
	const double z = 3. / 32.;

	// Boundary coefficients
	const double A = 5 / 4. - 8 * z;
	const double  B = 3 / 2. - 8 * z;
	const double  C = 3 / 4.;
	const double  D = 1 / 8.;
	const double  E = 5 / 4. - 4 * z;
	const double  F0 = 4 * z - 1 / 4.;
	const double  F1 = 1 - 4 * z;
	const double  F2 = 1. / 8.;

	template<typename T>
	using Vec = std::vector<T, std::allocator<T>>;

	/**
	 * Half box spline factor. 
	 */
	static double factor(int valence)
	{
		if (valence == 3) return 1.0 / 12.0;
		if (valence == 4) return 1.0 / 8.0;
		if (valence == 5) return 1.0 / 4.0 - 1. / (16. * std::sin(igl::PI * 2 * 0.2) * std::sin(igl::PI * 2 * 0.2));
		return 0.25;
	}
	using IndList = std::vector<int>;
	using CoeffList = std::vector<double>;
	/**
	 * Stencils for faces that are newly created at the center of old faces.
	 */
	void getEvenRegularStencil(int vertexValence, int location, IndList& inds, CoeffList& coeffs )
	{
		// We will handle this separately
		inds = { 0,1,2,3 };
		coeffs = Helpers::Constant(4, 1. / 16.);
	}
	/**
	 * Stencils for faces that move towards original vertex on subdividing.
	 */
	void getOddRegularStencil(int vertexValence, int location, IndList& inds, CoeffList& coeffs)
	{
		const double beta = factor(vertexValence);
		const double b8 = beta / 8.;
		std::vector<int> localInds;
		Helpers::intRange(0, vertexValence, localInds);
		CircularLookup<int, Vec> ca(&localInds);
		switch(vertexValence)
		{
		case 3:
			Helpers::assign(coeffs,
				{
					1. / 8. + beta * 0.5,
					3. / 4. - beta,
					1. / 8. + beta * 0.5
				}, 1. / 4.);
			ca.range(location - 1, location + 2, inds);
			break;
		case 4:
			Helpers::assign(coeffs,
				{
				1. / 8.,
				3. / 4. - beta,
				1. / 8.,
				beta
				}, 1. / 4.);
			ca.range(location - 1, location + 3, inds);
			break;
		default:
			Helpers::assign(coeffs,
				{
				beta *0.5,
				1./8.,
				3. / 4. - beta,
				1. / 8.,
				beta * 0.5
				}, 1. / 4.);
			ca.range(location - 2, location + 3, inds);
			break;
		}
	}
	// Even stencil for new central face.
	void getEvenBoundaryStencil(int neighboursCount,  IndList& inds, CoeffList& coeffs)
	{
		// First coefficient is for original face, others for neighbouring faces.
		switch(neighboursCount)
		{
		case 0:
			coeffs = { 1. / 4. };
			break;
		case 1:
			coeffs = { 2. / 12., 1. / 12.};
			break;
		case 2:
			Helpers::assign(coeffs, { A, 1./2. - A / 2., 1./2. - A / 2. }, 1./4.);
			break;
		case 3:
			coeffs = {}; // Not boundary
			break;
		default:
			break;
		}
	}
	void getOddBoundaryStencil(int vertexValence, int location, IndList& inds, CoeffList& coeffs)
	{
		const int fCount = vertexValence - 1;
		const int lastFace = fCount - 1;

		switch(fCount)
		{
		case 1:
			coeffs = { 1. / 4. };
			inds = { 0 };
			break;
		case 2:
			coeffs = { E / 4., 1 / 4. - E / 4. };
			if (location == 1)inds = { 1, 0 };
			else inds = {0, 1};
			break;
		default:
			if(location == 0 || location == lastFace)
			{
				coeffs = { C / 4., 
					D / 4., 
					1 / 4. - D / 4. - C / 4. };
				if (location == 0) Helpers::assignRange(0, 3, inds);
				else Helpers::assignDecreasingRange(lastFace -2, lastFace +1, inds);
			}
			else if(location == 1 || location == lastFace - 1)
			{
				if(fCount == 3)
				{
					coeffs = { 1 / 8. - B / 8., 
						B / 4., 
						1 / 8. - B / 8. };
					Helpers::assignRange(0, 3, inds);
				}
				else
				{
					coeffs = { 
						F0 / 4., 
						F1 / 4., // Target
						F2 / 4., 
						1 / 4. - F1 / 4. - F2 / 4. - F0 / 4. 
					};
					if (location == 1) Helpers::assignRange(0, 4, inds);
					else Helpers::assignDecreasingRange(lastFace - 3, lastFace + 1, inds);
				}
			}
			else
			{
				Helpers::assign(coeffs, 
					{
					1,
					1,
					4, // Target
					1,
					1,
					}, 1. / 32.);
				Helpers::intRange(location - 2, location + 3, inds);
			}
			break;
		}
	}
};
#endif