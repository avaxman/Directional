#ifndef DIRECTIONAL_HALFBOXSPLINECOEFFICIENTPROVIDER_H
#define DIRECTIONAL_HALFBOXSPLINECOEFFICIENTPROVIDER_H
#include <igl/PI.h>
#include <directional/Subdivision/RangeHelpers.h>

struct HalfboxSplineCoefficientProvider
{
	static int factor(int valence)
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
		const double beta = factor(3);
		switch(vertexValence)
		{
		case 3:
			coeffs = Helpers::Constant(3, 1./32. + beta * 1./8.);
			coeffs[location] = 3. / 16. - beta / 4.;
			inds = { 0,1,2 };
			break;
		case 4:
			coeffs = Helpers::Constant(3, 1. / 32.);
			coeffs[location] = 3. / 16. - beta / 4.;
			coeffs[location < 2 ? location + 2 : location - 2] = beta / 4.;
			inds = { 0,1,2,3 };
			break;
		case 5:
			coeffs = Helpers::Constant(5,beta / 8.);
			coeffs[location] = 3. / 16. - beta / 4.;
			coeffs[location == 0 ? vertexValence - 1 : location - 1] = 1. / 32.;
			coeffs[location == vertexValence - 1 ? 0 : location + 1] = 1. / 32.;
			Helpers::assignRange(0, vertexValence, inds);
			break;
		default:
			coeffs = { beta * 0.5, 1. / 8., 3. / 4. - beta, 1. / 8., beta * 0.5 };
			Helpers::assignWrappedRange(location - 2, location + 3, vertexValence, inds);
		}
	}
	void getEvenBoundaryStencil(int neighboursCount,  IndList& inds, CoeffList& coeffs)
	{
		switch(neighboursCount)
		{
		case 0:
			coeffs = { 1. / 4. };
			break;
		case 1:
			coeffs = { 5. / 16., -1. / 16. };
			break;
		case 2:
			coeffs = { 6. / 32., 1. / 32.,1. / 32. };
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
		switch(vertexValence)
		{
		case 2:
			coeffs = { 1. / 4. };
			inds = { 0 };
			break;
		case 3:
			coeffs = { 2. / 12., 1. / 12. };
			inds = { location, 1 - location };
			break;
		case 4:
			inds = { 0,1,2 };
			coeffs = Helpers::Constant(3, 1. / 24.);
			coeffs[location] = 4. / 24.;
			break;
		default:
			if(location == 0 || location == vertexValence-2)
			{
				coeffs = { 4. / 24., 1. / 24., 1. / 24. };
				if (location == 0) inds = { 0,1,2 };
				else inds = { location, location - 1, location - 2 };
			}
			else if(location == 1 || location == vertexValence - 3)
			{
				coeffs = { 5. / 96., 14. / 96., 2. / 96., 3. / 96. };
				if(location == 1)
				{
					inds = { 0,1,2,3 };
				}
				else inds = { location + 1, location, location - 1, location - 2 };
			}
			else
			{
				Helpers::assignWrappedRange(location - 2, location + 3, vertexValence, inds);
				coeffs = { 1. / 32., 1. / 32., 4. / 32., 1. / .32, 1. / 32. };
			}
			break;
		}
	}
};
#endif