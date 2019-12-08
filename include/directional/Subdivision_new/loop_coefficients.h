#ifndef DIRECTIONAL_LOOP_COEFFICIENTS_H
#define DIRECTIONAL_LOOP_COEFFICIENTS_H

namespace directional{
namespace subdivision{

	double loop_factor(int valence)
	{
		return valence == 3 ? 3. / 16. : 3. / (8. * valence);
	}

	void loop_coefficients(bool isBoundary, bool isEven, int valence, int location, std::vector<int>& inds, std::vector<double>& coeffs)
	{
		// Loop subdivision factor
		const double factor = loop_factor(valence);

		if(isEven)
		{
			// Boundary stencils
			if(isBoundary)
			{
				inds = { 1, 0, valence };
				coeffs = { 1.0 / 8.0, 3.0 / 4.0, 1.0 / 8.0 };
			}
			else
			{
				const double lFactor = factor;
				coeffs = std::vector<double>(valence+1, factor);
				coeffs[0] = 1. - valence * factor; // Set the central vertex coefficient.

				// Set the target indices.
				inds.resize(valence+1, 0);
				for(int i =0; i < valence+1;i++) inds[i] = i;
			}
		}
		else
		{
			// Boundary stencils
			if(isBoundary)
			{
				//On boundary 
				if (location == 1 || location == valence)
				{
					inds = { 0, location };
					coeffs = { 0.5, 0.5 };
				}
				// "Internal"
				else
				{
					inds = { 0, location, location - 1, location + 1 };
					coeffs = {3./8., 3./8.,1./8.,1./8.};
				}
			}
			else
			{
				// Indices: central vertex, ''target'' vertex, circularly CW vertex, circularly CCW vertex (wrt ''normal'').
				inds = { 0, location, location == 1 ? valence : location - 1, location == valence ? 1 : location + 1 };

				// Always same coefficients
				coeffs = { 3. / 8., 3. / 8., 1. / 8., 1. / 8. };
			}
		}
	}
	
}
}

#endif