#ifndef DIRECTIONAL_HBSPLINE_COEFFICIENTS_H
#define DIRECTIONAL_HBSPLINE_COEFFICIENTS_H
#include <vector>
#include <igl/PI.h>
#include <math.h>
namespace directional{
    namespace subdivision{
        /**
         * The halfbox spline factor:
         * In:
         * - valence: The vertex valence of the ring in question for which to apply the subdivision
         * Out:
         * - Factor
         */
        inline double hbspline_factor(int valence)
        {
            if (valence == 3) return 1.0 / 12.0;
            if (valence == 4) return 1.0 / 8.0;
            if (valence == 5) return 1.0 / 4.0 - 1. / (16. * std::sin(igl::PI * 2 * 0.2) * std::sin(igl::PI * 2 * 0.2));
            return 0.25;
        }

        inline void hbspline_coefficients(bool isBoundary, bool isEven, int valence, int location, std::vector<int>& inds, std::vector<double>& coeffs){
            // Boundary parameter
            const double z = 3. / 32.;
            const int fCount = isBoundary ? valence-1: valence;
            // Boundary coefficients
            const double A = 5 / 4. - 8 * z; // = 1 / 2
            const double  B = 3 / 2. - 8 * z; // = 3 / 4
            const double  C = 3 / 4.;
            const double  D = 1 / 8.;
            const double  E = 5 / 4. - 4 * z; // = 7/8
            const double  F0 = 4 * z - 1 / 4.; // = 1/8
            const double  F1 = 1 - 4 * z; // = 5/8
            const double  F2 = 1. / 8.;
            if(isEven){
                if(isBoundary){
                    // First coefficient is for original face, others for neighbouring faces.
                    switch(valence)
                    {
                    case 0:
                        coeffs = { 1. / 4. };
                        inds = {0};
                        break;
                    case 1:
                        coeffs = { 2. / 12., 1. / 12.};
                        inds = {0, 1};
                        break;
                    case 2:
                        coeffs = {A/4., (1./2.-A/2.) / 4., (1./2. - A/2.) / 4.};
                        inds = {0,1,2};
                        break;
                    case 3:
                        coeffs = {}; // Not boundary
                        break;
                    default:
                        break;
                    }
                }
                else{
                    // We will handle this separately
                    inds = { 0,1,2,3 };
                    const double f= 1./16.;
                    coeffs = {f, f, f, f};
                }
            }
            else{
                if(isBoundary){
                    const int lastFace = fCount - 1;

                    switch(fCount)
                    {
                    case 1:
                        coeffs = { 1. / 4. };
                        inds = { 0 };
                        break;
                    case 2:
                        coeffs = { E / 4., 1 / 4. - E / 4. }; // 7/32, 1/32
                        if (location == 1)inds = { 1, 0 };
                        else inds = {0, 1};
                        break;
                    default:
                        if(location == 0 || location == lastFace)
                        {
                            coeffs = { C / 4., // 3/16
                                D / 4.,  // 1/32
                                1 / 4. - D / 4. - C / 4. }; // 1/32
                            if (location == 0) inds = {0,1,2};
                            else inds={lastFace,lastFace-1,lastFace-2};
                        }
                        else if(location == 1 || location == lastFace - 1)
                        {
                            if(fCount == 3)
                            {
                                coeffs = { 1 / 8. - B / 8., // 1/32
                                    B / 4., // 3/16
                                    1 / 8. - B / 8. }; // 1/32
                                inds = {0,1,2};
                            }
                            else
                            {
                                coeffs = { 
                                    F0 / 4., // 1/32
                                    F1 / 4., // Target, 5/32
                                    F2 / 4., // 1/32
                                    1 / 4. - F1 / 4. - F2 / 4. - F0 / 4. // 1/32
                                };
                                inds={0,1,2,3};
                                if (location != 1)
                                    for(int i =0; i < 4 ;i++)inds[i]=lastFace-inds[i];
                            }
                        }
                        else
                        {
                            coeffs = {
                                1. / 32.,
                                1./ 32.,
                                4./ 32.,
                                1./ 32.,
                                1./ 32.
                            };
                            inds = {location-2,location-1, location, location+1, location+2};
                        }
                        break;
                    }
                }
                else{
                    const double beta = hbspline_factor(valence);
                    const double b8 = beta / 8.;
                    switch(valence)
                    {
                    case 3:
                        coeffs = {
                            (1. / 8. + beta * 0.5) / 4.,
                            (3. / 4. - beta) / 4.,
                            (1. / 8. + beta * 0.5) / 4.
                        };
                        inds = {location-1, location, location+1};
                        break;
                    case 4:
                        coeffs = {
                            (1. / 8.) /4.,
                            (3. / 4. - beta) /4.,
                            (1. / 8.) /4.,
                            beta / 4.
                        };
                        inds = {location-1, location, location+1, location + 2};
                        break;
                    default:
						coeffs = {
							(beta *0.5) / 4.,
							(1. / 8.) / 4.,
							(3. / 4. - beta) / 4.,
							(1. / 8.) / 4.,
							(beta * 0.5) / 4.
						};
                        inds = {location-2, location-1, location, location+1, location+2};
                        break;
                    }
                    // Make sure no negative indices are present
                    for(int i = 0; i <inds.size();i++) inds[i] =(inds[i]+fCount) % fCount;
                }
            }
        }
    }
}
#endif