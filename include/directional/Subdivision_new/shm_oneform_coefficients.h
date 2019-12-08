#ifndef DIRECTIONAL_SHM_ONEFORM_COEFFICIENTS_H
#define DIRECTIONAL_SHM_ONEFORM_COEFFICIENTS_H
#include <vector>
#include "hbspline_coefficients.h"
#include "loop_coefficients.h"
namespace directional{
    namespace subdivision{
        void shm_oneform_coefficients(bool isBoundary, bool isEven, int valence, int location, std::vector<int>& inds, std::vector<double>& coeffs){
            const int maxEdgeIndex = 2 * valence - 1 - isBoundary;
            const double z = 3. / 32.;
            const int eCount = 2 * valence - isBoundary;

            // Helper functions
            // Adds the range [min, maxExcl) to the given vector of indices
            auto assignRange = [](std::vector<int>& target, int min, int maxExcl){
                for(int i = min; i < maxExcl;i++)target.push_back(i);
            };
            // Mirrors the indices and coefficients. Scales coefficients appropriately, assuming
            // they are specified as [spoke coeff, ring coeff, ....]
            auto mirror = [&inds, &coeffs,&maxEdgeIndex](){
                for(int i =0; i < inds.size(); i++) inds[i] = maxEdgeIndex-inds[i];
                // std::reverse(inds.begin(), inds.end());
                // Flip coefficients for ring edges due to mirroring of stencil
                for (int i = 1; i < coeffs.size(); i += 2) coeffs[i] *= -1;
            };
            auto makeIndsValid = [&inds, &eCount](){
                for(int i = 0; i < inds.size(); i++) inds[i] = (inds[i] + eCount) % eCount;
            };
            // Sets the coefficients to the provided coefficients, scaled by the given scale
            auto setCoeffs = [&coeffs](const std::vector<double>& coeffsIn, double scale){
                for(double c : coeffsIn) coeffs.push_back(c*scale);
            };

            // The coefficients
            if(isEven){
                if(isBoundary){
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
                        switch(valence){
                        case 2: throw std::runtime_error("Invalid valence on boundary");
                        case 3:
                            coeffs = {
                                2. * z - 3 / 16.,  // = 0
                                2 * z - 3 / 16.,  // = 0
                                3 / 4. - 4. * z,  // = 3/8
                                3 / 16. - 2. * z,  //=0
                                2. * z - 3 / 16. // = 0
                            };
                            assignRange(inds,0, 5);
                            if (isMirrored) {
                                mirror();
                            }
                            break;
                        case 4:
                            coeffs = {
                                z - 1 / 16.,//= 1/32
                                z - 1 / 16., //= 1/32
                                17 / 32. - 2 * z, //= 11/32
                                3 / 32. - z, //= 0
                                z,    //= 3/32
                                -1 / 32.,    
                                -3 / 32.
                            };
                            assignRange(inds,0,7);
                            if(isMirrored) mirror();
                            break;
                        default:
                            coeffs = {
                                z - 1 / 16., //= 1/32
                                z - 1 / 16., //= 1/32
                                17 / 32. - 2 * z, //= 11/32 = Target
                                3 / 32. - z, //= 0
                                z,  //= 3/32
                                -1 / 32., 
                                1 / 32., 
                                -1 / 8. // On boundary
                            };
                            assignRange(inds, 0, 7);

                            if(isMirrored){
                                mirror();
                                // Last element should not be compensated.
                                coeffs.back() *= -1;
                                inds.push_back(0);
                            }
                            else{
                                inds.push_back(maxEdgeIndex);
                            }
                            makeIndsValid();
                            break;
                        }
                    }
                    else
                    {
                        setCoeffs({
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
                        assignRange(inds, location - 4, location + 7);
                        // Set proper indices for boundary elements
                        inds[9] = 0;
                        inds[10] = maxEdgeIndex;
                    }
                }
                else{
                    // Loop subdivision factor
                    const double alfa = loop_factor(valence);
                    // Halfbox spline subdivision factor
                    const double beta = hbspline_factor(valence);

                    const double b8 = beta / 8.;

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
                        assignRange(inds, location - 2, location + 3);
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
                        assignRange(inds, location - 3, location + 5);
                        break;
                    default:
                        coeffs = {};
                        inds = {};
                        // Higher valence
                        coeffs.reserve(9 + valence);
                        inds.reserve(9 + valence);
                        setCoeffs({
                            b8,
                            b8,
                            1. / 8.,
                            b8,
                            3. / 8. - 2. * b8, // Target edge
                            -b8,
                            1. / 8.,
                            -b8,
                            b8
                            },1.0);
                        assignRange(inds, location - 4, location + 5);
                        for(int i = 0; i < valence; i++)
                        {
                            coeffs.push_back(-alfa);
                            inds.push_back(2 * i);
                        }
                        break;
                    }
                    // Make sure all indices are valid
                    makeIndsValid();
                }
            }
            else{
                if(isBoundary){
                    const bool doMirror = location > maxEdgeIndex / 2;
                    // Location relative to closest boundary element in the ring
                    const int genLoc = doMirror ? maxEdgeIndex - location: location;

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
                                z - 1. / 4., // = -5/32
                                z + 1. / 8., // Target = 7/32
                                3. / 8. - 2. * z, //= 3/16
                                1./8. - z, //= 1/32
                                z //= 3/32
                            };
                            assignRange(inds, 0, 5);

                            if(doMirror)
                            {
								for (int i = 0; i < inds.size(); i++) inds[i] = maxEdgeIndex - inds[i];
								for (int i = 0; i < coeffs.size(); i += 2) coeffs[i] *= -1;
                            }
                        }
                    }
                    else
                    {
                        setCoeffs({
                            -3.,
                            1.,
                            -3.,
                            6., // Target
                            3.,
                            1.,
                            3.
                            },1./32.); // Common factor
                        assignRange(inds, location - 3, location + 4);
						makeIndsValid();
                    }
                }
                else{
                    if(valence == 3)
                    {
                        //This is actually a special case of the higher valence one, with th
                        //other coefficients ''folded'' to the same element, which gives a coefficient of 3/32 + -3/32 = 0.
                        setCoeffs({
                            1., 
                            -3.,
                            6.,
                            3.,
                            1.
                            }, 1/32.);
						assignRange(inds, location - 2, location + 3);
                    }
                    else
                    {
                        setCoeffs({
                            -3.,
                            1.,
                            -3.,
                            6., //Center
                            3.,
                            1.,
                            3.
                            }, 1 / 32.);
                        assignRange(inds, location - 3, location + 4);
                    }
                    // Make indices valid. Some may be outside of the vector indices due to circularity.
                    makeIndsValid();
                }
            }
        }
    }
}
#endif