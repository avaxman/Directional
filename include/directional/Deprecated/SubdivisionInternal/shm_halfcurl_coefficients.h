#ifndef DIRECTIONAL_SHM_HALFCURL_COEFFICIENTS_H
#define DIRECTIONAL_SHM_HALFCURL_COEFFICIENTS_H
#include <vector>
namespace directional{
    namespace subdivision{
        void shm_halfcurl_coefficients(bool isBoundary, bool isEven, int valence, int location, std::vector<int>& inds, std::vector<double>& coeffs){
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
            };
            auto makeIndsValid = [&inds, &eCount](){
                for(int i = 0; i < inds.size(); i++) inds[i] = (inds[i] + eCount) % eCount;
            };
            // Sets the coefficients to the provided coefficients, scaled by the given scale
            auto setCoeffs = [&coeffs](const std::vector<double>& coeffsIn, double scale){
                for(double c : coeffsIn) coeffs.push_back(c*scale);
            };


            if(isEven){
                if(isBoundary) 	{
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
                            setCoeffs({ -1, 1, 8, 1, -1 }, 1. / 32.);
                            assignRange(inds, 0, 5);
                        }
                        else
                        {
                            setCoeffs({ 
                                z - 5. / 32.,//=-1/16
                                z - 3. / 32., //=0
                                7. / 32., // Target
                                1. / 8. - z, //=1/32
                                3. / 32. - z, //=0
                                1. / 32.,
                                1. / 32. },1.0);
                            assignRange(inds, 0, 7);
                            if(location != 2) mirror();
                        }
                    }
                    else
                    {
                        setCoeffs({ 1.,1.,-1.,0, 6., 0, -1., 1., 1. }, 1. / 32.);
                        // Range always is within allowed elements
                        assignRange(inds, location - 4, location + 5);
                    }
                }
                else{
                    switch(valence)
                    {
                    case 3:
                        setCoeffs({ -1.,1.,11.,1.,-1., 1. }, 1. / 48.);
                        assignRange(inds, location - 2, location + 4);
                        break;
                    case 4:
                        setCoeffs({ -2, 1, 2, 1 ,10, 1, 2, 1 }, 1. / 64.);
                        assignRange(inds, location - 4, location + 4);
                        break;
                    case 5:
                    {
                        const double fact = 1. / (16. * (std::sqrt(5) + 5.));//1. / (8. * (2. * std::sqrt(5) + 10.));
                            // zeta = 2 / (sqrt(5) + 5) -> fact = 1/32 * zeta = 1/(16 * (sqrt(5) + 5))
                        setCoeffs({
                            1. / 32. - fact,
                            1. / 32. - fact,
                            -1. / 32.,
                            fact,
                            3. / 16. + fact * 2.,
                            fact,
                            -1. / 32.,
                            1. / 32. - fact,
                            1. / 32. - fact
                            },1.0);
                        assignRange(inds, location - 4, location + 5);
                    }
                        break;
                    case 6:
                        setCoeffs({ 1.,1., 5., 1. }, 1. / 32.);
                        inds = {
                            location+6,
                            location-3,
                            location,
                            location+3
                            };
                        break;
                    default:
                        setCoeffs({ 1.,1.,-1, 0, 6, 0, -1, 1, 1 }, 1. / 32.);
                        assignRange(inds, location - 4, location + 5);
                        break;
                    }
                    makeIndsValid();
                }
            }
            else{
                if(isBoundary){   
                    if(location == 1 || location == maxEdgeIndex-1)
                    {
                        if(valence == 2)
                        {
                            inds = { 1 };
                            coeffs = { 1. / 4. };
                        }
                        else
                        {
                            coeffs = { 
                                3. / 32. - z, //=0
                                9. / 32. - z, //=3/16
                                0, 
                                z - 3. / 32., //= 0
                                z - 1. / 32. //= 1/16
                            };
                            assignRange(inds, 0, 5);
                            if (location != 1) mirror();
                        }
                    }
                    else
                    {
                        setCoeffs({ 1.,2.,1. }, 1. / 16.);
                        inds = { location - 3,location, location + 3 };
                    }
                }
                else{
                    if (valence == 3)
                    {
                        coeffs = { 1. / 8.,1. / 8. };
                        inds = { location, location - 3};
                    }
                    else
                    {
                        coeffs = { 1. / 16.,1. / 8., 1. / 16. };
                        inds = { location - 3,location,location + 3 };
                    }
                    makeIndsValid();
                }
            }
        }
    }
}
#endif