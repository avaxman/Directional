#ifndef DIRECTIONAL_GET_W_H
#define DIRECTIONAL_GET_W_H
#include <igl/igl_inline.h>
#include <Eigen/Sparse>
#include "block_diag.h"
#include <igl/cat.h>

namespace directional
{
    /**
      * Computes the W operator that transforms a gamma2 field to the mean-curl representation.
      * 
      */
    IGL_INLINE void get_W(const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& EV,
        const Eigen::MatrixXi& FE,
        const Eigen::MatrixXi& EF,
        Eigen::SparseMatrix<double>& W)
    {
        // Construct gamma2 to gamma3 converter, easier to construct the oneform and curl converters
        // when usign gamma3 representation.
        Eigen::SparseMatrix<double> G2_To_G3(3 * FE.rows(), 2 * FE.rows());
        using namespace Eigen;
        using namespace std;
        std::vector<Triplet<double>> trips;
        trips.reserve(4 * FE.rows());
        for (int f = 0; f < FE.rows(); f++)
        {
            RowVector3d signs;
            for(int c = 0; c < 3; ++c)
            {
                for(int i = 0; i < 3; ++i)
                {
                    if(F(f,i) == EV(FE(f,c),0))
                    {
                        signs(c) = EV(FE(f, c), 1) == F(f, (i + 1) % 3) ? 1.0 : -1.0;
                        break;
                    }
                }
            }


            trips.emplace_back(3 * f, 2 * f, 1.);
            trips.emplace_back(3 * f + 1, 2 * f + 1, 1.);
            // Reconstruct last halfedge form from zero sum constraint.
            trips.emplace_back(3 * f + 2, 2 * f, -signs(0)*signs(2));
            trips.emplace_back(3 * f + 2, 2 * f + 1, -signs(1)*signs(2));
        }
        G2_To_G3.setFromTriplets(trips.begin(), trips.end());

        Eigen::SparseMatrix<double> GammaToOneform(EF.rows(), F.rows()* 3);
        trips.clear();
        trips.reserve(3 * F.rows());
        for (int e = 0; e < EF.rows(); e++)
        {
            const int isBoundary = EF(e, 0) == -1 || EF(e, 1) == -1;
            const double coeff = isBoundary ? 1.0 : 0.5;
            if (EF(e, 0) != -1)
            {
                int ind = -1;
                for (int c = 0; c < 3; ++c) if (e == FE(EF(e, 0), c)) { ind = c; break; }
                trips.emplace_back(e, 3 * EF(e, 0) + ind, coeff);
            }
            if (EF(e, 1) != -1)
            {
                int ind = -1;
                for (int c = 0; c < 3; ++c) if (e == FE(EF(e, 1), c)) { ind = c; break; }
                trips.emplace_back(e, 3 * EF(e, 1) + ind, coeff);
            }
        }
        GammaToOneform.setFromTriplets(trips.begin(), trips.end());
        GammaToOneform = GammaToOneform * G2_To_G3;
        // Construct half curl operator
        Eigen::SparseMatrix<double> GammaToHalfCurl(EF.rows(), F.rows() * 3);
        trips.clear();
        trips.reserve(3 * F.rows());
        for (int e = 0; e < EF.rows(); e++)
        {
            // No curl on boundary
            if (EF(e, 0) == -1 || EF(e, 1) == -1) continue;

            // Right gamma minus left gamma, relative to global edge orientation
            int ind1 = -1, ind0 = -1;
            for (int c = 0; c < 3; ++c) {
                if (e == FE(EF(e, 1), c)) { ind1 = c; }
                if (e == FE(EF(e, 0), c)) { ind0 = c; }
            }
            trips.emplace_back(e, 3 * EF(e, 1) + ind1, 0.5);
            trips.emplace_back(e, 3 * EF(e, 0) + ind0, -0.5);
        }
        GammaToHalfCurl.setFromTriplets(trips.begin(), trips.end());
        GammaToHalfCurl = GammaToHalfCurl * G2_To_G3;

        igl::cat(1, GammaToOneform, GammaToHalfCurl, W);

    }
}
#endif