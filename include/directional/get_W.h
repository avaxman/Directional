#ifndef DIRECTIONAL_GET_W_H
#define DIRECTIONAL_GET_W_H
#include <igl/igl_inline.h>
#include <Eigen/Sparse>
#include "block_diag.h"
#include <igl/cat.h>
#include "get_U.h"

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
        Eigen::SparseMatrix<double> U;
        get_U(F, EV, FE, EF,1, U);
        using namespace std;
        using namespace Eigen;
        std::vector<Triplet<double>> trips;

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
        GammaToOneform = GammaToOneform * U;
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
        GammaToHalfCurl = GammaToHalfCurl * U;

        igl::cat(1, GammaToOneform, GammaToHalfCurl, W);
    }
    /**
      * Computes the W operator that transforms a gamma2 field to the mean-curl representation.
      *
      */
    IGL_INLINE void get_W(const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& EV,
        const Eigen::MatrixXi& FE,
        const Eigen::MatrixXi& EF,
        const Eigen::VectorXi& matching,
        int N,
        Eigen::SparseMatrix<double>& W)
    {
        // Construct gamma2 to gamma3 converter, easier to construct the oneform and curl converters
        // when usign gamma3 representation.
        Eigen::SparseMatrix<double> UN;
        get_U(F, EV, FE, EF, N, UN);
        using namespace std;
        using namespace Eigen;

        std::vector<Triplet<double>> trips;

        const int edgeCount = EV.rows();
        const int gammaCount = 3 * FE.rows();// Construct half curl operator
        SparseHelper gammaToOneformConstruct(N * edgeCount, N * gammaCount, 2 * N * gammaCount);
        SparseHelper gammaToHalfCurlConstruct(N * edgeCount, N * gammaCount, 2 * N * gammaCount);
        for (int e = 0; e < edgeCount; e++)
        {
            // Find indices of edges in faces
            int ind0 = -1, ind1 = -1;
            for(int i =0; i < 3; ++i)
            {
                if(FE(EF(e,0),i) == e)
                {
                    ind0 = i;
                }
                if (FE(EF(e, 1), i) == e)
                {
                    ind1 = i;
                }
            }
            // Indices of left and right gamma's without levels
            const int gL = 3 * EF(e, 0) + ind0;
            const int gR = 3 * EF(e, 1) + ind1;

            for (int n = 0; n < N; n++)
            {
                gammaToOneformConstruct.addCoeff(e + n * edgeCount, gL + n * gammaCount, 0.5);
                // Follow matching to get to right element
                const int rLevel = modulo(n + matching(e), N);
                gammaToOneformConstruct.addCoeff(e + n * edgeCount, gR + rLevel * gammaCount, 0.5);

                gammaToHalfCurlConstruct.addCoeff(e + n * EF.rows(), gL + n * gammaCount, -1);
                //Right gamma has to be compensated for.
                gammaToHalfCurlConstruct.addCoeff(e + n * EF.rows(), gR + rLevel * gammaCount, 1);
            }
        }
        Eigen::SparseMatrix<double> GammaToOneform = gammaToOneformConstruct.toMat() * UN;
        Eigen::SparseMatrix<double> GammaToHalfCurl = gammaToHalfCurlConstruct.toMat() * UN;

        igl::cat(1, GammaToOneform, GammaToHalfCurl, W);
    }
}
#endif