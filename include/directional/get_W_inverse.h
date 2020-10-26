#ifndef DIRECTIONAL_GET_W_INVERSE_H
#define DIRECTIONAL_GET_W_INVERSE_H
#include <igl/igl_inline.h>
#include <Eigen/Sparse>
#include "block_diag.h"
#include <igl/cat.h>
#include "get_U_inverse.h"

namespace directional
{
    IGL_INLINE void get_W_inverse(const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& EV,
        const Eigen::MatrixXi& FE,
        const Eigen::MatrixXi& EF,
        Eigen::SparseMatrix<double>& W_inverse)
    {

        using SMat = Eigen::SparseMatrix<double>;
        // Matrix to convert oneform to gamma3 elements
        SMat OneFormToG3(3 * FE.rows(), EF.rows()), CToG3(3 * FE.rows(), EF.rows());

        std::vector<Eigen::Triplet<double>> tripsCBack, tripsOneFormBack;
        tripsCBack.reserve(2 * EF.rows());
        tripsOneFormBack.reserve(2 * EF.rows());
        for (int e = 0; e < EF.rows(); e++)
        {
            int ind0 = -1, ind1 = -1;
            for(int c = 0; c < 3; ++c)
            {
                if (EF(e, 0) != -1 && F(EF(e, 0), c) == e) ind0 = c;
                if (EF(e, 1) != -1 && F(EF(e, 1), c) == e) ind0 = c;
            }
            const bool isBoundary = EF(e, 0) == -1 || EF(e, 1) == -1;
            if (isBoundary)
            {
                const int s = EF(e, 0) == -1 ? 1 : 0;
                tripsOneFormBack.emplace_back(3 * EF(e, s) + (ind0 == -1 ? ind1 : ind0), e, 1.);
            }
            else
            {
                tripsOneFormBack.emplace_back(3 * EF(e, 0) + ind0, e, 1.);
                tripsOneFormBack.emplace_back(3 * EF(e, 1) + ind1, e, 1.);
                tripsCBack.emplace_back(3 * EF(e, 0) + ind0, e, -1.);
                tripsCBack.emplace_back(3 * EF(e, 1) + ind1, e, 1.);
            }
        }
        OneFormToG3.setFromTriplets(tripsOneFormBack.begin(), tripsOneFormBack.end());
        CToG3.setFromTriplets(tripsCBack.begin(), tripsCBack.end());
        SMat U_inverse;
        get_U_inverse(F, EV, FE, EF,1, U_inverse);
        OneFormToG3 = U_inverse * OneFormToG3;
        CToG3 = U_inverse * CToG3;
        igl::cat(2, OneFormToG3, CToG3, W_inverse);
    }

    IGL_INLINE void get_W_inverse(const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& EV,
        const Eigen::MatrixXi& FE,
        const Eigen::MatrixXi& EF,
        const Eigen::VectorXi& matching,
        int N,
        Eigen::SparseMatrix<double>& W_inverse)
    {

        using SMat = Eigen::SparseMatrix<double>;
        
        const int faceCount = FE.rows();
        const int gammaCount = 3 * faceCount;
        const int edgeCount = EF.rows();
        SparseHelper curlToGamma3Construct(3 * faceCount * N, edgeCount * N, 2 * edgeCount * N);
        SparseHelper meanToGamma3Construct(3 * faceCount * N, edgeCount * N, 2 * edgeCount * N);
        for (int e = 0; e < edgeCount; e++)
        {
            int ind0 = -1, ind1 = -1;
            for (int c = 0; c < 3; ++c)
            {
                if (EF(e, 0) != -1 && F(EF(e, 0), c) == e) ind0 = c;
                if (EF(e, 1) != -1 && F(EF(e, 1), c) == e) ind0 = c;
            }
            const int lGamma = 3 * EF(e, 0) + ind0;
            const int rGamma = 3 * EF(e, 1) + ind1;
            for (int n = 0; n < N; n++)
            {
                const int level = modulo(n + matching(e), N);
                meanToGamma3Construct.addCoeff(lGamma + n * gammaCount, e + n * edgeCount, 1);
                meanToGamma3Construct.addCoeff(rGamma + level * gammaCount, e + n * edgeCount, 1);

                curlToGamma3Construct.addCoeff(lGamma + n * gammaCount, e + n * edgeCount, -1);
                curlToGamma3Construct.addCoeff(rGamma + level * gammaCount, e + n * edgeCount, 1);
            }
        }
        SMat UInv;
        get_U_inverse(F, EV, FE, EF, N, UInv);
        SMat OneFormToG3 = UInv * meanToGamma3Construct.toMat();
        SMat CToG3 = UInv * curlToGamma3Construct.toMat();

        igl::cat(2, OneFormToG3, CToG3, W_inverse);
    }

}
#endif