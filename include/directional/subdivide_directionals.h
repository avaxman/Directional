#ifndef DIRECTIONAL_SUBDIVIDE_DIRECTIONALS_H
#define DIRECTIONAL_SUBDIVIDE_DIRECTIONALS_H
#include <Eigen/Eigen>
#include <directional/SubdivisionInternal/build_directional_subdivision_operators.h>
#include <directional/SubdivisionInternal/shm_edge_topology.h>
#include <directional/SubdivisionInternal/shm_halfcurl_coefficients.h>
#include <directional/SubdivisionInternal/shm_oneform_coefficients.h>
#include <directional/SubdivisionInternal/Sc_directional_triplet_provider.h>
#include <directional/SubdivisionInternal/Se_directional_triplet_provider.h>
#include <directional/SubdivisionInternal/Gamma_suite.h>
#include <directional/SubdivisionInternal/Sv_triplet_provider.h>
#include <directional/SubdivisionInternal/build_subdivision_operators.h>
#include <directional/SubdivisionInternal/DirectionalGamma_Suite.h>

namespace directional
{
        inline void columnfunction_to_rawfunction(const Eigen::VectorXd& columnFunction, int N, int elementsPerN, Eigen::MatrixXd& rawFunction)
        {
            rawFunction.resize(columnFunction.size() / (N * elementsPerN), N * elementsPerN);
            const int rows = rawFunction.rows();
            for (int n = 0; n < N; ++n)
            {
                for (int i = 0; i < rows; ++i)
                {
                    for(int j = 0; j < elementsPerN; ++j)
                    {
                        rawFunction(i, n * elementsPerN + j) = columnFunction(n * elementsPerN * rows + elementsPerN * i + j);
                    }
                }
            }
        }
        /**
     */
        inline void rawfunction_to_columnfunction(const Eigen::MatrixXd& rawFunction, int N, int elementsPerN, Eigen::VectorXd& columnFunction)
        {
            columnFunction.resize(rawFunction.rows() * rawFunction.cols());
            const int rowCount = rawFunction.rows();
            for (int r = 0; r < rowCount; ++r)
            {
                for (int n = 0; n < N; ++n)
                {
                    for(int j = 0; j < elementsPerN; ++j)
                    {
                        columnFunction(n * elementsPerN * rowCount + elementsPerN * r) = rawFunction(r, n * elementsPerN);
                    }
                }
            }
        }
    /**
     * \brief Reshapes a column directional, given as [x1y1z1x2y2z2...]' per directional, stacked vertically, to 
     * a raw form with  [x1y1z1_0, x1y1z1_1,..] per row so that the row has all directionals per face
     * \param columnDirectional N * F * 3 x 1 matrix of directional field in column form
     * \param N Number of directional fields
     * \param rawField F x N * 3 marix containing raw field representation.
     */
    inline void columndirectional_to_rawfield(const Eigen::VectorXd& columnDirectional, int N, Eigen::MatrixXd& rawField)
    {
        rawField.resize(columnDirectional.size() / (N * 3), N * 3);
        const int fCount = rawField.rows();
        for(int n = 0; n < N; ++n)
        {
            for(int f = 0; f < fCount; ++f)
            {
                rawField(f, n * 3) = columnDirectional(n * 3 * fCount + 3 * f);
                rawField(f, n * 3 + 1) = columnDirectional(n * 3 * fCount + 3 * f + 1);
                rawField(f, n * 3 + 2) = columnDirectional(n * 3 * fCount + 3 * f + 2);
            }
        }
    }

    inline void rawfield_to_columndirectional(const Eigen::MatrixXd& rawField, int N, Eigen::VectorXd& columnDirectional)
    {
        columnDirectional.resize(rawField.rows() * rawField.cols());
        const int fCount = rawField.rows();
        for (int f = 0; f < fCount; ++f)
        {
            for (int n = 0; n < N; ++n)
            {
                columnDirectional(n * 3 * fCount + 3 * f) = rawField(f, n * 3);
                columnDirectional(n * 3 * fCount + 3 * f + 1) = rawField(f, n * 3 + 1);
                columnDirectional(n * 3 * fCount + 3 * f + 2) = rawField(f, n * 3 + 2);
            }
        }
    }
     
    inline void subdivide_directionals(const Eigen::MatrixXd& V, 
        const Eigen::MatrixXi& F, 
        const Eigen::MatrixXi& EV,
        const Eigen::MatrixXi& EF,
        const Eigen::MatrixXd& rawField,
        const Eigen::VectorXi& matching, 
        int targetLevel, 
        Eigen::MatrixXd& Vk,
        Eigen::MatrixXi& Fk,
        Eigen::MatrixXi& EVK,
        Eigen::MatrixXi& EFK,
        Eigen::MatrixXd& rawFieldK, 
        Eigen::VectorXi& matchingK)
    {
        Eigen::MatrixXi EI, SFE, EI_fine, SFE_fine;
        shm_edge_topology(F, EV, EF, EI,SFE);
        const int N = rawField.cols() / 3;
        std::vector<int> initialSizes = std::vector<int>({ (int)(N * EV.rows()), (int)(N * EV.rows()) });
        std::vector<Eigen::SparseMatrix<double>> out, svOut;
        using coeffProv = coefficient_provider_t;

        auto Sv_provider = triplet_provider_wrapper<coeffProv>(subdivision::loop_coefficients, subdivision::Sv_triplet_provider<coeffProv>);
        auto Sc_directional_provider = directional_triplet_provider_wrapper<coeffProv>(subdivision::shm_halfcurl_coefficients, subdivision::Sc_directional_triplet_provider<coeffProv>);
        auto Se_directional_provider = directional_triplet_provider_wrapper<coeffProv>(subdivision::shm_oneform_coefficients, subdivision::Se_directional_triplet_provider<coeffProv>);
        build_directional_subdivision_operators(V, F, EV, EF, EI, SFE, matching, initialSizes, targetLevel, N,
            Fk, EVK, EFK, EI_fine, SFE_fine, matchingK, out, Se_directional_provider, Sc_directional_provider);

        // Construct regular vertex subdivision
        build_subdivision_operators(V, F, EV, EF, EI, SFE, std::vector<int>({(int)V.rows()}), targetLevel, Fk, EVK, EFK, EI_fine, SFE_fine, svOut, Sv_provider);

        // Get fine level vertices
        Vk = svOut[0] * V;

        Eigen::SparseMatrix<double> G2_To_Decomp_0, Gamma2_To_PCVF_K, Matched_Gamma2_To_PCVF_K, S_Gamma_directional, S_Decomp, Decomp_To_G2K, columnDirectional_To_G2;
        // Construct fine gamma operator
        directional::Matched_Gamma2_To_AC(EI, EF, SFE, matching, N, G2_To_Decomp_0);
        directional::Matched_AC_To_Gamma2(EFK, SFE_fine, EI_fine, matchingK, N, Decomp_To_G2K);
        directional::Gamma2_reprojector(Vk, Fk, EVK, SFE_fine, EFK, Gamma2_To_PCVF_K);
        // Construct the full reprojection for all directionals. Since gammas are face local,
        // the matching is not needed
        {
            std::vector<Eigen::SparseMatrix<double>*> base(N, &Gamma2_To_PCVF_K);
            directional::block_diag(base, Matched_Gamma2_To_PCVF_K);
        }

        directional::block_diag({ &out[0],&out[1] }, S_Decomp);

        Eigen::VectorXd columnDirectional, fineDirectional;

        // Convert rawfield to column directional for applying subdivision
        rawfield_to_columndirectional(rawField, N, columnDirectional);
        // Get matrix to convert column directional to gamma2 elements
        directional::columndirectional_to_gamma2_matrix(V, F, EV, SFE, EF, N, columnDirectional_To_G2);

        // The directional gamma subdivision operator
        S_Gamma_directional = Decomp_To_G2K * S_Decomp*G2_To_Decomp_0;

        Eigen::VectorXd subdivAc = S_Decomp * G2_To_Decomp_0 * columnDirectional_To_G2 * columnDirectional;
        std::cout << "Max C after subdivision:" << subdivAc.block(subdivAc.rows() / 2, 0, subdivAc.rows() / 2, 0).maxCoeff() << std::endl;
        // Compute the fine directional
        fineDirectional = Matched_Gamma2_To_PCVF_K * S_Gamma_directional * columnDirectional_To_G2 * columnDirectional;

        // Convert resulting column directional in fine level back to rawfield format
        columndirectional_to_rawfield(fineDirectional, N, rawFieldK);
    }
}
#endif 