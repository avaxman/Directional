#ifndef DIRECTIONAL_SUBDIVIDE_DIRECTIONALS_H
#define DIRECTIONAL_SUBDIVIDE_DIRECTIONALS_H
#include <Eigen/Eigen>
#include "build_directional_subdivision_operators.h"
#include "shm_edge_topology.h"
#include "shm_halfcurl_coefficients.h"
#include "shm_oneform_coefficients.h"
#include "Sc_directional_triplet_provider.h"
#include "Se_directional_triplet_provider.h"
#include "Sv_triplet_provider.h"
#include "build_subdivision_operators.h"
#include <directional/DirectionalGamma_Suite.h>

namespace directional
{
    namespace subdivision
    {

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
            for(int f = 0; f < rawField.size(); ++f)
            {
                rawField(f, n * 3) = columnDirectional(n*fCount + 3 * f);
                rawField(f, n * 3 + 1) = columnDirectional(n*fCount + 3 * f  + 1);
                rawField(f, n * 3 + 2) = columnDirectional(n*fCount + +3 * f + 2);
            }
        }
    }
    inline void rawfield_to_columndirectional(const Eigen::MatrixXd& rawField, int N, Eigen::VectorXd& columnDirectional)
    {
        columnDirectional.resize(rawField.rows() * rawField.cols());
        const int fCount = rawField.rows();
        for (int f = 0; f < rawField.size(); ++f)
        {
            for (int n = 0; n < N; ++n)
            {
                columnDirectional(n * 3 * fCount) = rawField(f, n * 3);
                columnDirectional(n * 3 * fCount + 1) = rawField(f, n * 3 + 1);
                columnDirectional(n * 3 * fCount + 2) = rawField(f, n * 3 + 2);
            }
        }
    }
     
    inline void subdivide_directionals(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& rawField, 
        const Eigen::VectorXi& matching, int targetLevel, Eigen::MatrixXd& Vk, Eigen::MatrixXi& Fk, Eigen::MatrixXd& rawFieldK, Eigen::VectorXi& matchingK,
        Eigen::MatrixXi& EFK, Eigen::MatrixXi& EVK)
    {

        const int N = rawField.cols() / 3;
        Eigen::MatrixXi E, EF, EI, SFE, EIK, SFEK;
        shm_edge_topology(F, V.rows(), E, EF, EI, SFE);
        std::vector<int> initialSizes = std::vector<int>({ (int)N * E.rows(), (int)N * E.rows() });
        std::vector<Eigen::SparseMatrix<double>> out, svOut;
        using coeffProv = coefficient_provider_t;

        auto Sv_provider = triplet_provider_wrapper<coeffProv>(loop_coefficients, Sv_triplet_provider<coeffProv>);
        auto Sc_directional_provider = directional_triplet_provider_wrapper<coeffProv>(shm_halfcurl_coefficients, Sc_directional_triplet_provider<coeffProv>);
        auto Se_directional_provider = directional_triplet_provider_wrapper<coeffProv>(shm_halfcurl_coefficients, Sc_directional_triplet_provider<coeffProv>);
        build_directional_subdivision_operators(V, F, E, EF, EI, SFE, matching, initialSizes, targetLevel, N,
            Fk, EVK, EFK, EIK, SFEK, matchingK, out, Sc_directional_provider, Se_directional_provider);
        // Construct regular vertex subdivision
        build_subdivision_operators(V, F, E, EF, EI, SFE, std::vector<int>({(int)V.rows()}), targetLevel, Fk, EVK, EFK, EIK, SFEK, svOut, Sv_provider);

        // Get fine level vertices
        Vk = svOut[0] * V;

        Eigen::SparseMatrix<double> G2_To_Decomp_0, Gamma2_To_PCVF_K, Matched_Gamma2_To_PCVF_K, S_Gamma_directional, S_Decomp, Decomp_To_G2K, columnDirectional_To_G2;
        // Construct fine gamma operator
        directional::Matched_Gamma2_To_AC(EI, EF, SFE, matching, N, G2_To_Decomp_0);
        directional::Matched_AC_To_Gamma2(EFK, SFEK, EIK, matchingK, N, Decomp_To_G2K);
        directional::Gamma2_reprojector(Vk, Fk, EVK, SFEK, EFK, Gamma2_To_PCVF_K);
        // Construct the full reprojection for all directionals. Since gammas are face local,
        // the matching is not needed
        igl::repmat(Gamma2_To_PCVF_K, N, 1, Matched_Gamma2_To_PCVF_K);
        directional::block_diag({ &out[1],&out[0] }, S_Decomp);
        S_Gamma_directional = Decomp_To_G2K * S_Decomp*G2_To_Decomp_0;

        Eigen::VectorXd columnDirectional, fineDirectional;
        directional::columndirectional_to_gamma2_matrix(V, F, E, SFE, EF, N, columnDirectional_To_G2);

        rawfield_to_columndirectional(rawField, N, columnDirectional);

        fineDirectional = Matched_Gamma2_To_PCVF_K * S_Gamma_directional * G2_To_Decomp_0 * columnDirectional_To_G2 * columnDirectional;

        columndirectional_to_rawfield(fineDirectional, N, rawFieldK);
    }
    }
}
#endif 