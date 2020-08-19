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
#include <directional/rawfield_to_columndirectional.h>
#include <directional/columndirectional_to_rawfield.h>

namespace directional
{
    /**
     * Subdivides a raw field directional on a coarse mesh defined by V,F to a raw field directional in subdivision level 'targetLevel', 
     * on the mesh as given by output V_fine, F_fine. Assumes a matching is given, this matching will be fixed during subdivision.
     * Input:
     * - V |V| x 3 matrix of vertex coordinates
     * - F |F| x  3 matrix of face to vertex connectivity, given in CCW order relative to the normal
     * - EV |E| x  2 matrix of edge to vertex connectivity, such that edge e is between vertices EV(e,0) and EV(e,1).
     * - EF |E| x  2 matrix of edge to face connectivity, such that face EF(e,0) is to the left of e and EF(e,1) is to the right of e.
     * - rawField |F| x (3 * N) matrix containing the N-directional raw field representation
     * - matching |E| x 1 vector describing the directional matching over edge e such that directional k in face EF(e,0) matches to 
     * directional (matching(e)+k)% N in face EF(e,1).
     * - targetLevel The target subdivision level to subdivide to
     * Output:
     * - V_fine |V_fine| x 3 matrix of fine mesh vertex coordinates
     * - F_fine |F_fine| x 3 matrix of face to vertex connectivity of fine mesh
     * - EV_fine |E| x  2 matrix of edge to vertex connectivity for fine mesh.
     * - EF_fine |E| x  2 matrix of edge to face connectivity fine mesh.
     * - rawField_fine |F_fine| x (3 * N) matrix containing the fine level N-directional raw field
     * - matching_fine |E_fine| x 1 vector containing the fine matching.
     */
    inline void subdivide_directionals(const Eigen::MatrixXd& V, 
        const Eigen::MatrixXi& F, 
        const Eigen::MatrixXi& EV,
        const Eigen::MatrixXi& EF,
        const Eigen::MatrixXd& rawField,
        const Eigen::VectorXi& matching, 
        int targetLevel, 
        Eigen::MatrixXd& V_fine,
        Eigen::MatrixXi& F_fine,
        Eigen::MatrixXi& EV_fine,
        Eigen::MatrixXi& EF_fine,
        Eigen::MatrixXd& rawField_fine, 
        Eigen::VectorXi& matching_fine)
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
            F_fine, EV_fine, EF_fine, EI_fine, SFE_fine, matching_fine, out, Se_directional_provider, Sc_directional_provider);

        // Construct regular vertex subdivision
        build_subdivision_operators(V, F, EV, EF, EI, SFE, std::vector<int>({(int)V.rows()}), targetLevel, 
        F_fine, EV_fine, EF_fine, EI_fine, SFE_fine, svOut, Sv_provider);

        // Get fine level vertices
        V_fine = svOut[0] * V;

        Eigen::SparseMatrix<double> G2_To_Decomp_0, Gamma2_To_PCVF_K, Matched_Gamma2_To_PCVF_K, S_Gamma_directional, S_Decomp, Decomp_To_G2K, columnDirectional_To_G2;
        // Construct fine gamma operator
        directional::Matched_Gamma2_To_AC(EI, EF, SFE, matching, N, G2_To_Decomp_0);
        directional::Matched_AC_To_Gamma2(EF_fine, SFE_fine, EI_fine, matching_fine, N, Decomp_To_G2K);
        directional::Gamma2_reprojector(V_fine, F_fine, EV_fine, SFE_fine, EF_fine, Gamma2_To_PCVF_K);
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

        // Compute the fine directional
        fineDirectional = Matched_Gamma2_To_PCVF_K * S_Gamma_directional * columnDirectional_To_G2 * columnDirectional;

        // Convert resulting column directional in fine level back to rawfield format
        columndirectional_to_rawfield(fineDirectional, N, rawField_fine);
    }

    /**
     * Subdivides a raw field directional on a coarse mesh defined by V,F to a raw field directional in subdivision level 'targetLevel', 
     * on the mesh as given by output V_fine, F_fine. Determines the matching for the directional field by applying directional::curl_matching
     * on the coarse mesh and field, and fixing this during subdivision.
     * Input:
     * - V |V| x 3 matrix of vertex coordinates
     * - F |F| x  3 matrix of face to vertex connectivity, given in CCW order relative to the normal
     * - rawField |F| x (3 * N) matrix containing the N-directional raw field representation
     * - targetLevel The target subdivision level to subdivide to
     * Output:
     * - V_fine |V_fine| x 3 matrix of fine mesh vertex coordinates
     * - F_fine |F_fine| x 3 matrix of face to vertex connectivity of fine mesh
     * - rawField_fine |F_fine| x (3 * N) matrix containing the fine level N-directional raw field
     */
    inline void subdivide_directionals(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXd& rawField,
        int targetLevel,
        Eigen::MatrixXd& V_fine,
        Eigen::MatrixXi& F_fine,
        Eigen::MatrixXd& rawField_fine)
    {
        // Compute internal edge topology
        Eigen::MatrixXi EV, EF, EI, SFE, EI_et, FE_et, EI_fine, SFE_fine;
        shm_edge_topology(F, V.rows(), EV, EF, EI, SFE);
        shm_edge_topology_to_igledgetopology(F, EV, EF, SFE, EI_et, FE_et);
        // Compute curl matching
        Eigen::VectorXi matching, matching_fine;
        {
            Eigen::VectorXd effort, curlNorm;
            directional::curl_matching(V, F, EV, EF, FE_et, rawField, matching, effort, curlNorm);
        }

        subdivide_directionals(V, F, EV, EF, rawField, matching, targetLevel, V_fine, F_fine, EV_fine, EF_fine, rawField_fine, matching_fine);
    }
}
#endif 