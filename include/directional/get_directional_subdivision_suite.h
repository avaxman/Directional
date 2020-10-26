#ifndef DIRECIONAL_GET_DIRECTIONAL_SUBDIVISION_SUITE_H
#define DIRECIONAL_GET_DIRECTIONAL_SUBDIVISION_SUITE_H
#include <Eigen/Core>
#include <Eigen/Sparse>
#include "SubdivisionInternal/shm_edge_topology.h"
#include "SubdivisionInternal/Sv_triplet_provider.h"
#include "SubdivisionInternal/Sc_triplet_provider.h"
#include "SubdivisionInternal/Se_triplet_provider.h"
#include "SubdivisionInternal/Sc_directional_triplet_provider.h"
#include "SubdivisionInternal/Se_directional_triplet_provider.h"
#include "SubdivisionInternal/Sf_triplet_provider.h"
#include "SubdivisionInternal/build_directional_subdivision_operators.h"
#include "SubdivisionInternal/shm_halfcurl_coefficients.h"
#include "SubdivisionInternal/shm_oneform_coefficients.h"
#include "SubdivisionInternal/hbspline_coefficients.h"
#include "get_W.h"
#include "get_P.h"
#include "get_P_inverse.h"
#include "get_W_inverse.h"

namespace directional{
    

    /**
     * Input:
     * - Vcoarse |V| x 3 matrix of vertex coordinates (row per coordinate)
     * - Fcoarse |F| x 3, connectivity between faces and vertices, where each row gives the 3 vertices of the face, in CCW order relative to the face normal.
     * - EVcoarse |E| x 2, connectivity between edge and vertex, where edge e is directed from EVcoarse(e,0) to EVcoarse(e,1)
     */
    void get_pcvf_subdivision_suite(
        const Eigen::MatrixXd& VCoarse,
        const Eigen::MatrixXi& FCoarse,
        const Eigen::MatrixXi& EVCoarse,
        int subdivisionLevel,
        Eigen::SparseMatrix<double>& S_epsstar,
        Eigen::SparseMatrix<double>& S_0,
        Eigen::SparseMatrix<double>& S_1,
        Eigen::SparseMatrix<double>& S_2,
        Eigen::SparseMatrix<double>& WCoarse,
        Eigen::SparseMatrix<double>& PCoarse,
        Eigen::MatrixXi& EVFine,
        Eigen::MatrixXi& FFine,
        Eigen::SparseMatrix<double>& WInvFine,
        Eigen::SparseMatrix<double>& PInvFine
    ) {
        Eigen::MatrixXi EICoarse, SFECoarse, EFCoarse, EIFine, SFEFine, EFFine;
        shm_edge_topology(FCoarse, EVCoarse, std::ref(EFCoarse), EICoarse, SFECoarse);
        // Initial sizes of the subdivision operators
        std::vector<int> initialSizes = std::vector<int>({ (int)VCoarse.rows(), (int)(EVCoarse.rows()), (int)(EVCoarse.rows()), (int)FCoarse.rows() });
        std::vector<Eigen::SparseMatrix<double>> subdivisionOperators;

        using coeffProv = coefficient_provider_t;

        // Providers for all subdivision operators
        auto Sv_provider = triplet_provider_wrapper<coeffProv>(subdivision::loop_coefficients, subdivision::Sv_triplet_provider<coeffProv>);
        auto Sc_provider = triplet_provider_wrapper<coeffProv>(subdivision::shm_halfcurl_coefficients, subdivision::Sc_triplet_provider<coeffProv>);
        auto Se_provider = triplet_provider_wrapper<coeffProv>(subdivision::shm_oneform_coefficients, subdivision::Se_triplet_provider<coeffProv>);
        auto Sf_provider = triplet_provider_wrapper<coeffProv>(subdivision::hbspline_coefficients, subdivision::Sf_triplet_provider<coeffProv>);

        // Construct regular vertex subdivision
        build_subdivision_operators(VCoarse, FCoarse, EVCoarse, EFCoarse, EICoarse, SFECoarse, std::vector<int>({ (int)VCoarse.rows() }), subdivisionLevel,
            FFine, EVFine, EFFine, EIFine, SFEFine, subdivisionOperators, Sv_provider, Se_provider, Sc_provider, Sf_provider);
        S_0 = subdivisionOperators[0];
        S_1 = subdivisionOperators[1];
        S_epsstar = subdivisionOperators[2];
        S_2 = subdivisionOperators[3];

        Eigen::VectorXd VCoarseVec;
        VCoarseVec.resize(3 * VCoarse.rows());
        for (int v = 0; v < VCoarse.rows(); ++v) VCoarseVec.middleRows(3 * v, 3) = VCoarse.row(v).transpose();
        Eigen::VectorXd VFineVec = S_0 * VCoarseVec;
        Eigen::MatrixXd VFine;
        VFine.setConstant(VCoarseVec.rows() / 3, 3, -1);
        for (int v = 0; v < VFine.rows(); ++v)
        {
            VFine.row(v) = VFineVec.middleRows(3 * v, 3).transpose();
        }

        get_W(FCoarse, EVCoarse, SFECoarse.leftCols(3), EFCoarse, WCoarse);
        get_P(VCoarse, FCoarse, EVCoarse, SFECoarse.leftCols(3), PCoarse);
        get_P_inverse(VFine, FFine, EVFine, SFEFine.leftCols(3), PInvFine);
        get_W_inverse(VFine, FFine, EVFine, SFEFine.leftCols(3), WInvFine);
    }


    /**
     * Input:
     * - Vcoarse |V| x 3 matrix of vertex coordinates (row per coordinate)
     * - Fcoarse |F| x 3, connectivity between faces and vertices, where each row gives the 3 vertices of the face, in CCW order relative to the face normal.
     * - EVcoarse |E| x 2, connectivity between edge and vertex, where edge e is directed from EVcoarse(e,0) to EVcoarse(e,1)
     * - matching |E| x 1, the matching between the left and right face of edge e, where directional at k in the left face maps to (k+ matching(e)) % N in the right face
     * - N, the number of directionals per face.
     */
    void get_directional_subdivision_suite(
        const Eigen::MatrixXd& VCoarse,
        const Eigen::MatrixXi& FCoarse,
        const Eigen::MatrixXi& EVCoarse,
        const Eigen::VectorXi& matchingCoarse,
        int N,
        int subdivisionLevel,
        Eigen::SparseMatrix<double>& S_epsstar,
        Eigen::SparseMatrix<double>& S_0,
        Eigen::SparseMatrix<double>& S_1,
        Eigen::SparseMatrix<double>& S_2,
        Eigen::SparseMatrix<double>& WCoarse,
        Eigen::SparseMatrix<double>& PCoarse,
        Eigen::VectorXi& matchingFine,
        Eigen::MatrixXi& EVFine,
        Eigen::MatrixXi& FFine,
        Eigen::SparseMatrix<double>& WInvFine,
        Eigen::SparseMatrix<double>& PInvFine
    ){
        if(N==1)
        {
            get_pcvf_subdivision_suite(VCoarse, FCoarse, EVCoarse, subdivisionLevel, S_epsstar, S_0, S_1, S_2, WCoarse, PCoarse, EVFine, FFine, WInvFine, PInvFine);
            return;
        }

        Eigen::MatrixXi EICoarse, SFECoarse, EFCoarse, EIFine, SFEFine, EFFine;
        shm_edge_topology(FCoarse, EVCoarse, std::ref(EFCoarse), EICoarse, SFECoarse);
        std::vector<int> initialSizes = std::vector<int>({ (int)(N * EVCoarse.rows()), (int)(N * EVCoarse.rows()) });
        std::vector<Eigen::SparseMatrix<double>> out, svOut;
        using coeffProv = coefficient_provider_t;

        auto Sv_provider = triplet_provider_wrapper<coeffProv>(subdivision::loop_coefficients, subdivision::Sv_triplet_provider<coeffProv>);
        auto Sc_directional_provider = directional_triplet_provider_wrapper<coeffProv>(subdivision::shm_halfcurl_coefficients, subdivision::Sc_directional_triplet_provider<coeffProv>);
        auto Se_directional_provider = directional_triplet_provider_wrapper<coeffProv>(subdivision::shm_oneform_coefficients, subdivision::Se_directional_triplet_provider<coeffProv>);
        //auto Sf_directional_provider = directional_triplet_provider_wrapper<coeffProv>(subdivision::hbspline_coefficients, subdivision::Sf_directional_triplet_provider<coeffProv>);
        build_directional_subdivision_operators(VCoarse, FCoarse, EVCoarse, EFCoarse, EICoarse, SFECoarse, matchingCoarse, initialSizes, subdivisionLevel, N,
            FFine, EVFine, EFFine, EIFine, SFEFine, matchingFine, out, Se_directional_provider, Sc_directional_provider);// , Sf_directional_provider);
        S_1 = out[0];
        S_epsstar = out[1];

        // Construct regular vertex subdivision
        build_subdivision_operators(VCoarse, FCoarse, EVCoarse, EFCoarse, EICoarse, SFECoarse, std::vector<int>({ (int)VCoarse.rows() }), subdivisionLevel,
            FFine, EVFine, EFFine, EIFine, SFEFine, svOut, Sv_provider);
        S_0 = svOut[0];

        Eigen::VectorXd VCoarseVec;
        VCoarseVec.resize(3 * VCoarse.rows());
        for (int v = 0; v < VCoarse.rows(); ++v) VCoarseVec.middleRows(3 * v, 3) = VCoarse.row(v).transpose();
        Eigen::VectorXd VFineVec = S_0 * VCoarseVec;
        Eigen::MatrixXd VFine;
        VFine.setConstant(VCoarseVec.rows() / 3, 3, -1);
        for (int v = 0; v < VFine.rows(); ++v)
        {
            VFine.row(v) = VFineVec.middleRows(3 * v, 3).transpose();
        }

        get_W(FCoarse, EVCoarse, SFECoarse.leftCols(3), EFCoarse, WCoarse);
        get_P(VCoarse, FCoarse, EVCoarse, SFECoarse.leftCols(3), PCoarse);
        get_P_inverse(VFine, FFine, EVFine, SFEFine.leftCols(3), PInvFine);
        get_W_inverse(VFine, FFine, EVFine, SFEFine.leftCols(3), WInvFine);
    }

}
#endif