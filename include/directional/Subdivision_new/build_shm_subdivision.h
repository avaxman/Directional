#ifndef DIRECTIONAL_BUILD_SHM_SUBDIVISION_H
#define DIRECTIONAL_BUILD_SHM_SUBDIVISION_H
#include "build_subdivision_operators.h"
#include "Sv_triplet_provider.h"
#include "Se_triplet_provider.h"
#include "Sc_triplet_provider.h"
#include "Sf_triplet_provider.h"
// Coefficient provideres
#include "loop_coefficients.h"
#include "shm_halfcurl_coefficients.h"
#include "shm_oneform_coefficients.h"
#include "hbspline_coefficients.h"
#include <vector>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
namespace directional{
    namespace subdivision{
        void build_shm_subdivision(
            const Eigen::MatrixXi& F0,
            const Eigen::MatrixXd& V0,
            const Eigen::MatrixXi& E0,
            const Eigen::MatrixXi& EF0,
            const Eigen::MatrixXi& EI0,
            const Eigen::MatrixXi& SFE0,
            int level,
            Eigen::MatrixXi& FK,
            Eigen::MatrixXi& EK,
            Eigen::MatrixXi& EFK,
            Eigen::MatrixXi& EIK,
            Eigen::MatrixXi& SFEK,
            Eigen::SparseMatrix<double>& Sv,
            Eigen::SparseMatrix<double>& Se,
            Eigen::SparseMatrix<double>& Sc,
            Eigen::SparseMatrix<double>& Sf
        )
        {
            std::vector<int> initialSizes = {(int)V0.rows(), (int)E0.rows(), (int)E0.rows(), (int)SFE0.rows()};
            std::vector<Eigen::SparseMatrix<double>> output;
            using coeffProv =coefficient_provider_t;
            auto Sv_provider = triplet_provider_wrapper<coeffProv>(loop_coefficients, Sv_triplet_provider<coeffProv>);
            auto Se_provider = triplet_provider_wrapper<coeffProv>(shm_oneform_coefficients, Se_triplet_provider<coeffProv>);
            auto Sc_provider = triplet_provider_wrapper<coeffProv>(shm_halfcurl_coefficients, Sc_triplet_provider<coeffProv>);
            auto Sf_provider = triplet_provider_wrapper<coeffProv>(hbspline_coefficients, Sf_triplet_provider<coeffProv>);
            build_subdivision_operators(F0,V0,E0,EF0,EI0,SFE0,initialSizes,level, FK, EK, EFK, EIK, SFEK, output, Sv_provider, Se_provider, Sc_provider, Sf_provider);
            Sv = output[0];
            Se = output[1];
            Sc = output[2];
            Sf = output[3];
        }
    }
}
#endif