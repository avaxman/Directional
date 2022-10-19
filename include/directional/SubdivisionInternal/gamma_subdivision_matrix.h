#ifndef DIRECTIONAL_GAMMA_SUBDIVISION_MATRIX_H
#define DIRECTIONAL_GAMMA_SUBDIVISION_MATRIX_H
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include "build_shm_subdivision.h"
#include "shm_edge_topology.h"
#include <directional/Gamma_suite.h>
#include "directional/block_diag.h"

namespace directional
{
    /**
     * \brief Constructs the gamma subdivision matrix. 
     * 
     * Subdivides the packed gamma representation, where the edges corresponding to the field are given by the edges
     * opposite the first two corners of the face, i.e. opposite vertices F(f,0) and F(f,1) for face f.
     * \param[in] V Input coarse vertices
     * \param[in] F Input coarse face-to-vertex connectivity.
     * \param[in] level The target level of subdivision
     * \param[out] gammaSubdivisionMatrix The packed gamma subdivision operator, given as a sparse matrix.
     */
    inline void gamma_subdivision_matrix(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int level, Eigen::SparseMatrix<double>& gammaSubdivisionMatrix)
    {
        // Coarse level topology
        Eigen::MatrixXi E, EF, SFE, EI;
        // Fine level topology
        Eigen::MatrixXi FK, EK, EFK, SFEK, EIK;
        // Subdivision matrices
        Eigen::SparseMatrix<double> Se, Sc, G2ToOneform, G2ToG3, G3ToDecomp, DecompToG3K, G3ToG2K, S_Decomp;
        directional::shm_edge_topology(F, V.rows(), E, EF, EI, SFE);

        // Compute all subdivision operators. We don't really need Sv and Sf, but convenient for now.
        {
            Eigen::SparseMatrix<double> Sv, Sf;
            directional::subdivision::build_shm_subdivision(F, V, E, EF, EI, SFE, level, FK, EK, EFK, EIK, SFEK, Sv, Se, Sc, Sf);
        }
        // Compute the matrices for conversion to average-halfcurl representation.
        directional::Gamma3_To_Decomp(EF, EI, F.rows(), G3ToDecomp);
        directional::Gamma2_To_Gamma3(SFE, G2ToG3);
        // Back projection operators
        directional::Decomp_To_Gamma3(EFK, EIK, FK.rows(), DecompToG3K);
        directional::Gamma3_To_Gamma2(SFEK, G3ToG2K);
        // Compute subdivision operator on average-halfcurl representation
        directional::block_diag(std::vector<Eigen::SparseMatrix<double>*>{&Se, &Sc}, S_Decomp);

        // Return the full operator for subdividing packed gamma fields
        gammaSubdivisionMatrix = G3ToG2K * DecompToG3K * S_Decomp * G3ToDecomp * G2ToG3;
    }
    /**
     * \brief Constructs the gamma subdivision matrix.
     *
     * Subdivides the packed gamma representation, where the edges corresponding to the field are given by the edges
     * opposite the first two corners of the face, i.e. opposite vertices F(f,0) and F(f,1) for face f.
     * \param[in] V Input coarse vertices
     * \param[in] F Input coarse face-to-vertex connectivity.
     * \param[in] level The target level of subdivision
     * \param[out] gammaSubdivisionMatrix The packed gamma subdivision operator, given as a sparse matrix.
     */
    inline void gamma_subdivision_matrix(   const Eigen::MatrixXd& V,
                                            const Eigen::MatrixXi& F,
                                            const Eigen::MatrixXi& E,
                                            const Eigen::MatrixXi& EF,
                                            const Eigen::MatrixXi& SFE,
                                            const Eigen::MatrixXi& EI,
                                            int level,
                                            Eigen::SparseMatrix<double>& gammaSubdivisionMatrix)
    {
        // Fine level topology
        Eigen::MatrixXi FK, EK, EFK, SFEK, EIK;
        // Subdivision matrices
        Eigen::SparseMatrix<double> Se, Sc, G2ToOneform, G2ToG3, G3ToDecomp, DecompToG3K, G3ToG2K, S_Decomp;

        // Compute all subdivision operators. We don't really need Sv and Sf, but convenient for now.
        {
            Eigen::SparseMatrix<double> Sv, Sf;
            directional::subdivision::build_shm_subdivision(F, V, E, EF, EI, SFE, level, FK, EK, EFK, EIK, SFEK, Sv, Se, Sc, Sf);
        }
        // Compute the matrices for conversion to average-halfcurl representation.
        directional::Gamma3_To_Decomp(EF, EI, F.rows(), G3ToDecomp);
        directional::Gamma2_To_Gamma3(SFE, G2ToG3);
        // Back projection operators
        directional::Decomp_To_Gamma3(EFK, EIK, FK.rows(), DecompToG3K);
        directional::Gamma3_To_Gamma2(SFEK, G3ToG2K);
        // Compute subdivision operator on average-halfcurl representation
        directional::block_diag(std::vector<Eigen::SparseMatrix<double>*>{&Se, &Sc}, S_Decomp);

        // Return the full operator for subdividing packed gamma fields
        gammaSubdivisionMatrix = G3ToG2K * DecompToG3K * S_Decomp * G3ToDecomp * G2ToG3;
    }
    inline void gamma_subdivision_matrix(const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& EF,
        const Eigen::MatrixXi& SFE,
        const Eigen::MatrixXi& EI,
        int level,
        Eigen::SparseMatrix<double>& gammaSubdivisionMatrix,
        Eigen::SparseMatrix<double>& vertexSubdivisionMatrix,
        Eigen::MatrixXi&FK, 
        Eigen::MatrixXi& EK, 
        Eigen::MatrixXi& EFK, 
        Eigen::MatrixXi& SFEK, 
        Eigen::MatrixXi& EIK
    )
    {
        // Subdivision matrices
        Eigen::SparseMatrix<double> Se, Sc, G2ToOneform, G2ToG3, G3ToDecomp, DecompToG3K, G3ToG2K, S_Decomp;

        // Compute all subdivision operators. We don't really need Sv and Sf, but convenient for now.
        {
            Eigen::SparseMatrix<double> Sf;
            directional::subdivision::build_shm_subdivision(F, V, E, EF, EI, SFE, level, FK, EK, EFK, EIK, SFEK, vertexSubdivisionMatrix, Se, Sc, Sf);
        }
        // Compute the matrices for conversion to average-halfcurl representation.
        directional::Gamma3_To_Decomp(EF, EI, F.rows(), G3ToDecomp);
        directional::Gamma2_To_Gamma3(SFE, G2ToG3);
        // Back projection operators
        directional::Decomp_To_Gamma3(EFK, EIK, FK.rows(), DecompToG3K);
        directional::Gamma3_To_Gamma2(SFEK, G3ToG2K);
        // Compute subdivision operator on average-halfcurl representation
        directional::block_diag(std::vector<Eigen::SparseMatrix<double>*>{&Se, &Sc}, S_Decomp);

        // Return the full operator for subdividing packed gamma fields
        gammaSubdivisionMatrix = G3ToG2K * DecompToG3K * S_Decomp * G3ToDecomp * G2ToG3;
    }
}
#endif