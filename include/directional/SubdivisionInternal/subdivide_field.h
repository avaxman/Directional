#ifndef DIRECTIONAL_SUBDIVIDE_FIELD_H
#define DIRECTIONAL_SUBDIVIDE_FIELD_H
#include <Eigen/Eigen>
#include "gamma_subdivision_matrix.h"
#include <directional/SubdivisionInternal/Gamma_suite.h>
namespace directional
{
    /**
     * \brief Computes a subdivided field from a coarse level field
     * 
     * 
     * \param V Coarse mesh vertices
     * \param F Coarse mesh face-to-vertex connectivity
     * \param pcvf The PCVF, given as a column vector of [xyzxyz...]' with a vector for every face in the mesh.
     * \param level The number of subdivisions to apply
     * \param[out] subdividedPcvf The PCVF in the fine level
     * \param[out] VK The vertices of the subdivision surface
     * \param[out] FK The face-to-vertex connectivity for the subdivision surface
     */
    inline void subdivide_field(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd& pcvf, int level, Eigen::MatrixXd& subdividedPcvf, Eigen::MatrixXd& VK, Eigen::MatrixXi& FK)
    {
        Eigen::SparseMatrix<double> gammaSubdivision, pcvfToG3, G3ToG2, Sv, G2KToPcvf;
        Eigen::VectorXd columnField(pcvf.rows() * pcvf.cols()), columnFieldK;
        for(int f = 0; f < pcvf.rows(); ++f)
        {
            for(int i = 0; i < 3; ++i)
            {
                columnField(3 * f + i) = pcvf(f, i);
            }
        }

        // Coarse level topology
        Eigen::MatrixXi E, EF, SFE, EI, EK, EFK, SFEK, EIK;
        directional::shm_edge_topology(F, V.rows(), E, EF, EI, SFE);
        directional::Gamma3_projector(V, F, E, SFE, EF, pcvfToG3);
        directional::Gamma3_To_Gamma2(SFE, G3ToG2);
        directional::gamma_subdivision_matrix(V, F, E, EF, SFE, EI, level, gammaSubdivision, Sv, FK, EK, EFK, SFEK, EIK);
        // Subdivide vertices
        VK = Sv * V;
        directional::Gamma2_reprojector(VK, FK, EK, SFEK, EFK, G2KToPcvf);
        columnFieldK  = G2KToPcvf * gammaSubdivision * G3ToG2 * pcvfToG3 * columnField;
        subdividedPcvf.resize(columnFieldK.size() / 3, 3);
        for (int f = 0; f < subdividedPcvf.rows(); ++f)
        {
            for (int i = 0; i < 3; ++i)
            {
                subdividedPcvf(f, i) = columnFieldK(3 * f + i);
            }
        }
    }
}
#endif 