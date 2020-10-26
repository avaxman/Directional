#ifndef DIRECTIONAL_GET_U_INVERSE_H
#define DIRECTIONAL_GET_U_INVERSE_H
#include <igl/igl_inline.h>
#include <Eigen/Sparse>

namespace directional
{

    IGL_INLINE void get_U_inverse(const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& EV,
        const Eigen::MatrixXi& FE,
        const Eigen::MatrixXi& EF,
        int N,
        Eigen::SparseMatrix<double>& U_inverse)
    {
        if(N==1)
        {
            U_inverse = Eigen::SparseMatrix<double>(2 * FE.rows(), 3 * FE.rows());
            using namespace Eigen;
            using namespace std;
            vector<Triplet<double>> trips;
            trips.reserve(2 * FE.rows());
            for (int f = 0; f < FE.rows(); f++)
            {
                trips.emplace_back(2 * f, 3 * f, 1.);
                trips.emplace_back(2 * f + 1, 3 * f + 1, 1.);
                // Drop third halfedge form element per face
            }
            U_inverse.setFromTriplets(trips.begin(), trips.end());
        }
        else
        {
            Eigen::SparseMatrix<double> UInv1;
            get_U_inverse(F, EV, FE, EF,1, UInv1);
            std::vector<Eigen::SparseMatrix<double>*> copies(N, &UInv1);
            directional::block_diag(copies, U_inverse);
        }
    }
}
#endif