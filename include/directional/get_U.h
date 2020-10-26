#ifndef DIRECTIONAL_GET_U_H
#define DIRECTIONAL_GET_U_H
#include <igl/igl_inline.h>
#include <Eigen/Sparse>

namespace directional
{
    IGL_INLINE bool isCCWOriented(const Eigen::RowVector3i& faceVerts, const Eigen::RowVector2i& edge)
    {
        for(int i = 0; i < 3; ++i)
        {
            if(faceVerts(i) == edge(0))
            {
                if (faceVerts((i + 1) % 3) == edge(1)) return true;
                return false;
            }
        }
        throw std::runtime_error("Edge and face do not share any vertex");
    }

    IGL_INLINE void get_U(const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& EV,
        const Eigen::MatrixXi& FE,
        const Eigen::MatrixXi& EF,
        int N,
        Eigen::SparseMatrix<double>& U)
    {
        if(N == 1)
        {
            U = Eigen::SparseMatrix<double>(3 * FE.rows(), 2 * FE.rows());
            using namespace Eigen;
            using namespace std;
            vector<Triplet<double>> trips;
            trips.reserve(4 * FE.rows());
            for (int f = 0; f < FE.rows(); f++)
            {
                RowVector3d signs;
                for (int c = 0; c < 3; ++c) signs(c) = isCCWOriented(F.row(f), EV.row(FE(f, c))) ? 1.0 : -1.0;

                trips.emplace_back(3 * f, 2 * f, 1.);
                trips.emplace_back(3 * f + 1, 2 * f + 1, 1.);
                // Reconstruct last halfedge form from zero sum constraint.
                trips.emplace_back(3 * f + 2, 2 * f, -signs(0)*signs(2));
                trips.emplace_back(3 * f + 2, 2 * f + 1, -signs(1)*signs(2));
            }
            U.setFromTriplets(trips.begin(), trips.end());
        }
        else
        {
            Eigen::SparseMatrix<double> U1;
            get_U(F, EV, FE, EF, 1, U1);
            std::vector<Eigen::SparseMatrix<double>*> copies(N, &U1);
            directional::block_diag(copies, U);
        }
    }
}
#endif