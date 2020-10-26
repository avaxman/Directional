#ifndef DIRECTIONAL_GET_P_INVERSE_H
#define DIRECTIONAL_GET_P_INVERSE_H
#include <igl/igl_inline.h>
#include <Eigen/Sparse>
#include <igl/doublearea.h>
#include <igl/per_face_normals.h>

namespace directional
{
    /**
     * Construct the P^{-1} operator that converts gamma2 fields to PCVF fields.
     * Input:
     * V - |V| x 3 matrix of vertex locations
     * F - |F| x 3 matrix of vertex indices per face
     * EV - |E| x 2 matrix of edge to vertex mapping. Indices in the first column refer to start vertices, column 1 end vertices.
     * sFE - |F| x 6 matrix of face to edge mappings, including orientation. For face f, the edge opposite the ith corner is stored
     *			at sFE(f,i). Its orientation with respect to the face normal is located at sFE(f,i+3), where 0 denotes CCW and 1 CW.
     * EF - |E| x 2 matrix of edge to face mapping. Elements in column 0 refere to face to the left of the edge (in global orientation),
     *			for column 1 the faces to the right. For missing faces due to boundary, a -1 is used.
     * Output:
     * P - 3|F| x 2|F| sparse matrix mapping PCVFs to gamma3 elements. PCVF is assumed to be in [x1;y1;z1;x2;y2;z2;...etc] order in Matlab notation (column vec of xyz values)
     *
     */
    IGL_INLINE void get_P_inverse(const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& EV,
        const Eigen::MatrixXi& FE,
        Eigen::SparseMatrix<double>& P_inverse)
    {
        using namespace Eigen;
        P_inverse = Eigen::SparseMatrix<double>(3 * F.rows(), 2 * F.rows());
        std::vector<Triplet<double>> trips;
        trips.reserve(9 * F.rows());
        Eigen::VectorXd DA;
        Eigen::MatrixXd N;
        igl::doublearea(V, F, DA);
        igl::per_face_normals(V, F, N);
        for (int f = 0; f < F.rows(); f++)
        {
            RowVector3d localN = N.row(f);
            RowVector3d center = (V.row(F(f, 0)) + V.row(F(f, 1)) + V.row(F(f, 2))) / 3.0;
            const int e0ind = FE(f, 0), e1ind = FE(f, 1);
            // Compute gamma signs (could be easier when we can assume something about FE...)
            const double s0 = (V.row(EV(e0ind, 0))).cross(V.row(EV(e0ind, 1)) - center).dot(localN) > 0 ? 1.0 : -1.0;
            const double s1 = (V.row(EV(e1ind, 0))).cross(V.row(EV(e1ind, 1)) - center).dot(localN) > 0 ? 1.0 : -1.0;
            // Get the edges
            RowVector3d e0 = V.row(EV(e0ind, 1)) - V.row(EV(e0ind, 0));
            RowVector3d e1 = V.row(EV(e1ind, 1)) - V.row(EV(e1ind, 0));
            // Compute basis elements
            RowVector3d vals0 = localN.cross(e0) / DA(f);
            RowVector3d vals1 = -localN.cross(e1) / DA(f);
            const double signProd = s0 * s1;
            for (int i = 0; i < 3; i++)
            {
                trips.emplace_back(3 * f + i, 2 * f, signProd * vals1(i));
                trips.emplace_back(3 * f + i, 2 * f + 1, signProd * vals0(i));
            }
        }
        P_inverse.setFromTriplets(trips.begin(), trips.end());
    }
}
#endif