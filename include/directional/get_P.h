#ifndef DIRECTIONAL_GET_P_H
#define DIRECTIONAL_GET_P_H
#include <igl/igl_inline.h>
#include <Eigen/Sparse>
namespace directional
{
    /**
      * Projects a PCVF into the gamma2 space by taking the innerproduct of the globally directed edges per face with the
      * vector field value per face, omitting the third edge per face. Projections are stored at the index of the corner opposite
      * the edge within the face.
      * Input:
      * V - |V| x 3 matrix of vertex locations
      * F - |F| x 3 matrix of vertex indices per face
      * EV - |E| x 2 matrix of edge to vertex mapping. Indices in the first column refer to start vertices, column 1 end vertices.
      * sFE - |F| x 6 matrix of face to edge mappings, including orientation. For face f, the edge opposite the ith corner is stored
      *			at sFE(f,i). Its orientation with respect to the face normal is located at sFE(f,i+3), where 0 denotes CCW and 1 CW.
      * EF - |E| x 2 matrix of edge to face mapping. Elements in column 0 refere to face to the left of the edge (in global orientation),
      *			for column 1 the faces to the right. For missing faces due to boundary, a -1 is used.
      * Output:
      * P - 2|F| x 3|F| sparse matrix mapping PCVFs to gamma3 elements. PCVF is assumed to be in [x1;y1;z1;x2;y2;z2;...etc] order in Matlab notation (column vec of xyz values)
      * The gamma's are at the edges specified by FE(f,0) and FE(f,1) for face f.
      *
      */
    IGL_INLINE void get_P(const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& EV,
        const Eigen::MatrixXi& FE,
        Eigen::SparseMatrix<double>& P)
    {
        // Projects FEM vector field into Gamma3 space.
        P = Eigen::SparseMatrix<double>(2 * F.rows(), 3 * F.rows());

        std::vector<Eigen::Triplet<double>> trips;
        trips.reserve(9 * F.rows());
        for (int f = 0; f < F.rows(); f++) {
            // Iterate over corners
            for (int c = 0; c < 2; c++) {
                //Get edge opposite corner in global orientation
                const int eI = FE(f, c);
                Eigen::RowVector3d edge = V.row(EV(eI, 1)) - V.row(EV(eI, 0));
                for (int i = 0; i < 3; i++) {
                    trips.emplace_back(2 * f + c, 3 * f + i, edge[i]);
                }
            }
        }
        P.setFromTriplets(trips.begin(), trips.end());
    }
}
#endif