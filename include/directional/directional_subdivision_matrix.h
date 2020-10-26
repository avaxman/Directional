#ifndef DIRECTIONAL_DIRECTIONAL_SUBDIVISION_MATRIX_H
#define DIRECTIONAL_DIRECTIONAL_SUBDIVISION_MATRIX_H
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <igl/edge_topology.h>

namespace directional{
    void subdivide_mesh_topology(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::MatrixXi& F_subdivided, int& newVertexCount)
    {
        // Reserve space
        F_subdivided.setConstant(F.rows() * 4, F.cols(), -1);

        Eigen::MatrixXi EV, FE, EF;

        // Compute the edge topology
        igl::edge_topology(V,F, EV, FE, EF);

        
        for(int f = 0; f < F.rows(); ++f){
            // Let F(f,c) be the vertex at corner c of face f. For the subdivided topology
            // we will demand that face 4 * f + c is neighbouring the vertex F(f,c) in the original, 
            // with that vertex in the same corner as the original.
            for(int c = 0; c < 3; ++ c){
                F_subdivided(4*f+c, c) = F(f,c);
            }
        }
        // Determine via edges to which faces the odd vertex elements connect
        for(int e = 0; e < EF.rows(); ++e)
        {
            // Left and right face
            for(int f = 0; f < 2; ++f){
                // Find corners to which the edge connects on this face
                int cs[2] = {-1,-1};
                // The vertex in the face that is not part of the edge
                int nonConnectedC = -1;
                for(int c = 0; c < 3; ++c){
                    if(F(f,c) != EV(e,0) && F(f,c)!=EV(e,1)){
                       nonConnectedC = c;
                       cs[0] = (c+1)%3; 
                       cs[1] = (c+2)%3;
                       break;
                    }
                }
                // Assign the odd element vertex to the correct corner of the face for the even face elements
                F_subdivided(4*f+cs[0], (cs[0] + 1)%3 ) = V.rows() + e;
                F_subdivided(4*f+cs[1], (cs[1] + 2)%3 ) = V.rows() + e;
                // Assign the edge to the odd face element, storing the odd vertex at the corner with 
                // the same corner index as the corner opposite it in the original face
                F_subdivided(4*f+3, nonConnectedC) = V.rows() + e;
            }
        }
        newVertexCount = V.rows() + EV.rows();
    }

    void directional_subdivision_matrix(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXi& matching, 
        int directionalsPerFace, Eigen::SparseMatrix<double>& output)
    {
        const int vertexCount = V.rows();
        Eigen::SparseMatrix<double> S1, S_epsstar,W;

        directional::get_S_1(vertexCount, F, matching, directionalsPerFace, S1);
        directional::get_S_epsstar(vertexCount, F, matching, directionalsPerFace, S_epsstar);
        directional::get_W(vertexCount, F, matching, directionalsPerFace, W);
        Eigen::MatrixXi F_subdivded;
        int newVertexCount=-1;
        subdivide_mesh_topology(V, F, F_subdivded, newVertexCount);
        Eigen::VectorXi matching_subdivided;
        directional::subdivide_matching(V,F, matching, matching_subdivided);
        directional::get_W_inv(newVertexCount, F_subdivded, matching_subdivided, directionalsPerFace, W);
        Eigen::SparseMatrix<double> blockMat;
        directional::block_diag({&S1, &S_epsstar}, blockMat);
        output = P_inv * W_inv * blockMat * W * P;
    }
}

#endif