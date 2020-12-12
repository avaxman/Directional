#include "catch2/catch.hpp"
#include "directional/get_directional_subdivision_suite.h"
#include <igl/read_triangle_mesh.h>
#include "test_includes.h"
#include <igl/adjacency_matrix.h>

// Subdivision level to check
const int baseSubdivisionLevel = 1;

void adjacencyMatrix(const Eigen::MatrixXi& adjacency, Eigen::SparseMatrix<double>& adjacencyMatrix)
{
    adjacencyMatrix = Eigen::SparseMatrix<double>(adjacency.rows(), adjacency.maxCoeff() + 1);
    std::vector<Eigen::Triplet<double>> trips;
    trips.reserve(adjacency.rows() * adjacency.cols());
    for(auto i = 0; i < adjacency.rows(); ++i)
    {
        for(auto j = 0; j < adjacency.cols(); ++j)
        {
            if (adjacency(i, j) < 0) continue;
            trips.emplace_back(i, adjacency(i, j), 1.0);
        }
    }
    adjacencyMatrix.setFromTriplets(trips.begin(), trips.end());
}
void dec_D1(const TriangleMesh& data, Eigen::SparseMatrix<double>& d1)
{
    d1 = Eigen::SparseMatrix<double>(data.F.rows(), data.E.rows());
    std::vector<Eigen::Triplet<double>> trips;
    for(auto f = 0; f < data.FE.rows(); ++f)
    {
        for(auto j = 0; j < data.FE.cols(); ++j)
        {
            const auto e = data.FE(f, j);
            double v = 0;
            for(int c = 0; c < 3; ++c)
            {
                if(data.F(f,c) == data.E(e,0))
                {
                    if(data.F(f,(c+1)%3) == data.E(e,1))
                    {
                        v = 1.0;
                    }
                    else
                    {
                        v = -1.0;
                    }
                }
            }
            trips.emplace_back(f, e, v);
        }
    }
    d1.setFromTriplets(trips.begin(), trips.end());
}

std::vector<std::string> testCaseFiles = {"bimba.off","cathead.off","chipped-torus.obj","half-torus.obj","horser.off","tester-sphere.off","torus.obj","bunny2.off" };

TEST_CASE("SHM commutation","[shm_commutation]")
{
    for(int i = 0; i < testCaseFiles.size(); ++i)
    {
        TriangleMesh coarseMesh, fineMesh;
        // Subdivision operators
        Eigen::SparseMatrix<double> S_0, S_1, S_2, S_epsstar, S_gamma;
        // Operators that convert between PCVF and Gamma, and Gamma and mean-curl decomposition.
        Eigen::SparseMatrix<double> WCoarse, PCoarse, PInvCoarse, WInvFine, PInvFine, PFine;

        FEM_operators coarseOperators, fineOperators;

        // Load mesh via IGL
        coarseMesh.read(testCaseFiles[i]);
        coarseMesh.compute_edge_topology();

        directional::get_pcvf_subdivision_suite(coarseMesh.V, coarseMesh.F, coarseMesh.E, baseSubdivisionLevel, S_epsstar, S_0, S_1, S_2, S_gamma, WCoarse, PCoarse, 
            fineMesh.E, fineMesh.F, WInvFine, PInvFine);
        fineMesh.V = S_0 * coarseMesh.V;
        // Compute edge topology
        fineMesh.compute_edge_topology_fixed_E();

        // Compute FEM operators
        coarseMesh.FEM_suite(coarseOperators);
        // Compute fine FEM operators
        fineMesh.FEM_suite(fineOperators);

        directional::get_P(fineMesh.V, fineMesh.F, fineMesh.E, fineMesh.FE, 1, PFine);
        directional::get_P_inverse(coarseMesh.V, coarseMesh.F, coarseMesh.E, coarseMesh.FE, 1, PInvCoarse);

        SECTION("S0_SGamma_Gv_Commutation["+testCaseFiles[i]+"]")
        {
            auto lhs = PFine * fineOperators.Gv * S_0;
            auto rhs = S_gamma * PCoarse * coarseOperators.Gv;
            REQUIRE(lhs.rows() == rhs.rows());
            REQUIRE(lhs.cols() == rhs.cols());
            // Difference should be zero
            Eigen::SparseMatrix<double> diff = lhs - rhs;
            REQUIRE(diff.coeffs().maxCoeff() == Approx(0).margin(1e-10));
        }
        SECTION("S1_S2_d1_Commutation[" + testCaseFiles[i] + "]")
        {
            Eigen::SparseMatrix<double> d1Coarse, d1Fine;
            dec_D1(coarseMesh, d1Coarse);
            dec_D1(fineMesh, d1Fine);

            auto lhs = S_2 * d1Coarse;
            auto rhs = d1Fine * S_1;
            REQUIRE(lhs.rows() == rhs.rows());
            REQUIRE(lhs.cols() == rhs.cols());
            // Difference should be zero
            Eigen::SparseMatrix<double> diff = lhs - rhs;
            REQUIRE(diff.coeffs().maxCoeff() == Approx(0).margin(1e-10));
        }
        SECTION("Sepsstar_SGamma_C_Commutation[" + testCaseFiles[i] + "]")
        {
            auto lhs = fineOperators.C * PInvFine * S_gamma;
            auto rhs = S_epsstar * coarseOperators.C * PInvCoarse;
            REQUIRE(S_gamma.coeffs().maxCoeff() > 0);
            REQUIRE(S_epsstar.coeffs().maxCoeff() > 0);
            REQUIRE(lhs.rows() == rhs.rows());
            REQUIRE(lhs.cols() == rhs.cols());
            // Difference should be zero
            Eigen::SparseMatrix<double> diff = lhs - rhs;
            REQUIRE(diff.coeffs().maxCoeff() == Approx(0).margin(1e-10));
        }
        SECTION("Sepsstar_S2_A_E_To_F_Commutation[" + testCaseFiles[i] + "]")
        {
            Eigen::SparseMatrix<double> A_E_To_F_fine, A_E_To_F_coarse;
            adjacencyMatrix(coarseMesh.FE, A_E_To_F_coarse);
            adjacencyMatrix(fineMesh.FE, A_E_To_F_fine);

            auto lhs = S_2 * A_E_To_F_coarse;
            auto rhs = A_E_To_F_fine * S_epsstar;
            REQUIRE(lhs.rows() == rhs.rows());
            REQUIRE(lhs.cols() == rhs.cols());
            // Difference should be zero
            Eigen::SparseMatrix<double> diff = lhs - rhs;
            REQUIRE(diff.coeffs().maxCoeff() == Approx(0).margin(1e-10));
        }
    }
}