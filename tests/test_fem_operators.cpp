#include "catch2/catch.hpp"
#include "directional/get_directional_subdivision_suite.h"
#include <igl/read_triangle_mesh.h>
#include "test_includes.h"
#include <igl/adjacency_matrix.h>

std::vector<std::string> femOpTestCaseFiles = { "bimba.off","cathead.off","chipped-torus.obj","half-torus.obj","horser.off","tester-sphere.off","torus.obj","bunny2.off" };

TEST_CASE("FEM operators", "[fem_operators]")
{
    for (int i = 0; i < femOpTestCaseFiles.size(); ++i)
    {
        TriangleMesh coarseMesh, fineMesh;
        // Operators that convert between PCVF and Gamma, and Gamma and mean-curl decomposition.
        Eigen::SparseMatrix<double> WCoarse, WInvCoarse, PCoarse, PInvCoarse;

        FEM_operators coarseOperators;

        // Load mesh via IGL
        coarseMesh.read(femOpTestCaseFiles[i]);
        coarseMesh.compute_edge_topology();
        // Compute FEM operators
        coarseMesh.FEM_suite(coarseOperators);

        directional::get_P(coarseMesh.V, coarseMesh.F, coarseMesh.E, coarseMesh.FE, 1, PCoarse);
        directional::get_P_inverse(coarseMesh.V, coarseMesh.F, coarseMesh.E, coarseMesh.FE, 1, PInvCoarse);
        directional::get_W(coarseMesh.F, coarseMesh.E, coarseMesh.FE, coarseMesh.EF, WCoarse);
        directional::get_W_inverse(coarseMesh.F, coarseMesh.E, coarseMesh.FE, coarseMesh.EF, WInvCoarse);

        SECTION("C Grad = 0 operator[" + femOpTestCaseFiles[i] + "]")
        {
            Eigen::SparseMatrix<double> result = coarseOperators.C * coarseOperators.Gv;

            REQUIRE(result.coeffs().maxCoeff() == Approx(0).margin(1e-10));
        }
        SECTION("Div J Grad_e = 0 operator[" + femOpTestCaseFiles[i] + "]")
        {
            Eigen::SparseMatrix<double> result = coarseOperators.D * coarseOperators.J * coarseOperators.Ge;

            REQUIRE(result.coeffs().maxCoeff() == Approx(0).margin(1e-10));
        }
        SECTION("P P^-1 = Id[" + femOpTestCaseFiles[i] + "]")
        {
            Eigen::SparseMatrix<double> prod = PCoarse * PInvCoarse;
            Eigen::SparseMatrix<double> id(prod.rows(), prod.cols());
            id.setIdentity();

            Eigen::SparseMatrix<double> diff = prod - id;

            REQUIRE(diff.coeffs().maxCoeff() == Approx(0).margin(1e-10));
        }
        SECTION("P^-1 P  = Id for PCVFs[" + femOpTestCaseFiles[i] + "]")
        {
            Eigen::MatrixXd basis;
            coarseMesh.basis(basis);
            // Make sure the edge vectors are columns
            basis.transposeInPlace();
            // Create an expansion operator that takes per edge 2 coefficients and returns a PCVF 
            // using the basis.
            Eigen::SparseMatrix<double> basisExpander;
            toBlocks(basis, 2, 2, basisExpander);

            Eigen::SparseMatrix<double> prod = basisExpander.transpose() * PInvCoarse * PCoarse * basisExpander;
            Eigen::SparseMatrix<double> id(prod.rows(), prod.cols());
            id.setIdentity();

            Eigen::SparseMatrix<double> diff = prod - id;

            REQUIRE(diff.coeffs().maxCoeff() == Approx(0).margin(1e-10));
        }
        SECTION("W^-1 W = Id[" + femOpTestCaseFiles[i] + "]")
        {
            Eigen::SparseMatrix<double> prod = WInvCoarse * WCoarse;
            Eigen::SparseMatrix<double> id(prod.rows(), prod.cols());
            id.setIdentity();

            Eigen::SparseMatrix<double> diff = prod - id;

            REQUIRE(diff.coeffs().maxCoeff() == Approx(0).margin(1e-10));
        }
       
    }
}