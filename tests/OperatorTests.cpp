//#define BOOST_TEST_MODULE Operators
//#include <boost/test/included/unit_test.hpp>
//#include <boost/test/data/test_case.hpp>
//#include <boost/test/data/monomorphic.hpp>
#include "DCELFixture.h"
#include "Helpers.h"
#include <Eigen/Eigen>
#include <gtest/gtest.h>
#include <igl/read_triangle_mesh.h>
#include <directional/Subdivision/subdivision.h>
#include <directional/Subdivision/DCEL.h>
#include <igl/max.h>
#include <igl/cat.h>
#define SIMPLE_FILES {"bimba.off", "horser.off", "chipped-torus.obj"}

#include <directional/FEM_suite.h>
#include <directional/Gamma_suite.h>

double tolerance = 1e-7;

using SparseMat = Eigen::SparseMatrix<double>;

std::string description(const Eigen::SparseMatrix<double>& mat)
{
	std::stringstream str;
	str << "Cols:" << mat.cols() << ", rows: " << mat.rows();
	return str.str();
}

using namespace directional_fixtures;

/**
 * Test for verifying that the DCEL works properly
 */
TEST_P(DCELTestFixture, DCELTest)
{
	// Test the initial 
	DCEL dcel(ED);
	dcel.iterateRings(*this);
	// Check all nodes and edges are visited
	this->checkVisit();

	EdgeData EDL;
	Eigen::MatrixXi E0ToEk;
	ED.quadrisect(m_mesh.V.rows(), E0ToEk, EDL);
	// Replace the edge data locally

	EdgeData localED = ED;
	this->ED = EDL;

	// Run again on quadrisected stuff
	DCEL dcel2(EDL);
	dcel2.iterateRings(*this);
	checkVisit();

	this->ED = localED;

	// Check that boundary is preserved
	Eigen::VectorXi vIsBoundary;
	ED.boundaryVerts(vIsBoundary);
	Eigen::VectorXi vKIsBoundary;
	EDL.boundaryVerts(vKIsBoundary);
	for(int e = 0; e < E0ToEk.rows(); e++)
	{
		const int bCount = vIsBoundary(ED.E(e, 0)) + vIsBoundary(ED.E(e, 1));
		const int sBound = vIsBoundary(ED.E(e, 0));
		const int eBound = vIsBoundary(ED.E(e, 1));
		const int lBound = ED.EF(e, 0) == -1 ? -1 : vIsBoundary(ED.F(ED.EF(e,0),ED.EI(e, 0)));
		const int rBound = ED.EF(e, 1) == -1 ? -1 : vIsBoundary(ED.F(ED.EF(e, 1), ED.EI(e, 1)));

		// Iterate over all subdivided edge types: per edge flap, the initial edge is split into a ''start'' and ''end'' edge
		// In addition, the left and right faces are subdivided, giving the ''left'' and ''right'' edges that are ''parallel'' to the
		// initial edge
		for(int j = 0; j < 4; j++)
		{
			if (E0ToEk(e, j) != -1) {
				const int eIn = E0ToEk(e, j);
				int bCountInner = vKIsBoundary(EDL.E(eIn, 0)) + vKIsBoundary(EDL.E(eIn, 1));
				const int fullBound = lBound == -1 || rBound == -1;
				switch(j)
				{
					// ''Start'' edge
				case 0:
					EXPECT_TRUE(bCountInner == sBound + fullBound) << "Invalid boundary at start edge: " << bCountInner << " vs " << sBound + fullBound;
					break;
					// ''End'' edge
				case 1:
					EXPECT_TRUE(bCountInner == eBound + fullBound) << "Invalid boundary at end edge: " << bCountInner << " vs " << eBound + fullBound;
					break;
					// ''Left'' edge
				case 2:
					EXPECT_TRUE(bCountInner == (sBound & lBound) + (eBound & lBound)) << "Invalid boundary at left edge: " << bCountInner << " vs " << (sBound & lBound) + (eBound & lBound);
					break;
					// ''Right'' edge
				case 3:
					EXPECT_TRUE(bCountInner == (sBound & rBound) + (eBound & rBound)) << "Invalid boundary at right edge: " << bCountInner << " vs " << (sBound & rBound) + (eBound & rBound);
					break;
				default:
					break;
				}
			}
		}
	}
}

TEST_P(MeshTestFixture, DECTests)
{
	using SMat = Eigen::SparseMatrix<double>;
	EdgeData ED(m_mesh.F);
	SMat d0, d1;
	directional::DEC_d0(m_mesh.V, m_mesh.F, ED.E, ED.sFE, ED.EF, d0);
	directional::DEC_d1(m_mesh.V, m_mesh.F, ED.E, ED.sFE, ED.EF, d1);
	SMat prod = d1 * d0;
	auto els = helpers::getDifferenceFromZero(prod, 1e-7, 100);
	EXPECT_TRUE(els.size() == 0) << std::string("d1 d0 is not zero") + helpers::tripletsToString(els);
}

TEST_P(MeshTestFixture, GammaInverseTest)
{
	using SMat = Eigen::SparseMatrix<double>;

	// Construct gamma operators
	helpers::Gamma2_Ops ops_Gamma(m_mesh);

	SMat possiblyId = ops_Gamma.Chi_To_Gamma2 * ops_Gamma.Gamma2_To_Chi;
	auto diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	EXPECT_TRUE(diffs.size() == 0) << std::string("P P^Inv is not identity for Gamma2:") + helpers::tripletsToString(diffs);
	// Note that P^Inv P is not identity by definition but a normal component remover.
}

TEST_P(MeshTestFixture, G2G3ConversionTest)
{
	using SMat = Eigen::SparseMatrix<double>;

	// Construct gamma operators
	helpers::Gamma2_Ops ops_Gamma(m_mesh);
	SMat G2_To_G3, G3_To_G2;
	EdgeData ED = m_mesh.getED();
	directional::Gamma2_To_Gamma3(ED.sFE, G2_To_G3);
	directional::Gamma3_To_Gamma2(ED.sFE, G3_To_G2);

	/*SMat possiblyId = G2_To_G3 * G3_To_G2;
	std::vector<std::tuple<int, int, double>> diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	BOOST_TEST(diffs.size() == 0, std::string("G3->G2->G3 is not identity for Gamma2:") + helpers::tripletsToString(diffs));*/

	SMat possiblyId = G3_To_G2 * G2_To_G3;
	std::vector < std::tuple<int, int, double>> diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	EXPECT_EQ(diffs.size(), 0) << std::string("G2->G3->G2 is not identity for Decomp:") + helpers::tripletsToString(diffs);
}


TEST_P(MeshTestFixture, DecompTest)
{
	using SMat = Eigen::SparseMatrix<double>;

	// Construct gamma operators
	helpers::Gamma2_Ops ops_Gamma(m_mesh);

	EdgeData ED(m_mesh.F);
	SMat G3To1Form, C, DCToG3, G3ToDC;
	directional::Gamma3_To_Decomp(ED.EF, ED.EI, ED.sFE.rows(), G3ToDC);
	directional::Decomp_To_Gamma3(ED.EF, ED.EI, ED.sFE.rows(), DCToG3);

	SMat possiblyId = ops_Gamma.G3ToG2 * ops_Gamma.G2ToG3;
	std::vector < std::tuple<int, int, double>> diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	EXPECT_TRUE(diffs.size() == 0) << std::string("G2->G3->G2 is not identity for Decomp:") + helpers::tripletsToString(diffs);

	possiblyId = ops_Gamma.G3ToG2 * DCToG3 * G3ToDC * ops_Gamma.G2ToG3;
	diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	EXPECT_TRUE(diffs.size() == 0) << std::string("G2->Decomp->G2 is not identity for Gamma2:") + helpers::tripletsToString(diffs);

	possiblyId = ops_Gamma.Decomp_To_G2 * ops_Gamma.G2_To_Decomp;
	diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	EXPECT_TRUE(diffs.size() == 0) << std::string("G2->Decomp->G2 from ops struct is not identity for Gamma2:") + helpers::tripletsToString(diffs);

	// Note that this is identity up to some non-null sum removing factors in hte matrix.
	/*possiblyId = G3ToDC * ops_Gamma.G2ToG3 * ops_Gamma.G3ToG2*DCToG3 ;
	diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	BOOST_TEST(diffs.size() == 0, std::string("Decomp->G2->Decomp is not identity for Decomp:") + helpers::tripletsToString(diffs));*/
}

TEST_P(MeshTestFixture, Decomp3Test)
{
	using SMat = Eigen::SparseMatrix<double>;

	EdgeData ED = m_mesh.getED();

	// Construct gamma operators
	SMat G3To1Form, C, DCToG3, G3ToDC;
	directional::Gamma3_To_Decomp(ED.EF, ED.EI, ED.sFE.rows(), G3ToDC);
	//igl::cat(2, SMat(2.0 * G3To1Form.transpose()), SMat(C.transpose()), DCToG3);
	directional::Decomp_To_Gamma3(ED.EF, ED.EI, ED.sFE.rows(), DCToG3);


	SMat possiblyId = DCToG3 * G3ToDC;
	std::vector<std::tuple<int, int, double>> diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	EXPECT_TRUE(diffs.size() == 0) << std::string("G3->Decomp->G3 is not identity for Gamma2:") + helpers::tripletsToString(diffs);

	possiblyId = G3ToDC * DCToG3;
	diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	EXPECT_TRUE(diffs.size() == 0) << std::string("Decomp->G3->Decomp is not identity for Decomp:") + helpers::tripletsToString(diffs);
}

TEST_P(MeshTestFixture, OperatorTests)
{
	using SMat = Eigen::SparseMatrix<double>;

	// Construct gamma operators
	helpers::Gamma2_Ops ops_Gamma(m_mesh);
	helpers::FEM_Ops ops_FEM(m_mesh);

	ASSERT_TRUE(ops_Gamma.Chi_To_Gamma2.cols() == ops_FEM.D.cols()) << "Operators Chi to Gamma2 and Divergence of chi do not have same column number";
	ASSERT_TRUE(ops_Gamma.D.rows() == ops_FEM.D.rows()) << "Operators D_Gamma and D_Chi to not map to same space";

	SMat diff = (ops_Gamma.D * ops_Gamma.Chi_To_Gamma2 - ops_FEM.D);
	double val = helpers::maxAbs(diff);
	EXPECT_TRUE(val < tolerance) << std::string("D_Gamma * P differs from D_Chi: ") +std::to_string(val);
	SMat diff2 = (ops_Gamma.D  - ops_FEM.D * ops_Gamma.Gamma2_To_Chi);
	val = helpers::maxAbs(diff2);
	EXPECT_TRUE(val < tolerance) << std::string("D_Gamma differs from D_Chi * PInv: ") + std::to_string(val);
}

INSTANTIATE_TEST_SUITE_P(OperatorTests, MeshTestFixture,
	testing::ValuesIn(SIMPLE_FILES));

int main(int argc, char **argv) {

	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}