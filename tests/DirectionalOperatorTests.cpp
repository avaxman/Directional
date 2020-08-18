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

//TODO!!! Incorporate below into testcases
//{// Debug check partial result
//Eigen::VectorXd partial = G2_To_Decomp_0 * columnDirectional_To_G2 * columnDirectional;
//// Verify that nullsum is correct
//Eigen::SparseMatrix<double> d1, avgToFace;
//directional::Matched_A_C_To_F(SFE, EF.rows(), N, matching, avgToFace);
//directional::Matched_D1(SFE, EF.rows(), N, matching, d1);
//auto diff = d1 * partial.segment(0, partial.size() / 2) - avgToFace * partial.segment(partial.size() / 2, partial.size() / 2);
//std::cout << "Null sum coarse max diff" << diff.maxCoeff() << std::endl;
//        }
//
//{// Debug check partial result
//    Eigen::VectorXd partial = S_Decomp * G2_To_Decomp_0 * columnDirectional_To_G2 * columnDirectional;
//    // Verify that nullsum is correct
//    Eigen::SparseMatrix<double> d1, avgToFace;
//    directional::Matched_A_C_To_F(SFEK, EFK.rows(), N, matchingK, avgToFace);
//    directional::Matched_D1(SFEK, EFK.rows(), N, matchingK, d1);
//    auto diff = d1 * partial.segment(0, partial.size() / 2) - avgToFace * partial.segment(partial.size() / 2, partial.size() / 2);
//    std::cout << "Null sum fine max diff" << diff.maxCoeff() << std::endl;
//    std::cout << "Count diff:" << diff.size() << ", F * N:" << SFEK.rows() * N << std::endl;
//    failedNullsumFaces.setConstant(diff.size(), 0);
//    for (int i = 0; i < diff.size(); ++i)
//    {
//        failedNullsumFaces(i) = diff(i) > 1e-6 ? 1 : 0;
//    }
//}
//S_Gamma_directional = Decomp_To_G2K * S_Decomp*G2_To_Decomp_0;
//
//
//
//fineDirectional = Matched_Gamma2_To_PCVF_K * S_Gamma_directional * columnDirectional_To_G2 * columnDirectional;
//
//columndirectional_to_rawfield(fineDirectional, N, rawFieldK);
//for (int i = 0; i < N; ++i)
//{
//    Eigen::MatrixXd block = rawField.block(0, 3 * i, rawField.rows(), 3);
//    std::cout << "Coarse max norm for " << i << ":" << block.rowwise().norm().maxCoeff() << std::endl;
//}
//for (int i = 0; i < N; ++i)
//{
//    Eigen::MatrixXd block = rawFieldK.block(0, 3 * i, rawFieldK.rows(), 3);
//    std::cout << "Fine max norm for " << i << ":" << block.rowwise().norm().maxCoeff() << std::endl;
//}
//{
//    Eigen::VectorXd testField;
//    rawfield_to_columndirectional(rawField, N, testField);
//    Eigen::SparseMatrix<double> Matched_Gamma2_To_PCVF_0;
//    { //Construct Gamma2 to PCVF for directionals
//        Eigen::SparseMatrix<double> Gamma2_To_PCVF_0;
//        directional::Gamma2_reprojector(V, F, E, SFE, EF, Gamma2_To_PCVF_0);
//        std::vector<Eigen::SparseMatrix<double>*> base(N, &Gamma2_To_PCVF_0);
//        directional::block_diag(base, Matched_Gamma2_To_PCVF_0);
//    }
//    {
//        Eigen::MatrixXd testOut;
//        columndirectional_to_rawfield(testField, N, testOut);
//        Eigen::MatrixXd diff = testOut - rawField;
//        std::cout << "Columndirectional conversion error : " << diff.maxCoeff() << std::endl;
//    }
//    {
//
//        Eigen::MatrixXd testOut;
//        columndirectional_to_rawfield(Matched_Gamma2_To_PCVF_0 * columnDirectional_To_G2 * testField, N, testOut);
//        Eigen::MatrixXd diff2 = testOut - rawField;
//        std::cout << "Conversion via G2 and column direction conversion error : " << diff2.maxCoeff() << std::endl;
//    }
//    {
//
//        Eigen::MatrixXd testOut;
//        Eigen::SparseMatrix<double> Decomp_To_G20;
//        directional::Matched_AC_To_Gamma2(EF, SFE, EI, matching, N, Decomp_To_G20);
//        columndirectional_to_rawfield(Matched_Gamma2_To_PCVF_0 * Decomp_To_G20 * G2_To_Decomp_0 * columnDirectional_To_G2 * testField, N, testOut);
//        Eigen::MatrixXd diff2 = testOut - rawField;
//        std::cout << "Conversion via AC, G2 and column direction conversion error : " << diff2.maxCoeff() << std::endl;
//    }
//    {
//        Eigen::VectorXd testFieldK;
//        rawfield_to_columndirectional(rawFieldK, N, testFieldK);
//        Eigen::MatrixXd testOut;
//        Eigen::SparseMatrix<double> G2_To_Decomp_K, columnDirectional_To_G2K;
//        directional::Matched_Gamma2_To_AC(EIK, EFK, SFEK, matchingK, N, G2_To_Decomp_K);
//        directional::columndirectional_to_gamma2_matrix(Vk, Fk, EVK, SFEK, EFK, N, columnDirectional_To_G2K);
//
//        columndirectional_to_rawfield(Matched_Gamma2_To_PCVF_K * Decomp_To_G2K * G2_To_Decomp_K * columnDirectional_To_G2K * testFieldK, N, testOut);
//        Eigen::MatrixXd diff2 = testOut - rawFieldK;
//        std::cout << "Conversion via AC, G2 and column direction conversion error in fine: " << diff2.maxCoeff() << std::endl;
//    }
//    /*{
//        Eigen::SparseMatrix<double> Decomp_To_G20;
//        Eigen::SparseMatrix<double> Decomp_To_G2PCVF, DCToG3, PCVFG2ToG3;
//        directional::Matched_AC_To_Gamma2(EF, SFE, EI, matchingK, 1, Decomp_To_G20);
//        directional::Decomp_To_Gamma3(EF, EI, F.rows(), DCToG3);
//        directional::Gamma2_To_Gamma3(SFE, PCVFG2ToG3);
//        auto diff =
//    }*/
//}


//double tolerance = 1e-7;
//
//using SparseMat = Eigen::SparseMatrix<double>;
//
//std::string description(const Eigen::SparseMatrix<double>& mat)
//{
//	std::stringstream str;
//	str << "Cols:" << mat.cols() << ", rows: " << mat.rows();
//	return str.str();
//}
//
//using namespace directional_fixtures;
//
///**
// * Test for verifying that the DCEL works properly
// */
//TEST_P(DCELTestFixture, DCELTest)
//{
//	// Test the initial 
//	DCEL dcel(ED);
//	dcel.iterateRings(*this);
//	// Check all nodes and edges are visited
//	this->checkVisit();
//
//	EdgeData EDL;
//	Eigen::MatrixXi E0ToEk;
//	ED.quadrisect(m_mesh.V.rows(), E0ToEk, EDL);
//	// Replace the edge data locally
//
//	EdgeData localED = ED;
//	this->ED = EDL;
//
//	// Run again on quadrisected stuff
//	DCEL dcel2(EDL);
//	dcel2.iterateRings(*this);
//	checkVisit();
//
//	this->ED = localED;
//
//	// Check that boundary is preserved
//	Eigen::VectorXi vIsBoundary;
//	ED.boundaryVerts(vIsBoundary);
//	Eigen::VectorXi vKIsBoundary;
//	EDL.boundaryVerts(vKIsBoundary);
//	for(int e = 0; e < E0ToEk.rows(); e++)
//	{
//		const int bCount = vIsBoundary(ED.E(e, 0)) + vIsBoundary(ED.E(e, 1));
//		const int sBound = vIsBoundary(ED.E(e, 0));
//		const int eBound = vIsBoundary(ED.E(e, 1));
//		const int lBound = ED.EF(e, 0) == -1 ? -1 : vIsBoundary(ED.F(ED.EF(e,0),ED.EI(e, 0)));
//		const int rBound = ED.EF(e, 1) == -1 ? -1 : vIsBoundary(ED.F(ED.EF(e, 1), ED.EI(e, 1)));
//
//		// Iterate over all subdivided edge types: per edge flap, the initial edge is split into a ''start'' and ''end'' edge
//		// In addition, the left and right faces are subdivided, giving the ''left'' and ''right'' edges that are ''parallel'' to the
//		// initial edge
//		for(int j = 0; j < 4; j++)
//		{
//			if (E0ToEk(e, j) != -1) {
//				const int eIn = E0ToEk(e, j);
//				int bCountInner = vKIsBoundary(EDL.E(eIn, 0)) + vKIsBoundary(EDL.E(eIn, 1));
//				const int fullBound = lBound == -1 || rBound == -1;
//				switch(j)
//				{
//					// ''Start'' edge
//				case 0:
//					EXPECT_TRUE(bCountInner == sBound + fullBound) << "Invalid boundary at start edge: " << bCountInner << " vs " << sBound + fullBound;
//					break;
//					// ''End'' edge
//				case 1:
//					EXPECT_TRUE(bCountInner == eBound + fullBound) << "Invalid boundary at end edge: " << bCountInner << " vs " << eBound + fullBound;
//					break;
//					// ''Left'' edge
//				case 2:
//					EXPECT_TRUE(bCountInner == (sBound & lBound) + (eBound & lBound)) << "Invalid boundary at left edge: " << bCountInner << " vs " << (sBound & lBound) + (eBound & lBound);
//					break;
//					// ''Right'' edge
//				case 3:
//					EXPECT_TRUE(bCountInner == (sBound & rBound) + (eBound & rBound)) << "Invalid boundary at right edge: " << bCountInner << " vs " << (sBound & rBound) + (eBound & rBound);
//					break;
//				default:
//					break;
//				}
//			}
//		}
//	}
//}
//
//TEST_P(MeshTestFixture, DECTests)
//{
//	using SMat = Eigen::SparseMatrix<double>;
//	EdgeData ED(m_mesh.F);
//	SMat d0, d1;
//	directional::DEC_d0(m_mesh.V, m_mesh.F, ED.E, ED.sFE, ED.EF, d0);
//	directional::DEC_d1(m_mesh.V, m_mesh.F, ED.E, ED.sFE, ED.EF, d1);
//	SMat prod = d1 * d0;
//	auto els = helpers::getDifferenceFromZero(prod, 1e-7, 100);
//	EXPECT_TRUE(els.size() == 0) << std::string("d1 d0 is not zero") + helpers::tripletsToString(els);
//}
//
//TEST_P(MeshTestFixture, GammaInverseTest)
//{
//	using SMat = Eigen::SparseMatrix<double>;
//
//	// Construct gamma operators
//	helpers::Gamma2_Ops ops_Gamma(m_mesh);
//
//	SMat possiblyId = ops_Gamma.Chi_To_Gamma2 * ops_Gamma.Gamma2_To_Chi;
//	auto diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
//	EXPECT_TRUE(diffs.size() == 0) << std::string("P P^Inv is not identity for Gamma2:") + helpers::tripletsToString(diffs);
//	// Note that P^Inv P is not identity by definition but a normal component remover.
//}
//
//TEST_P(MeshTestFixture, G2G3ConversionTest)
//{
//	using SMat = Eigen::SparseMatrix<double>;
//
//	// Construct gamma operators
//	helpers::Gamma2_Ops ops_Gamma(m_mesh);
//	SMat G2_To_G3, G3_To_G2;
//	EdgeData ED = m_mesh.getED();
//	directional::Gamma2_To_Gamma3(ED.sFE, G2_To_G3);
//	directional::Gamma3_To_Gamma2(ED.sFE, G3_To_G2);
//
//	/*SMat possiblyId = G2_To_G3 * G3_To_G2;
//	std::vector<std::tuple<int, int, double>> diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
//	BOOST_TEST(diffs.size() == 0, std::string("G3->G2->G3 is not identity for Gamma2:") + helpers::tripletsToString(diffs));*/
//
//	SMat possiblyId = G3_To_G2 * G2_To_G3;
//	std::vector < std::tuple<int, int, double>> diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
//	EXPECT_EQ(diffs.size(), 0) << std::string("G2->G3->G2 is not identity for Decomp:") + helpers::tripletsToString(diffs);
//}
//
//
//TEST_P(MeshTestFixture, DecompTest)
//{
//	using SMat = Eigen::SparseMatrix<double>;
//
//	// Construct gamma operators
//	helpers::Gamma2_Ops ops_Gamma(m_mesh);
//
//	EdgeData ED(m_mesh.F);
//	SMat G3To1Form, C, DCToG3, G3ToDC;
//	directional::Gamma3_To_Decomp(ED.EF, ED.EI, ED.sFE.rows(), G3ToDC);
//	directional::Decomp_To_Gamma3(ED.EF, ED.EI, ED.sFE.rows(), DCToG3);
//
//	SMat possiblyId = ops_Gamma.G3ToG2 * ops_Gamma.G2ToG3;
//	std::vector < std::tuple<int, int, double>> diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
//	EXPECT_TRUE(diffs.size() == 0) << std::string("G2->G3->G2 is not identity for Decomp:") + helpers::tripletsToString(diffs);
//
//	possiblyId = ops_Gamma.G3ToG2 * DCToG3 * G3ToDC * ops_Gamma.G2ToG3;
//	diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
//	EXPECT_TRUE(diffs.size() == 0) << std::string("G2->Decomp->G2 is not identity for Gamma2:") + helpers::tripletsToString(diffs);
//
//	possiblyId = ops_Gamma.Decomp_To_G2 * ops_Gamma.G2_To_Decomp;
//	diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
//	EXPECT_TRUE(diffs.size() == 0) << std::string("G2->Decomp->G2 from ops struct is not identity for Gamma2:") + helpers::tripletsToString(diffs);
//
//	// Note that this is identity up to some non-null sum removing factors in hte matrix.
//	/*possiblyId = G3ToDC * ops_Gamma.G2ToG3 * ops_Gamma.G3ToG2*DCToG3 ;
//	diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
//	BOOST_TEST(diffs.size() == 0, std::string("Decomp->G2->Decomp is not identity for Decomp:") + helpers::tripletsToString(diffs));*/
//}
//
//TEST_P(MeshTestFixture, Decomp3Test)
//{
//	using SMat = Eigen::SparseMatrix<double>;
//
//	EdgeData ED = m_mesh.getED();
//
//	// Construct gamma operators
//	SMat G3To1Form, C, DCToG3, G3ToDC;
//	directional::Gamma3_To_Decomp(ED.EF, ED.EI, ED.sFE.rows(), G3ToDC);
//	//igl::cat(2, SMat(2.0 * G3To1Form.transpose()), SMat(C.transpose()), DCToG3);
//	directional::Decomp_To_Gamma3(ED.EF, ED.EI, ED.sFE.rows(), DCToG3);
//
//
//	SMat possiblyId = DCToG3 * G3ToDC;
//	std::vector<std::tuple<int, int, double>> diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
//	EXPECT_TRUE(diffs.size() == 0) << std::string("G3->Decomp->G3 is not identity for Gamma2:") + helpers::tripletsToString(diffs);
//
//	possiblyId = G3ToDC * DCToG3;
//	diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
//	EXPECT_TRUE(diffs.size() == 0) << std::string("Decomp->G3->Decomp is not identity for Decomp:") + helpers::tripletsToString(diffs);
//}
//
//TEST_P(MeshTestFixture, OperatorTests)
//{
//	using SMat = Eigen::SparseMatrix<double>;
//
//	// Construct gamma operators
//	helpers::Gamma2_Ops ops_Gamma(m_mesh);
//	helpers::FEM_Ops ops_FEM(m_mesh);
//
//	ASSERT_TRUE(ops_Gamma.Chi_To_Gamma2.cols() == ops_FEM.D.cols()) << "Operators Chi to Gamma2 and Divergence of chi do not have same column number";
//	ASSERT_TRUE(ops_Gamma.D.rows() == ops_FEM.D.rows()) << "Operators D_Gamma and D_Chi to not map to same space";
//
//	SMat diff = (ops_Gamma.D * ops_Gamma.Chi_To_Gamma2 - ops_FEM.D);
//	double val = helpers::maxAbs(diff);
//	EXPECT_TRUE(val < tolerance) << std::string("D_Gamma * P differs from D_Chi: ") +std::to_string(val);
//	SMat diff2 = (ops_Gamma.D  - ops_FEM.D * ops_Gamma.Gamma2_To_Chi);
//	val = helpers::maxAbs(diff2);
//	EXPECT_TRUE(val < tolerance) << std::string("D_Gamma differs from D_Chi * PInv: ") + std::to_string(val);
//}
//
//INSTANTIATE_TEST_SUITE_P(OperatorTests, MeshTestFixture,
//	testing::ValuesIn(SIMPLE_FILES));
//
//int main(int argc, char **argv) {
//
//	::testing::InitGoogleTest(&argc, argv);
//	return RUN_ALL_TESTS();
//}