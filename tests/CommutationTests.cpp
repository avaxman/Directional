// Include this before anything else
#include "Helpers.h"
#include <directional/Subdivision/subdivision.h>
#include <igl/max.h>
#include <directional/Gamma_suite.h>
#include <directional/block_diag.h>
#define SIMPLE_FILES {"bimba.off", "horser.off","chipped-torus.obj"}


// Type aliases
using SparseMat = Eigen::SparseMatrix<double>;

/**
 *		TEST GLOBAL SETTINGS
 */
double tolerance = 1e-7;
int level = 1;

struct SubData
{
	using SMat = Eigen::SparseMatrix<double>;
	EdgeData ED0, EDL;
	SMat S_v, S_e, S_f, S_c;
	Eigen::MatrixXd VK;
	SubData(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int level)
	{
		subdivision_operators(F, V, level, S_v, S_f, S_e, S_c, ED0, EDL);
		VK = S_v * V;
	}
};
using namespace directional_fixtures;

// Tests whether the quadrisection step is correct or not.
TEST_P(MeshTestFixture, QuadrisectionCorrectness)
{
	EdgeData ED0, ED1;
	ED0.construct(m_mesh.F);

	ASSERT_TRUE(ED0.isConsistent(true));
	Eigen::MatrixXi E0ToEk;
	ED0.quadrisect(m_mesh.V.rows(), E0ToEk, ED1);
	
	// Check new element counts
	EXPECT_EQ(ED1.faceCount(), 4 * ED0.faceCount());
	EXPECT_EQ(ED1.edgeCount(), 2 * ED0.edgeCount() + 3 * ED0.faceCount());
	// Require consistency for the quadrisected element
	ASSERT_TRUE(ED1.isConsistent(true));

	// Count boundary edges
	int bCount = 0;
	for (int e = 0; e < ED1.EF.rows(); e++)
	{
		if (ED1.EF(e, 0) == -1 || ED1.EF(e, 1) == -1) bCount++;
	}
	// Check that boundary edges is correctly recorded and consistent with pre-quadrisection
	EXPECT_EQ(bCount, ED1.boundaryEdgeCount);
	EXPECT_EQ(2 * ED0.boundaryEdgeCount, ED1.boundaryEdgeCount);

	// Check that sFE is correct.
	for (int f = 0; f < ED0.sFE.rows(); f++)
	{
		for (int j = 0; j < 3; j++)
		{
			EXPECT_EQ(ED0.F(f, j), ED1.F(4 * f + j, j));

			const int e = ED0.sFE(f, j);
			const int c = ED0.sFE(f, j + 3);

			// Check the edge mapping
			EXPECT_EQ(ED1.E(E0ToEk(e, 0), 0), ED0.E(e, 0));
			EXPECT_EQ(ED1.E(E0ToEk(e, 1), 1), ED0.E(e, 1));
			EXPECT_EQ(ED1.E(E0ToEk(e, 0), 1), ED0.vertexCount() + e);
			EXPECT_EQ(ED1.E(E0ToEk(e, 1), 0), ED0.vertexCount() + e);
		}
	}
}

TEST_P(MeshTestFixture, DEC_SEQUENCE)
{
	using SMat = Eigen::SparseMatrix<double>;
	EdgeData ED0, EDL;

	// Construct subdivision operators
	SMat S_V, S_F, S_E, S_C;
	subdivision_operators(m_mesh.F, m_mesh.V, level, S_V, S_F, S_E, S_C, ED0, EDL);

	// Test sizes
	ASSERT_EQ(S_E.rows(), EDL.edgeCount());
	ASSERT_EQ(S_C.rows(), EDL.edgeCount());
	ASSERT_EQ(S_E.cols(), ED0.edgeCount());
	ASSERT_EQ(S_C.cols(), ED0.edgeCount());

	directional_fixtures::Mesh meshK;
	meshK.F = EDL.F;
	// Subdivided upwards.
	meshK.V = S_V * m_mesh.V;

	// Require new level edge data to be consistent
	ASSERT_TRUE(EDL.isConsistent());

	// Test d1 d0 = 0 relation
	SparseMat d0_0, d0_K, d1_K;
	directional::DEC_d0(meshK.V, meshK.F, EDL.E, EDL.sFE, EDL.EF, d0_K);
	directional::DEC_d1(meshK.V, meshK.F, EDL.E, EDL.sFE, EDL.EF, d1_K);

	auto diffEls2 = helpers::getDifferenceFromZero(SparseMat(d1_K* d0_K), tolerance, 100);
	EXPECT_TRUE(diffEls2.size() == 0) << "d1_K d0_K is not zero" << helpers::tripletsToString(diffEls2);
}

TEST_P(MeshTestFixture, SEC_D0_Commutation)
{
	SubData sd(m_mesh.V, m_mesh.F, level);

	// Test basic d0 commutation
	SparseMat d0_0, d0_K;
	directional::DEC_d0(m_mesh.V, m_mesh.F, sd.ED0.E, sd.ED0.sFE, sd.ED0.EF, d0_0);
	directional::DEC_d0(sd.VK, sd.EDL.F, sd.EDL.E, sd.EDL.sFE, sd.EDL.EF, d0_K);
	// Commutation with d0
	SparseMat d0_p1 = sd.S_e * d0_0;
	SparseMat d0_p2 = d0_K * sd.S_v;

	// Retrieve differences
	std::vector<std::tuple<int, int, double, double>> diffEls = helpers::getAbsDifference(d0_p1, d0_p2, tolerance, 100);
	std::stringstream msg;
	if (diffEls.size() > 0)
	{
		Eigen::VectorXi BVs;
		Eigen::VectorXi BVsK;
		Eigen::VectorXi valence;
		sd.ED0.boundaryVerts(BVs);
		sd.ED0.vertexValence(valence);
		sd.EDL.boundaryVerts(BVsK);
		const int v0Count = sd.ED0.vertexCount();
		msg << std::endl;
		msg << "v0 count:" << v0Count << std::endl;
		for (int i = 0; i < diffEls.size(); i++)
		{
			const int r = std::get<0>(diffEls[i]), c = std::get<1>(diffEls[i]);;
			double v1 = std::get<2>(diffEls[i]), v2 = std::get<3>(diffEls[i]);

			const int v0K = sd.EDL.E(r, 0), v1K = sd.EDL.E(r, 1);

			msg << "El at " << r << ',' << c << "| src is bnd:" << BVs(c) << ", deg:"
				<< valence(c) << ", vals:" << v1 << ", " << v2 << ", E bnd:[" << BVsK(sd.EDL.E(r, 0)) << ", " << BVsK(sd.EDL.E(r, 1))
				<< "], EV:[" << sd.EDL.E(r, 0) << ',' << sd.EDL.E(r, 1) << "]"
				<< std::endl;
		}
	}
	EXPECT_TRUE(diffEls.size() == 0) << "d0 does not commute via S_e/S_v" << msg.str();
}

/**
 * Second SEC commutation
 */
TEST_P(MeshTestFixture, SEC_D1_Commutation)
{
	using SMat = Eigen::SparseMatrix<double>;
	int level = 1;

	SubData sd(m_mesh.V, m_mesh.F, level);

	// Test basic d0 commutation
	SparseMat d1_0, d1_K;
	directional::DEC_d1(m_mesh.V, m_mesh.F, sd.ED0.E, sd.ED0.sFE, sd.ED0.EF, d1_0);
	// Construct fine level d1s
	directional::DEC_d1(sd.VK, sd.EDL.F, sd.EDL.E, sd.EDL.sFE, sd.EDL.EF, d1_K);
	// To compare
	SparseMat d1_p1 = sd.S_f * d1_0;
	SparseMat d1_p2 = d1_K * sd.S_e;

	std::vector<std::tuple<int, int, double, double>> diffEls = helpers::getAbsDifference(d1_p1, d1_p2, tolerance, 100);
	EXPECT_TRUE(diffEls.size() == 0) << "d1 does not commute via S_e/S_f" << helpers::tripletsToString(diffEls);
}

// Curl commutation via the averager.
TEST_P(MeshTestFixture, SHM_Curl_Commutation)
{
	using SMat = Eigen::SparseMatrix<double>;
	int level = 1;

	SubData sd(m_mesh.V, m_mesh.F, level);
	SMat A_EToF0, A_EToFK;
	directional::Edge_To_Face_Average(sd.ED0.sFE, sd.ED0.edgeCount(), A_EToF0);
	directional::Edge_To_Face_Average(sd.EDL.sFE, sd.EDL.edgeCount(), A_EToFK);

	// To compare
	SparseMat lhs = sd.S_f * A_EToF0;
	SparseMat rhs = A_EToFK * sd.S_c;

	std::vector<std::tuple<int, int, double, double>> diffEls = helpers::getAbsDifference(lhs, rhs, tolerance, 100);
	std::stringstream msg;
	msg << std::endl;
	Eigen::VectorXi BVs;
	Eigen::VectorXi valence;
	sd.ED0.boundaryVerts(BVs);
	sd.ED0.vertexValence(valence);
	for(int i = 0; i < diffEls.size(); i++)
	{
		const int r = std::get<0>(diffEls[i]), c = std::get<1>(diffEls[i]);;
		double v1 = std::get<2>(diffEls[i]), v2 = std::get<3>(diffEls[i]);
		const int f = r / 4;
		const int corner = r - 4 * f;
		bool inOrigF = false;
		int ind = -1;
		for(int i = 0; i < 3; i++)
		{
			if(sd.ED0.sFE(f,i) == c)
			{
				inOrigF = true;
				ind = i;
				break;
			}
		}
		msg << "At [" << r << ',' << c << "], vals [" << v1 << "," << v2 << "], original corner:" << corner << ", in orig f: " << inOrigF <<"," << ind << " ";
		if(corner < 3)
		{
			msg << "Vert val:" << valence(sd.ED0.F(f, corner)) << ", is bnd: " << BVs(sd.ED0.F(f, corner));
		}
		msg << std::endl;
	}
	EXPECT_TRUE(diffEls.size() == 0) << "Average_E_To_F does not commute via S_c/S_f" << msg.str();
}

TEST_P(MeshTestFixture, SHM_Gradient_Commutation)
{
	using SMat = Eigen::SparseMatrix<double>;
	int level = 1;

	SubData sd(m_mesh.V, m_mesh.F, level);
	// Coarse level operators
	helpers::Gamma2_Ops ops0(sd.ED0, m_mesh.V);
	// Fine level operators
	helpers::Gamma2_Ops opsK(sd.EDL, sd.VK);

	SMat S_Gamma;
	directional::block_diag({ &sd.S_e, &sd.S_c }, S_Gamma);
	S_Gamma = opsK.Decomp_To_G2 * S_Gamma * ops0.G2_To_Decomp;


	// To compare
	SparseMat lhs = S_Gamma * ops0.Gv;
	SparseMat rhs = opsK.Gv * sd.S_v;

	std::vector<std::tuple<int, int, double, double>> diffEls = helpers::getAbsDifference(lhs, rhs, tolerance, 100);
	EXPECT_TRUE(diffEls.size() == 0) << "d1 does not commute via S_e/S_f" << helpers::tripletsToString(diffEls);
}

TEST_P(MeshTestFixture, SHM_SGamma_SE_Commutation)
{
	using SMat = Eigen::SparseMatrix<double>;

	SubData sd(m_mesh.V, m_mesh.F, level);
	// Coarse level operators
	helpers::Gamma2_Ops ops0(sd.ED0, m_mesh.V);
	// Fine level operators
	helpers::Gamma2_Ops opsK(sd.EDL, sd.VK);

	SMat S_Gamma;
	directional::block_diag({ &sd.S_e, &sd.S_c }, S_Gamma);
	S_Gamma = opsK.Decomp_To_G2 * S_Gamma * ops0.G2_To_Decomp;


	// To compare
	SparseMat lhs = sd.S_e * ops0.Gamma2_To_Oneform;
	SparseMat rhs = opsK.Gamma2_To_Oneform * S_Gamma;

	std::vector<std::tuple<int, int, double, double>> diffEls = helpers::getAbsDifference(lhs, rhs, tolerance, 100);
	EXPECT_TRUE(diffEls.size() == 0)<<"Avg_Gamma_To_E does not commute via S_e/S_Gamma:"<< helpers::tripletsToString(diffEls);
}

INSTANTIATE_TEST_SUITE_P(Commutations, MeshTestFixture,
	testing::ValuesIn(SIMPLE_FILES));

int main(int argc, char **argv) {
	
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}