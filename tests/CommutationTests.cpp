// Include this before anything else
#include "Helpers.h"
#include <directional/Subdivision/subdivision.h>
#include <igl/max.h>
#include <directional/Gamma_suite.h>
#include <directional/block_diag.h>
#define SIMPLE_FILES {"bimba.off", "horser.off","chipped-torus.obj","cheburashka-subdivision.off"}

// New subdivision code
#include <directional/Subdivision_new/shm_edge_topology.h>
#include <directional/Subdivision_new/is_edgedata_consistent.h>
#include <directional/Subdivision_new/quadrisect.h>
#include <directional/Subdivision_new/build_shm_subdivision.h>
#include <directional/Subdivision_new/iterate_rings.h>
#include <directional/Subdivision_new/shm_oneform_coefficients.h>
#include <directional/Subdivision/OneFormCoefficientProvider.h>

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

struct SubData_New
{
	using SMat = Eigen::SparseMatrix<double>;
	Eigen::MatrixXi E0, EF0, EI0, SFE0, F0;
	Eigen::MatrixXd V0, VK;
	Eigen::MatrixXi EK, EFK, EIK, SFEK, FK;
	SMat S_v, S_e, S_f, S_c;
	SubData_New(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int level): V0(V), F0(F)
	{
		directional::shm_edge_topology(F0, V0.rows(), E0, EF0, EI0, SFE0);
		directional::subdivision::build_shm_subdivision(F0, V0, E0, EF0, EI0, SFE0, level, FK,EK, EFK, EIK, SFEK, S_v, S_e, S_c, S_f);
		VK = S_v * V;
	}
};
void vertex_valences(const Eigen::MatrixXi& E, int vCount, Eigen::VectorXi& valences)
{
	valences.setConstant(vCount, 0);
	for(int i = 0; i < E.rows(); i++)
	{
		valences(E(i, 0)) += 1;
		valences(E(i, 1)) += 1;
	}
}
void neighbor_count(const Eigen::MatrixXi& EF, const Eigen::MatrixXi& SFE, Eigen::VectorXi& neighCount)
{
	neighCount.setConstant(SFE.rows(), 0);
	for (int i = 0; i < SFE.rows(); i++)
	{
		for(int j = 0; j < 3; j++)
		{
			if (EF(SFE(i, j), 1 - SFE(i, j + 3)) != -1) neighCount(i)++;
		}
	}
}
int valenceForFace(const Eigen::VectorXi& valences, const Eigen::MatrixXi& F, int face, int corner)
{
	return valences(F(face, corner));
}
void boundary_es(const Eigen::MatrixXi& EF, Eigen::VectorXi& BEs)
{
	BEs.setConstant(EF.rows(), 0);
	for (int i = 0; i < EF.rows(); i++)
	{
		if (EF(i, 0) == -1 || EF(i, 1) == -1)
		{
			BEs(i) = 1;
		}
	}
}
void boundary_verts(const Eigen::MatrixXi& E, const Eigen::MatrixXi& EF, int vCount, Eigen::VectorXi& BVs)
{
	BVs.setConstant(vCount, 0);
	for (int i = 0; i < E.rows(); i++)
	{
		if(EF(i,0)==-1 || EF(i,1) == -1)
		{
			BVs(E(i, 0)) = 1;
			BVs(E(i, 1)) = 1;
		}
	}
}
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

TEST_P(MeshTestFixture, QuadrisectionCorrectness_New)
{
	using namespace directional;
	using namespace directional::subdivision;
	const int vCount = m_mesh.V.rows();
	Eigen::MatrixXi E, EF, EI, SFE;
	directional::shm_edge_topology(m_mesh.F, m_mesh.V.rows(), E, EF, EI, SFE);

	const int fCount0 = m_mesh.F.rows();
	const int eCount0 = E.rows();

	ASSERT_TRUE(is_edgedata_consistent(m_mesh.F, E, EF, EI, SFE, true));
	Eigen::MatrixXi F1, E1, EF1, EI1, SFE1, E0ToEk;
	quadrisect(m_mesh.F, vCount, E, SFE, EF, EI, E0ToEk, F1, E1, SFE1, EF1, EI1);

	const int fCount1 = F1.rows();
	const int eCount1 = E1.rows();
	
	// Check new element counts
	EXPECT_EQ(fCount1, 4 * fCount0);
	EXPECT_EQ(eCount1, 2 * eCount0 + 3 * fCount0);
	// Require consistency for the quadrisected element
	ASSERT_TRUE(is_edgedata_consistent(F1, E1, EF1, EI1, SFE1, true));

	// Count boundary edges
	int bCount0 = 0;
	for (int e = 0; e < EF.rows(); e++)
	{
		if (EF(e, 0) == -1 || EF(e, 1) == -1) bCount0++;
	}
	int bCount1 = 0;
	for (int e = 0; e < EF1.rows(); e++)
	{
		if (EF1(e, 0) == -1 || EF1(e, 1) == -1) bCount1++;
	}
	// Check that boundary edges is correctly recorded and consistent with pre-quadrisection
	//EXPECT_EQ(bCount, ED1.boundaryEdgeCount);
	EXPECT_EQ(2 * bCount0, bCount1);

	// Check that sFE is correct.
	for (int f = 0; f < SFE.rows(); f++)
	{
		for (int j = 0; j < 3; j++)
		{
			EXPECT_EQ(m_mesh.F(f, j), F1(4 * f + j, j));

			const int e = SFE(f, j);
			const int c = SFE(f, j + 3);

			// Check the edge mapping
			EXPECT_EQ(E1(E0ToEk(e, 0), 0), E(e, 0));
			EXPECT_EQ(E1(E0ToEk(e, 1), 1), E(e, 1));
			EXPECT_EQ(E1(E0ToEk(e, 0), 1), vCount + e);
			EXPECT_EQ(E1(E0ToEk(e, 1), 0), vCount + e);
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

TEST_P(MeshTestFixture, DEC_SEQUENCE_New)
{
	using SMat = Eigen::SparseMatrix<double>;
	
	Eigen::MatrixXi F0, E0, EF0, EI0, SFE0;
	int vCount0 = m_mesh.V.rows();
	F0 = m_mesh.F;
	directional::shm_edge_topology(F0, vCount0, E0, EF0, EI0, SFE0);

	// Subdivided connectivity
	Eigen::MatrixXi F1, E1, EF1, EI1, SFE1;

	// Construct subdivision operators
	SMat S_V, S_F, S_E, S_C;
	using namespace directional::subdivision;
	build_shm_subdivision(F0, m_mesh.V, E0, EF0, EI0, SFE0, level, F1, E1, EF1, EI1, SFE1, S_V, S_E, S_C, S_F);
	
	// Test sizes
	ASSERT_EQ(S_E.rows(), E1.rows());
	ASSERT_EQ(S_C.rows(), E1.rows());
	ASSERT_EQ(S_E.cols(), E0.rows());
	ASSERT_EQ(S_C.cols(), E0.rows());

	Eigen::MatrixXd V1 = S_V * m_mesh.V;

	// Require new level edge data to be consistent
	//ASSERT_TRUE(EDL.isConsistent());

	// Test d1 d0 = 0 relation
	SparseMat d0_0, d0_K, d1_K;
	directional::DEC_d0(V1, F1, E1, SFE1, EF1, d0_K);
	directional::DEC_d1(V1, F1, E1, SFE1, EF1, d1_K);

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

TEST_P(MeshTestFixture, Compare_Coeffs_Old_New)
{
	SubData_New sd(m_mesh.V, m_mesh.F, level);
	auto handler = [this](const std::vector<int>& edges, const std::vector<int>& edgeOrients)
	{
		
		const bool isBoundary = edges.size() & 1;
		const int valence = (edges.size() + isBoundary) / 2;
		OneFormCoefficientProvider p{};
		std::string endMsg = std::string("even, ") + std::to_string(isBoundary);
		for(int eI = 0; eI < edges.size(); eI += 2)
		{
			std::vector<int> indsOld, indsNew;
			std::vector<double> coeffsOld, coeffsNew;
			if (isBoundary) p.getEvenBoundaryStencil(valence, eI, indsOld, coeffsOld);
			else p.getEvenRegularStencil(valence, eI, indsOld, coeffsOld);
			directional::subdivision::shm_oneform_coefficients(isBoundary, true, valence, eI, indsNew, coeffsNew);
			ASSERT_TRUE(helpers::isSame(indsOld, indsNew)) << "Inds not same at valence" << valence << ",loc" << eI <<endMsg << helpers::printCombos(indsOld, indsNew);
			ASSERT_TRUE(helpers::isSame(coeffsOld, coeffsNew)) << "Coeffs not same at valence" << valence << ",loc" << eI <<endMsg << helpers::printCombos(coeffsOld,coeffsNew);
		}
		endMsg = std::string("odd, ") + std::to_string(isBoundary);
		for (int eI = 1; eI < edges.size(); eI += 2)
		{
			std::vector<int> indsOld, indsNew;
			std::vector<double> coeffsOld, coeffsNew;
			if (isBoundary) p.getOddBoundaryStencil(valence, eI, indsOld, coeffsOld);
			else p.getOddRegularStencil(valence, eI, indsOld, coeffsOld);
			directional::subdivision::shm_oneform_coefficients(isBoundary, false, valence, eI, indsNew, coeffsNew);
			ASSERT_TRUE(helpers::isSame(indsOld, indsNew)) << "Inds not same at valence" << valence << ",loc" << eI << endMsg << helpers::printCombos(indsOld, indsNew);
			ASSERT_TRUE(helpers::isSame(coeffsOld, coeffsNew)) << "Coeffs not same at valence" << valence << ",loc" << eI << endMsg << helpers::printCombos(coeffsOld, coeffsNew);
		}

		
		ASSERT_TRUE(true);
	};
	directional::iterate_rings(m_mesh.V.rows(), sd.E0, sd.EF0, sd.EI0, sd.SFE0, handler);
}

TEST_P(MeshTestFixture, Compare_Old_New)
{
	SubData sd(m_mesh.V, m_mesh.F, level);
	SubData_New sdN(m_mesh.V, m_mesh.F, level);

	EXPECT_TRUE(helpers::isSame(sd.ED0.E, sdN.E0));
	EXPECT_TRUE(helpers::isSame(sd.ED0.EF, sdN.EF0));
	EXPECT_TRUE(helpers::isSame(sd.ED0.sFE, sdN.SFE0));

	EXPECT_TRUE(helpers::isSame(sd.EDL.E, sdN.EK));
	EXPECT_TRUE(helpers::isSame(sd.EDL.EF, sdN.EFK));
	EXPECT_TRUE(helpers::isSame(sd.EDL.sFE, sdN.SFEK));

	std::vector<std::tuple<int, int, double, double>> diffEls = helpers::getAbsDifference(sd.S_v, sdN.S_v, tolerance, 20);
	std::stringstream msg;
	if(diffEls.size() > 0)
	{
		for(auto el : diffEls)
		{
			Eigen::VectorXi valence;
			vertex_valences(sdN.E0, sdN.V0.rows(), valence);
			int r, c;
			double v1, v2;
			std::tie(r, c, v1, v2) = el;

			msg << "Failed at " << r << ',' << c << ", vals: " << v1 << " vs. " << v2 << ", deg: " << valence(c) << std::endl;
		}
		
	}
	ASSERT_EQ(diffEls.size(), 0) << "S_V not the same in new:" << std::endl  << msg.str();
	diffEls = helpers::getAbsDifference(sd.S_c, sdN.S_c, tolerance, 20);
	ASSERT_EQ(diffEls.size(), 0) << "S_C not the same in new";
	diffEls = helpers::getAbsDifference(sd.S_e, sdN.S_e, tolerance, 20);
	if (diffEls.size() > 0)
	{
		for (auto el : diffEls)
		{
			Eigen::VectorXi valence;
			Eigen::VectorXi BE,BV;
			vertex_valences(sdN.E0, sdN.V0.rows(), valence);
			boundary_es(sdN.EF0, BE);
			boundary_verts(sdN.E0, sdN.EF0,sdN.V0.rows(), BV);
			int r, c;
			double v1, v2;
			std::tie(r, c, v1, v2) = el;
			
			msg << "Failed at " << r << ',' << c << ", vals: " << v1 << " vs. " << v2 << ", isbnd:" << BE(c) << ",vbnd:["
			<< BV(sdN.E0(c,0)) << ',' << BV(sdN.E0(c, 1)) << ']'
			<< std::endl;
		}
	}
	ASSERT_EQ(diffEls.size(), 0) << "S_E not the same in new:" << std::endl << msg.str();
	diffEls = helpers::getAbsDifference(sd.S_f, sdN.S_f, tolerance, 20);
	if (diffEls.size() > 0)
	{
		for (auto el : diffEls)
		{
			Eigen::VectorXi valence;
			Eigen::VectorXi neighbCount;
			vertex_valences(sdN.E0, sdN.V0.rows(), valence);
			neighbor_count(sdN.EF0, sdN.SFE0, neighbCount);
			int r, c;
			double v1, v2;
			std::tie(r, c, v1, v2) = el;
			bool isEven = ((r + 1) % 4 == 0);

			msg << "Failed at " << r << ',' << c << ", vals: " << v1 << " vs. " << v2 << ", even: " << isEven << ",";
			if(isEven)
			{
				msg << "neighbors= " << neighbCount(c);
			}
			else
			{
				msg << "valence = " << valenceForFace(valence, sdN.F0, c, r % 4);
			}
			msg << std::endl;
		}

	}
	ASSERT_EQ(diffEls.size(), 0) << "S_F not the same in new:" << msg.str();
}

TEST_P(MeshTestFixture, SEC_D0_Commutation_New)
{
	SubData_New sd(m_mesh.V, m_mesh.F, level);

	// Test basic d0 commutation
	SparseMat d0_0, d0_K;
	directional::DEC_d0(m_mesh.V, m_mesh.F, sd.E0, sd.SFE0, sd.EF0, d0_0);
	directional::DEC_d0(sd.VK, sd.FK, sd.EK, sd.SFEK, sd.EFK, d0_K);
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
		boundary_verts(sd.E0, sd.EF0, sd.V0.rows(), BVs);
		boundary_verts(sd.EK, sd.EFK, sd.VK.rows(), BVsK);
		vertex_valences(sd.E0, sd.V0.rows(), valence);
		const int v0Count = sd.V0.rows();
		msg << std::endl;
		msg << "v0 count:" << v0Count << std::endl;
		for (int i = 0; i < diffEls.size(); i++)
		{
			const int r = std::get<0>(diffEls[i]), c = std::get<1>(diffEls[i]);;
			double v1 = std::get<2>(diffEls[i]), v2 = std::get<3>(diffEls[i]);

			const int v0K = sd.EK(r, 0), v1K = sd.EK(r, 1);

			msg << "El at " << r << ',' << c << "| src is bnd:" << BVs(c) << ", deg:"
				<< valence(c) << ", vals:" << v1 << ", " << v2 << ", E bnd:[" << BVsK(sd.EK(r, 0)) << ", " << BVsK(sd.EK(r, 1))
				<< "], EV:[" << sd.EK(r, 0) << ',' << sd.EK(r, 1) << "]"
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

/**
 * Second SEC commutation
 */
TEST_P(MeshTestFixture, SEC_D1_Commutation_New)
{
	using SMat = Eigen::SparseMatrix<double>;
	int level = 1;

	SubData_New sd(m_mesh.V, m_mesh.F, level);

	// Test basic d0 commutation
	SparseMat d1_0, d1_K;
	directional::DEC_d1(m_mesh.V, m_mesh.F, sd.E0, sd.SFE0, sd.EF0, d1_0); // coarse level
	directional::DEC_d1(sd.VK, sd.FK, sd.EK, sd.SFEK, sd.EFK, d1_K); // Fine level
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

TEST_P(MeshTestFixture, SHM_Curl_Commutation_New)
{
	using SMat = Eigen::SparseMatrix<double>;
	int level = 1;

	SubData_New sd(m_mesh.V, m_mesh.F, level);
	SMat A_EToF0, A_EToFK;
	directional::Edge_To_Face_Average(sd.SFE0, sd.E0.rows(), A_EToF0);
	directional::Edge_To_Face_Average(sd.SFEK, sd.EK.rows(), A_EToFK);

	// To compare
	SparseMat lhs = sd.S_f * A_EToF0;
	SparseMat rhs = A_EToFK * sd.S_c;

	std::vector<std::tuple<int, int, double, double>> diffEls = helpers::getAbsDifference(lhs, rhs, tolerance, 100);
	std::stringstream msg;
	/*
	msg << std::endl;
	Eigen::VectorXi BVs;
	Eigen::VectorXi valence;
	sd.ED0.boundaryVerts(BVs);
	sd.ED0.vertexValence(valence);
	for (int i = 0; i < diffEls.size(); i++)
	{
		const int r = std::get<0>(diffEls[i]), c = std::get<1>(diffEls[i]);;
		double v1 = std::get<2>(diffEls[i]), v2 = std::get<3>(diffEls[i]);
		const int f = r / 4;
		const int corner = r - 4 * f;
		bool inOrigF = false;
		int ind = -1;
		for (int i = 0; i < 3; i++)
		{
			if (sd.ED0.sFE(f, i) == c)
			{
				inOrigF = true;
				ind = i;
				break;
			}
		}
		msg << "At [" << r << ',' << c << "], vals [" << v1 << "," << v2 << "], original corner:" << corner << ", in orig f: " << inOrigF << "," << ind << " ";
		if (corner < 3)
		{
			msg << "Vert val:" << valence(sd.ED0.F(f, corner)) << ", is bnd: " << BVs(sd.ED0.F(f, corner));
		}
		msg << std::endl;
	}*/
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

TEST_P(MeshTestFixture, SHM_Gradient_Commutation_New)
{
	using SMat = Eigen::SparseMatrix<double>;
	int level = 1;

	SubData_New sd(m_mesh.V, m_mesh.F, level);
	// Coarse level operators
	helpers::Gamma2_Ops_New ops0(m_mesh, sd.E0, sd.EF0, sd.EI0,sd.SFE0);
	directional_fixtures::Mesh meshK;
	meshK.F = sd.FK;
	meshK.V = sd.VK;
	// Fine level operators
	helpers::Gamma2_Ops_New opsK(meshK, sd.EK, sd.EFK, sd.EIK, sd.SFEK);
	
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

TEST_P(MeshTestFixture, SHM_SGamma_SE_Commutation_New)
{
	using SMat = Eigen::SparseMatrix<double>;

	SubData_New sd(m_mesh.V, m_mesh.F, level);
	// Coarse level operators
	helpers::Gamma2_Ops_New ops0(m_mesh, sd.E0, sd.EF0, sd.EI0, sd.SFE0);
	// Fine level operators
	directional_fixtures::Mesh meshK;
	meshK.F = sd.FK;
	meshK.V = sd.VK;
	helpers::Gamma2_Ops_New opsK(meshK, sd.EK, sd.EFK, sd.EIK, sd.SFEK);

	SMat S_Gamma;
	directional::block_diag({ &sd.S_e, &sd.S_c }, S_Gamma);
	S_Gamma = opsK.Decomp_To_G2 * S_Gamma * ops0.G2_To_Decomp;


	// To compare
	SparseMat lhs = sd.S_e * ops0.Gamma2_To_Oneform;
	SparseMat rhs = opsK.Gamma2_To_Oneform * S_Gamma;

	std::vector<std::tuple<int, int, double, double>> diffEls = helpers::getAbsDifference(lhs, rhs, tolerance, 100);
	EXPECT_TRUE(diffEls.size() == 0) << "Avg_Gamma_To_E does not commute via S_e/S_Gamma:" << helpers::tripletsToString(diffEls);
}

INSTANTIATE_TEST_SUITE_P(Commutations, MeshTestFixture,
	testing::ValuesIn(SIMPLE_FILES));

int main(int argc, char **argv) {
	
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}