#define BOOST_TEST_MODULE Commutations
#define NOMINMAX
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
//#define DIR_ASSERT(x) BOOST_REQUIRE(x)
#define DIR_ASSERT(x) BOOST_REQUIRE_MESSAGE((x), "Failed condition " #x "at" __FILE__);
#define DIR_ASSERT_M(x, m) BOOST_REQUIRE_MESSAGE(x, m)
//#define assert(x) BOOST_TEST(x)
#include <Eigen/Eigen>
#include <igl/read_triangle_mesh.h>
#include <directional/Subdivision/subdivision.h>
#include <igl/max.h>
#include <igl/cat.h>
#include "Helpers.h"
#include <directional/FEM_suite.h>
#include <directional/Gamma_suite.h>
#include <directional/block_diag.h>
#define SIMPLE_FILES {"bimba.off", "horser.off","chipped-torus.obj"}


// Type aliases
using SparseMat = Eigen::SparseMatrix<double>;

/**
 *		TEST SETTINGS
 */
double tolerance = 1e-7;
int level = 1;

std::string description(const Eigen::SparseMatrix<double>& mat)
{
	std::stringstream str;
	str << "Cols:" << mat.cols() << ", rows: " << mat.rows();
	return str.str();
}

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

// Test correctness of the quadrisection step.
BOOST_DATA_TEST_CASE(EdgeQuadrisect, helpers::MeshDataset(TEST_SHARED_PATH, SIMPLE_FILES))
{
	EdgeData ED0,ED1;
	ED0.construct(sample.F);

	BOOST_TEST(ED0.isConsistent(true));
	Eigen::MatrixXi E0ToEk;
	ED0.quadrisect(sample.V.rows(), E0ToEk, ED1);
	//BOOST_REQUIRE_EQUAL(ED0.boundaryEdgeCount, 0);
	BOOST_REQUIRE_EQUAL(ED1.faceCount(), 4 * ED0.faceCount());
	BOOST_REQUIRE_EQUAL(ED1.edgeCount(), 2 * ED0.edgeCount() +  3 * ED0.faceCount());
	BOOST_TEST(ED1.isConsistent(true));
	int bCount = 0;
	for(int e = 0; e < ED1.EF.rows(); e++)
	{
		if (ED1.EF(e, 0) == -1 || ED1.EF(e, 1) == -1) bCount++;
	}
	BOOST_TEST_MESSAGE(std::string("BCOUNt:") + std::to_string(bCount));
	BOOST_REQUIRE_EQUAL(bCount, ED1.boundaryEdgeCount);
	BOOST_REQUIRE_EQUAL(2 * ED0.boundaryEdgeCount, ED1.boundaryEdgeCount);

	for(int f = 0; f < ED0.sFE.rows(); f++)
	{
		for(int j = 0; j < 3; j++)
		{
			BOOST_REQUIRE_EQUAL(ED0.F(f, j), ED1.F(4*f+j,j));

			const int e = ED0.sFE(f, j);
			const int c = ED0.sFE(f, j + 3);
			
			BOOST_REQUIRE_EQUAL(ED1.E(E0ToEk(e, 0), 0), ED0.E(e, 0));
			BOOST_REQUIRE_EQUAL(ED1.E(E0ToEk(e, 1), 1), ED0.E(e, 1));
			BOOST_REQUIRE_EQUAL(ED1.E(E0ToEk(e, 0), 1), ED0.vertexCount() + e);
			BOOST_REQUIRE_EQUAL(ED1.E(E0ToEk(e, 1), 0), ED0.vertexCount() + e);
		}
	}
}

BOOST_DATA_TEST_CASE(FineOperators, helpers::MeshDataset(TEST_SHARED_PATH, SIMPLE_FILES))
{
	using SMat = Eigen::SparseMatrix<double>;
	EdgeData ED0, EDL;
	SMat S_V, S_F, S_E, S_C;
	subdivision_operators(sample.F, sample.V, level, S_V, S_F, S_E, S_C, ED0, EDL);

	BOOST_REQUIRE_EQUAL(S_E.rows(), EDL.edgeCount());
	BOOST_REQUIRE_EQUAL(S_C.rows(), EDL.edgeCount());
	BOOST_REQUIRE_EQUAL(S_E.cols(), ED0.edgeCount());
	BOOST_REQUIRE_EQUAL(S_C.cols(), ED0.edgeCount());

	helpers::Mesh meshK;
	meshK.F = EDL.F;
	// Subdivided upwards.
	meshK.V = S_V * sample.V;

	BOOST_TEST(EDL.isConsistent());

	// Test d1 d0 = 0 relation
	SparseMat d0_0, d0_K, d1_K;
	directional::DEC_d0(meshK.V, meshK.F, EDL.E, EDL.sFE, EDL.EF, d0_K);
	directional::DEC_d1(meshK.V, meshK.F, EDL.E, EDL.sFE, EDL.EF, d1_K);

	auto diffEls2 = helpers::getDifferenceFromZero(SparseMat(d1_K* d0_K), tolerance, 100);
	BOOST_TEST(diffEls2.size() == 0, std::string("d1_K d0_K is not zero") + helpers::tripletsToString(diffEls2));
}

/**
 * First SEC commutation
 */
BOOST_DATA_TEST_CASE(S_e_S_v_d0_commutation, helpers::MeshDataset(TEST_SHARED_PATH, SIMPLE_FILES))
{
	SubData sd(sample.V, sample.F, level);

	// Test basic d0 commutation
	SparseMat d0_0, d0_K;
	directional::DEC_d0(sample.V, sample.F, sd.ED0.E, sd.ED0.sFE, sd.ED0.EF, d0_0);
	directional::DEC_d0(sd.VK, sd.EDL.F, sd.EDL.E, sd.EDL.sFE, sd.EDL.EF, d0_K);
	SparseMat d0_p1 = sd.S_e * d0_0;
	SparseMat d0_p2 = d0_K * sd.S_v;

	std::vector<std::tuple<int, int, double, double>> diffEls = helpers::getAbsDifference(d0_p1, d0_p2, tolerance, 100);
	std::stringstream msg;
	if(diffEls.size() > 0)
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
		for(int i = 0; i < diffEls.size(); i++)
		{
			const int r = std::get<0>(diffEls[i]), c = std::get<1>(diffEls[i]);;
			double v1 = std::get<2>(diffEls[i]), v2 = std::get<3>(diffEls[i]);

			const int v0K = sd.EDL.E(r, 0), v1K = sd.EDL.E(r, 1);

			msg << "El at " << r << ',' << c << "| src is bnd:" << BVs(c) << ", deg:"
				<< valence(c) << ", vals:" << v1 << ", " << v2 << ", E bnd:[" << BVsK(sd.EDL.E(r,0)) << ", " << BVsK(sd.EDL.E(r, 1))
			<< "], EV:[" << sd.EDL.E(r,0) << ',' << sd.EDL.E(r,1) <<"]"
			<< std::endl;
		}
	}
	BOOST_TEST(diffEls.size() == 0, std::string("d0 does not commute via S_e/S_v") + msg.str());
}

/**
 * Second SEC commutation
 */
BOOST_DATA_TEST_CASE(S_e_S_f_d1_commutation, helpers::MeshDataset(TEST_SHARED_PATH, SIMPLE_FILES))
{
	using SMat = Eigen::SparseMatrix<double>;
	int level = 1;

	SubData sd(sample.V, sample.F, level);

	// Test basic d0 commutation
	SparseMat d1_0, d1_K;
	directional::DEC_d1(sample.V, sample.F, sd.ED0.E, sd.ED0.sFE, sd.ED0.EF, d1_0);
	// Construct fine level d1s
	directional::DEC_d1(sd.VK, sd.EDL.F, sd.EDL.E, sd.EDL.sFE, sd.EDL.EF, d1_K);
	// To compare
	SparseMat d1_p1 = sd.S_f * d1_0;
	SparseMat d1_p2 = d1_K * sd.S_e;

	std::vector<std::tuple<int, int, double, double>> diffEls = helpers::getAbsDifference(d1_p1, d1_p2, tolerance, 100);
	BOOST_TEST(diffEls.size() == 0, std::string("d1 does not commute via S_e/S_f") + helpers::tripletsToString(diffEls));
}

BOOST_DATA_TEST_CASE(S_c_S_f_AEToF_commutation, helpers::MeshDataset(TEST_SHARED_PATH, SIMPLE_FILES))
{
	using SMat = Eigen::SparseMatrix<double>;
	int level = 1;

	SubData sd(sample.V, sample.F, level);
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
	BOOST_TEST(diffEls.size() == 0, std::string("Average_E_To_F does not commute via S_c/S_f") + msg.str());
}

BOOST_DATA_TEST_CASE(S_GammaGrad_Comm, helpers::MeshDataset(TEST_SHARED_PATH, SIMPLE_FILES))
{
	using SMat = Eigen::SparseMatrix<double>;
	int level = 1;

	SubData sd(sample.V, sample.F, level);
	// Coarse level operators
	helpers::Gamma2_Ops ops0(sd.ED0, sample.V);
	// Fine level operators
	helpers::Gamma2_Ops opsK(sd.EDL, sd.VK);

	SMat S_Gamma;
	directional::block_diag({ &sd.S_e, &sd.S_c }, S_Gamma);
	S_Gamma = opsK.Decomp_To_G2 * S_Gamma * ops0.G2_To_Decomp;


	// To compare
	SparseMat lhs = S_Gamma * ops0.Gv;
	SparseMat rhs = opsK.Gv * sd.S_v;

	std::vector<std::tuple<int, int, double, double>> diffEls = helpers::getAbsDifference(lhs, rhs, tolerance, 100);
	BOOST_TEST(diffEls.size() == 0, std::string("d1 does not commute via S_e/S_f") + helpers::tripletsToString(diffEls));
}

BOOST_DATA_TEST_CASE(S_Gamma_S_E_A_Comm, helpers::MeshDataset(TEST_SHARED_PATH, SIMPLE_FILES))
{
	using SMat = Eigen::SparseMatrix<double>;

	SubData sd(sample.V, sample.F, level);
	// Coarse level operators
	helpers::Gamma2_Ops ops0(sd.ED0, sample.V);
	// Fine level operators
	helpers::Gamma2_Ops opsK(sd.EDL, sd.VK);

	SMat S_Gamma;
	directional::block_diag({ &sd.S_e, &sd.S_c }, S_Gamma);
	S_Gamma = opsK.Decomp_To_G2 * S_Gamma * ops0.G2_To_Decomp;


	// To compare
	SparseMat lhs = sd.S_e * ops0.Gamma2_To_Oneform;
	SparseMat rhs = opsK.Gamma2_To_Oneform * S_Gamma;

	std::vector<std::tuple<int, int, double, double>> diffEls = helpers::getAbsDifference(lhs, rhs, tolerance, 100);
	BOOST_TEST(diffEls.size() == 0, std::string("Avg_Gamma_To_E does not commute via S_e/S_Gamma") + helpers::tripletsToString(diffEls));
}
