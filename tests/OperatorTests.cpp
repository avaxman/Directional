#define BOOST_TEST_MODULE Operators
#ifdef WIN32
#define NOMINMAX
#endif
#include <boost/test/included/unit_test.hpp>
#include <boost/test/data/test_case.hpp>
#include <boost/test/data/monomorphic.hpp>
#define DIR_ASSERT(X) BOOST_REQUIRE(X)
#include <Eigen/Eigen>
#include <igl/read_triangle_mesh.h>
#include <directional/Subdivision/subdivision.h>
#include <directional/Subdivision/DCEL.h>
#include <igl/max.h>
#include <igl/cat.h>
#include "Helpers.h"
#define SIMPLE_FILES {"bimba.off", "horser.off", "chipped-torus.obj"}

#include <directional/FEM_suite.h>
#include <directional/Gamma_suite.h>

double tolerance = 1e-7;

struct FEM_Ops
{
	using SparseMat = Eigen::SparseMatrix<double>;
	SparseMat Gv, Ge, J, C, D;

	FEM_Ops(const helpers::Mesh& m)
	{
		EdgeData ED;
		directional::construct_edge_topology(m.F, ED.E, ED.EF, ED.sFE, ED.EI);
		directional::FEM_suite(m.V, m.F, ED.E, ED.sFE, ED.EF, Gv, Ge, J, C, D);
	}
	FEM_Ops() {}
	void constructOperators(const EdgeData& ED, const helpers::Mesh& m)
	{
		directional::FEM_suite(m.V, m.F, ED.E, ED.sFE, ED.EF, Gv, Ge, J, C, D);
	}
};
struct DCELTester
{
	const EdgeData& ED;
	Eigen::VectorXi Vvisited;
	Eigen::MatrixXi Evisited;
	Eigen::VectorXi BoundaryVs;
	Eigen::VectorXi valences;
	DCELTester(const EdgeData& ED): ED(ED)
	{
		Vvisited.setConstant(ED.vertexCount(), 0);
		Evisited.setConstant(ED.edgeCount(), 2, 0);
		ED.boundaryVerts(BoundaryVs);
		for(int i = 0; i < BoundaryVs.rows(); i++)
		{
			BOOST_REQUIRE(BoundaryVs(i) == 0 || BoundaryVs(i) == 1);
		}
		ED.vertexValence(valences);
	}
	int startVert(int ind, const std::vector<int>& edges, const std::vector<int>& edgeSides)
	{
		return ED.E(edges[ind], edgeSides[ind]);
	}
	int endVert(int ind, const std::vector<int>& edges, const std::vector<int>& edgeSides)
	{
		return ED.E(edges[ind], 1 - edgeSides[ind]);
	}
	int lface(int ind, const std::vector<int>& edges, const std::vector<int>& edgeSides)
	{
		return ED.EF(edges[ind], edgeSides[ind]);
	}
	int rface(int ind, const std::vector<int>& edges, const std::vector<int>& edgeSides)
	{
		return ED.EF(edges[ind], 1-edgeSides[ind]);
	}
	void handleRegularRing(const std::vector<int>& edges, const std::vector<int>& edgeSides)
	{
		BOOST_TEST(edges.size() == edgeSides.size(), "Sizes do not match");
		const int vCentral = ED.E(edges[0], edgeSides[0]);
		BOOST_REQUIRE_EQUAL(Vvisited(vCentral), 0);
		BOOST_TEST_CHECKPOINT(std::string("Central vert:") + std::to_string(vCentral));
		for(int i = 0; i < edges.size()/2; i++)
		{
			const int spoke = 2 * i;
			const int ring = 2 * i + 1;
			// Test connection to central vertex
			int v = startVert(spoke, edges, edgeSides);
			BOOST_REQUIRE_MESSAGE(v == vCentral, "Invalid spoke start vert, expected " << vCentral << " got " << v);
			// Test connection between spoke and ring
			v = endVert(spoke, edges, edgeSides);
			int v2 = startVert(ring, edges, edgeSides);
			BOOST_REQUIRE_MESSAGE(v == v2, "Spoke and ring do not meet in vertex: " << v << " vs " << v2);
			// Test common face
			int f1 = lface(spoke, edges, edgeSides), f2 = lface(ring, edges, edgeSides);
			BOOST_REQUIRE_MESSAGE(f1 == f2, "Spoke and ring do not share face at right side: " << f1 << " vs " << f2);
			// Check connection to next pair
			const int nextSpoke = (2 * i + 2) % edges.size();
			f2 = rface(nextSpoke, edges, edgeSides);
			BOOST_REQUIRE_MESSAGE(f1 == f2 ,"Next spoke is not properly related via face: " << f1 << " vs " << f2);
			Evisited(edges[spoke], edgeSides[spoke])+=1;
			Evisited(edges[ring], edgeSides[ring]) += 1;

			// Check face and vertices
			const int f = ED.EF(edges[spoke], edgeSides[spoke]);
			const int vs[3] = { ED.E(edges[spoke], edgeSides[spoke]), ED.E(edges[ring], edgeSides[ring]), ED.E(edges[ring], 1- edgeSides[ring]) };
			//BOOST_TEST_MESSAGE("Verts:" << vs[0] << ',' << vs[1] << ',' << vs[2] << ", row: " << ED.F.row(f));
			int start = -1;
			for(int i = 0; i < 3; i++)
			{
				if(vs[0] == ED.F(f,i))
				{

					start = i;
					break;
				}
			}
			//BOOST_TEST_MESSAGE("Start:" << start);
			for(int i = 0; i < 3; i++)
			{
				BOOST_REQUIRE_EQUAL(ED.F(f, (start + i) % 3), vs[i]);
			}

		}
		Vvisited(vCentral)++;
	}
	void handleBoundaryRing(const std::vector<int>& edges, const std::vector<int>& edgeSides)
	{
		//TODO
		if(ED.boundaryEdgeCount == 0)
			BOOST_REQUIRE_MESSAGE(false, "Not expecting boundary");
		else
		{
			// First element and last element are boundary
			BOOST_REQUIRE_EQUAL(ED.EF(edges[0], 1 - edgeSides[0]), -1);
			BOOST_REQUIRE_EQUAL(ED.EF(edges.back(), edgeSides.back()), -1);

			BOOST_TEST(edges.size() == edgeSides.size(), "Sizes do not match");
			const int vCentral = ED.E(edges[0], edgeSides[0]);
			BOOST_REQUIRE_EQUAL(Vvisited(vCentral), 0);
			BOOST_TEST_CHECKPOINT(std::string("Central vert boundary:") + std::to_string(vCentral));
			BOOST_REQUIRE_EQUAL(valences(vCentral), (edges.size() + 1)/2);
			// Update visisted edge at end of ring
			Evisited(edges.back(), edgeSides.back()) += 1;

			BoundaryVs(vCentral) -= 1;

			for (int i = 0; i < (edges.size() -1) / 2; i++)
			{
				const int spoke = 2 * i;
				const int ring = 2 * i + 1;
				
				// Test connection to central vertex
				int v = startVert(spoke, edges, edgeSides);
				BOOST_REQUIRE_MESSAGE(v == vCentral, "Invalid spoke start vert, expected " << vCentral << " got " << v);
				
				// Test connection between spoke and ring
				v = endVert(spoke, edges, edgeSides);
				int v2 = startVert(ring, edges, edgeSides);
				BOOST_REQUIRE_MESSAGE(v == v2, "Spoke and ring do not meet in vertex: " << v << " vs " << v2);

				// Test common face
				int f1 = lface(spoke, edges, edgeSides), f2 = lface(ring, edges, edgeSides);
				BOOST_REQUIRE_MESSAGE(f1 == f2, "Spoke and ring do not share face at right side: " << f1 << " vs " << f2);
				
				// Check connection to next pair
				const int nextSpoke = (2 * i + 2) % edges.size();
				f2 = rface(nextSpoke, edges, edgeSides);
				BOOST_REQUIRE_MESSAGE(f1 == f2, "Next spoke is not properly related via face: " << f1 << " vs " << f2);
				Evisited(edges[spoke], edgeSides[spoke]) += 1;
				Evisited(edges[ring], edgeSides[ring]) += 1;

				// Check face and vertices
				const int f = ED.EF(edges[spoke], edgeSides[spoke]);
				const int vs[3] = { ED.E(edges[spoke], edgeSides[spoke]), ED.E(edges[ring], edgeSides[ring]), ED.E(edges[ring], 1 - edgeSides[ring]) };
				//BOOST_TEST_MESSAGE("Verts:" << vs[0] << ',' << vs[1] << ',' << vs[2] << ", row: " << ED.F.row(f));
				int start = -1;
				for (int i = 0; i < 3; i++)
				{
					if (vs[0] == ED.F(f, i))
					{

						start = i;
						break;
					}
				}
				//BOOST_TEST_MESSAGE("Start:" << start);
				for (int i = 0; i < 3; i++)
				{
					BOOST_REQUIRE_EQUAL(ED.F(f, (start + i) % 3), vs[i]);
				}

			}
			Vvisited(vCentral)++;
		}
	}
	void checkVisit()
	{
		// Check all vertices were visited
		for(int i = 0; i < Vvisited.rows(); i++)
		{
			BOOST_TEST(Vvisited(i) == 1, std::string("Missed vertex ") + std::to_string(i));
		}
		for(int i = 0; i < Evisited.rows(); i++)
		{
			if(ED.EF(i,0) == -1)
			{
				BOOST_TEST(Evisited(i, 0) == 1, "Expected to see edge " << i << " once at orient 0, but saw " << Evisited(i, 0));
			}
			else
			{
				BOOST_TEST(Evisited(i, 0) == 2, "Expected to see edge " << i << " twice at orient 0, but saw " << Evisited(i, 0));
			}
			if (ED.EF(i, 1) == -1)
			{
				BOOST_TEST(Evisited(i, 1) == 1, "Expected to see edge " << i << " once at orient 1, but saw " << Evisited(i, 1));
			}
			else
			{
				BOOST_TEST(Evisited(i, 1) == 2, "Expected to see edge " << i << " twice at orient 1, but saw " << Evisited(i, 1));
			}
			
		}
		for(int i = 0; i < BoundaryVs.rows(); i++)
		{
			BOOST_TEST(BoundaryVs(i) == 0, "Invalid boundary : at " << i << ": " << BoundaryVs(i));
		}
	}
};
using SparseMat = Eigen::SparseMatrix<double>;

std::string description(const Eigen::SparseMatrix<double>& mat)
{
	std::stringstream str;
	str << "Cols:" << mat.cols() << ", rows: " << mat.rows();
	return str.str();
}


BOOST_DATA_TEST_CASE(DCELTest, helpers::MeshDataset(TEST_SHARED_PATH, SIMPLE_FILES))
{
	EdgeData ED(sample.F);
	DCEL dcel(ED);
	DCELTester test(ED);
	dcel.iterateRings(test);
	test.checkVisit();

	EdgeData EDL;
	Eigen::MatrixXi E0ToEk;
	ED.quadrisect(sample.V.rows(), E0ToEk, EDL);
	DCEL dcel2(EDL);
	DCELTester test2(EDL);
	dcel2.iterateRings(test2);
	test2.checkVisit();

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
		for(int j = 0; j < 4; j++)
		{
			if (E0ToEk(e, j) != -1) {
				const int eIn = E0ToEk(e, j);
				int bCountInner = vKIsBoundary(EDL.E(eIn, 0)) + vKIsBoundary(EDL.E(eIn, 1));
				const int fullBound = lBound == -1 || rBound == -1;
				switch(j)
				{
				case 0:
					BOOST_REQUIRE_MESSAGE(bCountInner == sBound + fullBound, "Invalid boundary at start edge: " << bCountInner << " vs " << sBound + fullBound);
					break;
				case 1:
					BOOST_REQUIRE_MESSAGE(bCountInner == eBound + fullBound, "Invalid boundary at end edge: " << bCountInner << " vs " << eBound + fullBound);
					break;
				case 2:
					BOOST_REQUIRE_MESSAGE(bCountInner == (sBound & lBound) + (eBound & lBound), "Invalid boundary at left edge: " << bCountInner << " vs " << (sBound & lBound) + (eBound & lBound));
					break;
				case 3:
					BOOST_REQUIRE_MESSAGE(bCountInner == (sBound & rBound) + (eBound & rBound), "Invalid boundary at right edge: " << bCountInner << " vs " << (sBound & rBound) + (eBound & rBound));
					break;
				default:
					break;
				}
			}
		}
	}
}

BOOST_DATA_TEST_CASE(DECTests, helpers::MeshDataset(TEST_SHARED_PATH, SIMPLE_FILES))
{
	using SMat = Eigen::SparseMatrix<double>;
	EdgeData ED(sample.F);
	SMat d0, d1;
	directional::DEC_d0(sample.V, sample.F, ED.E, ED.sFE, ED.EF, d0);
	directional::DEC_d1(sample.V, sample.F, ED.E, ED.sFE, ED.EF, d1);
	SMat prod = d1 * d0;
	auto els = helpers::getDifferenceFromZero(prod, 1e-7, 100);
	BOOST_REQUIRE(els.size() == 0, std::string("d1 d0 is not zero") + helpers::tripletsToString(els));
}

BOOST_DATA_TEST_CASE(GammaInverseTest, helpers::MeshDataset(TEST_SHARED_PATH, SIMPLE_FILES))
{
	using SMat = Eigen::SparseMatrix<double>;

	// Construct gamma operators
	helpers::Gamma2_Ops ops_Gamma(sample);

	SMat possiblyId = ops_Gamma.Chi_To_Gamma2 * ops_Gamma.Gamma2_To_Chi;
	auto diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	BOOST_TEST(diffs.size() == 0, std::string("P P^Inv is not identity for Gamma2:") + helpers::tripletsToString(diffs));
	// Note that P^Inv P is not identity by definition but a normal component remover.
}

BOOST_DATA_TEST_CASE(G2G3ConversionTest, helpers::MeshDataset(TEST_SHARED_PATH, SIMPLE_FILES))
{
	using SMat = Eigen::SparseMatrix<double>;

	// Construct gamma operators
	helpers::Gamma2_Ops ops_Gamma(sample);
	SMat G2_To_G3, G3_To_G2;
	EdgeData ED(sample.F);
	directional::Gamma2_To_Gamma3(ED.sFE, G2_To_G3);
	directional::Gamma3_To_Gamma2(ED.sFE, G3_To_G2);

	/*SMat possiblyId = G2_To_G3 * G3_To_G2;
	std::vector<std::tuple<int, int, double>> diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	BOOST_TEST(diffs.size() == 0, std::string("G3->G2->G3 is not identity for Gamma2:") + helpers::tripletsToString(diffs));*/

	SMat possiblyId = G3_To_G2 * G2_To_G3;
	std::vector < std::tuple<int, int, double>> diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	BOOST_TEST(diffs.size() == 0, std::string("G2->G3->G2 is not identity for Decomp:") + helpers::tripletsToString(diffs));
}


BOOST_DATA_TEST_CASE(DecompTest, helpers::MeshDataset(TEST_SHARED_PATH, SIMPLE_FILES))
{
	using SMat = Eigen::SparseMatrix<double>;

	// Construct gamma operators
	helpers::Gamma2_Ops ops_Gamma(sample);

	EdgeData ED(sample.F);
	SMat G3To1Form, C, DCToG3, G3ToDC;
	directional::Gamma3_To_Decomp(ED.EF, ED.EI, ED.sFE.rows(), G3ToDC);
	directional::Decomp_To_Gamma3(ED.EF, ED.EI, ED.sFE.rows(), DCToG3);

	SMat possiblyId = ops_Gamma.G3ToG2 * ops_Gamma.G2ToG3;
	std::vector < std::tuple<int, int, double>> diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	BOOST_TEST(diffs.size() == 0, std::string("G2->G3->G2 is not identity for Decomp:") + helpers::tripletsToString(diffs));

	possiblyId = ops_Gamma.G3ToG2 * DCToG3 * G3ToDC * ops_Gamma.G2ToG3;
	diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	BOOST_TEST(diffs.size() == 0, std::string("G2->Decomp->G2 is not identity for Gamma2:") + helpers::tripletsToString(diffs));

	possiblyId = ops_Gamma.Decomp_To_G2 * ops_Gamma.G2_To_Decomp;
	diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	BOOST_TEST(diffs.size() == 0, std::string("G2->Decomp->G2 from ops struct is not identity for Gamma2:") + helpers::tripletsToString(diffs));

	// Note that this is identity up to some non-null sum removing factors in hte matrix.
	/*possiblyId = G3ToDC * ops_Gamma.G2ToG3 * ops_Gamma.G3ToG2*DCToG3 ;
	diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	BOOST_TEST(diffs.size() == 0, std::string("Decomp->G2->Decomp is not identity for Decomp:") + helpers::tripletsToString(diffs));*/
}

BOOST_DATA_TEST_CASE(Decomp3Test, helpers::MeshDataset(TEST_SHARED_PATH, SIMPLE_FILES))
{
	using SMat = Eigen::SparseMatrix<double>;

	EdgeData ED(sample.F);

	// Construct gamma operators
	SMat G3To1Form, C, DCToG3, G3ToDC;
	directional::Gamma3_To_Decomp(ED.EF, ED.EI, ED.sFE.rows(), G3ToDC);
	//igl::cat(2, SMat(2.0 * G3To1Form.transpose()), SMat(C.transpose()), DCToG3);
	directional::Decomp_To_Gamma3(ED.EF, ED.EI, ED.sFE.rows(), DCToG3);


	SMat possiblyId = DCToG3 * G3ToDC;
	std::vector<std::tuple<int, int, double>> diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	BOOST_TEST(diffs.size() == 0, std::string("G3->Decomp->G3 is not identity for Gamma2:") + helpers::tripletsToString(diffs));

	possiblyId = G3ToDC * DCToG3;
	diffs = helpers::getDifferenceFromIdentity(possiblyId, tolerance);
	BOOST_TEST(diffs.size() == 0, std::string("Decomp->G3->Decomp is not identity for Decomp:") + helpers::tripletsToString(diffs));
}

BOOST_DATA_TEST_CASE(OperatorTests, helpers::MeshDataset(TEST_SHARED_PATH, SIMPLE_FILES))
{
	using SMat = Eigen::SparseMatrix<double>;

	// Construct gamma operators
	helpers::Gamma2_Ops ops_Gamma(sample);
	FEM_Ops ops_FEM(sample);

	BOOST_REQUIRE(ops_Gamma.Chi_To_Gamma2.cols() == ops_FEM.D.cols(), "Operators Chi to Gamma2 and Divergence of chi do not have same column number");
	BOOST_REQUIRE(ops_Gamma.D.rows() == ops_FEM.D.rows(), "Operators D_Gamma and D_Chi to not map to same space");

	SMat diff = (ops_Gamma.D * ops_Gamma.Chi_To_Gamma2 - ops_FEM.D);
	double val = helpers::maxAbs(diff);
	BOOST_TEST(val < tolerance, std::string("D_Gamma * P differs from D_Chi: ") +std::to_string(val));
	SMat diff2 = (ops_Gamma.D  - ops_FEM.D * ops_Gamma.Gamma2_To_Chi);
	val = helpers::maxAbs(diff2);
	BOOST_TEST(val < tolerance, std::string("D_Gamma differs from D_Chi * PInv: ") + std::to_string(val));
}
