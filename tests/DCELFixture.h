#ifndef DCEL_FIXTURE_H
#define DCEL_FIXTURE_H
#include <directional/Subdivision/DCEL.h>
#include <gtest/gtest.h>
#include "Helpers.h"

namespace directional_fixtures
{

	class DCELTestFixture: public testing::TestWithParam<const char*>
	{
	protected:
		Mesh m_mesh;
		void SetUp() override
		{
			// Resolve against shared path
			m_mesh.fileName = std::string(TEST_SHARED_PATH) + "/";
			m_mesh.fileName = m_mesh.fileName + GetParam();
			// Read mesh
			bool success = igl::read_triangle_mesh(m_mesh.fileName, m_mesh.V, m_mesh.F);
			ASSERT_TRUE(success);
			ASSERT_GT(m_mesh.F.rows(), 0) << "No faces in mesh";
			ED = m_mesh.getED();

			// Setup members
			Vvisited.setConstant(ED.vertexCount(), 0);
			Evisited.setConstant(ED.edgeCount(), 2, 0);
			ED.boundaryVerts(BoundaryVs);
			for (int i = 0; i < BoundaryVs.rows(); i++)
			{
				ASSERT_TRUE(BoundaryVs(i) == 0 || BoundaryVs(i) == 1);
			}
			ED.vertexValence(valences);
		}
		EdgeData ED;
		Eigen::VectorXi Vvisited;
		Eigen::MatrixXi Evisited;
		Eigen::VectorXi BoundaryVs;
		Eigen::VectorXi valences;

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
			return ED.EF(edges[ind], 1 - edgeSides[ind]);
		}
	public:
		void handleRegularRing(const std::vector<int>& edges, const std::vector<int>& edgeSides)
		{
			ASSERT_TRUE(edges.size() == edgeSides.size()) << "Sizes do not match";
			const int vCentral = ED.E(edges[0], edgeSides[0]);
			ASSERT_EQ(Vvisited(vCentral), 0);
			SCOPED_TRACE(std::string("Central vert:") + std::to_string(vCentral));
			for (int i = 0; i < edges.size() / 2; i++)
			{
				const int spoke = 2 * i;
				const int ring = 2 * i + 1;
				// Test connection to central vertex
				int v = startVert(spoke, edges, edgeSides);
				EXPECT_TRUE(v == vCentral) << "Invalid spoke start vert, expected " << vCentral << " got " << v;
				// Test connection between spoke and ring
				v = endVert(spoke, edges, edgeSides);
				int v2 = startVert(ring, edges, edgeSides);
				EXPECT_TRUE(v == v2) << "Spoke and ring do not meet in vertex: " << v << " vs " << v2;
				// Test common face
				int f1 = lface(spoke, edges, edgeSides), f2 = lface(ring, edges, edgeSides);
				EXPECT_TRUE(f1 == f2) << "Spoke and ring do not share face at right side: " << f1 << " vs " << f2;
				// Check connection to next pair
				const int nextSpoke = (2 * i + 2) % edges.size();
				f2 = rface(nextSpoke, edges, edgeSides);
				EXPECT_TRUE(f1 == f2) << "Next spoke is not properly related via face: " << f1 << " vs " << f2;
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
					EXPECT_EQ(ED.F(f, (start + i) % 3), vs[i]);
				}

			}
			Vvisited(vCentral)++;
		}
		void handleBoundaryRing(const std::vector<int>& edges, const std::vector<int>& edgeSides)
		{
			if (ED.boundaryEdgeCount == 0) // handleBoundaryRing() should only be called when there is an actual boundary
				FAIL() << "Not expecting boundary";
			else
			{
				// First element and last element are boundary
				EXPECT_EQ(ED.EF(edges[0], 1 - edgeSides[0]), -1);
				EXPECT_EQ(ED.EF(edges.back(), edgeSides.back()), -1);

				// Descriptors for edges should be of same size
				EXPECT_TRUE(edges.size() == edgeSides.size()) << "Sizes do not match";

				const int vCentral = ED.E(edges[0], edgeSides[0]);
				// Expect the central vertex to not be visited yet
				EXPECT_EQ(Vvisited(vCentral), 0);
				SCOPED_TRACE(std::string("Central vert boundary:") + std::to_string(vCentral));
				// Check the valence 
				EXPECT_EQ(valences(vCentral), (edges.size() + 1) / 2);

				// Update visisted edge at end of ring
				Evisited(edges.back(), edgeSides.back()) += 1;

				BoundaryVs(vCentral) -= 1;

				for (int i = 0; i < (edges.size() - 1) / 2; i++)
				{
					const int spoke = 2 * i;
					const int ring = 2 * i + 1;

					// Test connection to central vertex
					int v = startVert(spoke, edges, edgeSides);
					EXPECT_TRUE(v == vCentral) << "Invalid spoke start vert, expected " << vCentral << " got " << v;

					// Test connection between spoke and ring
					v = endVert(spoke, edges, edgeSides);
					int v2 = startVert(ring, edges, edgeSides);
					EXPECT_TRUE(v == v2) << "Spoke and ring do not meet in vertex: " << v << " vs " << v2;

					// Test common face
					int f1 = lface(spoke, edges, edgeSides), f2 = lface(ring, edges, edgeSides);
					EXPECT_TRUE(f1 == f2) << "Spoke and ring do not share face at right side: " << f1 << " vs " << f2;

					// Check connection to next pair
					const int nextSpoke = (2 * i + 2) % edges.size();
					f2 = rface(nextSpoke, edges, edgeSides);
					EXPECT_TRUE(f1 == f2) << "Next spoke is not properly related via face: " << f1 << " vs " << f2;
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
						EXPECT_EQ(ED.F(f, (start + i) % 3), vs[i]);
					}

				}
				Vvisited(vCentral)++;
			}
		}
		void checkVisit()
		{
			// Check all vertices were visited
			for (int i = 0; i < Vvisited.rows(); i++)
			{
				EXPECT_TRUE(Vvisited(i) == 1) << std::string("Missed vertex ") + std::to_string(i);
			}
			for (int i = 0; i < Evisited.rows(); i++)
			{
				if (ED.EF(i, 0) == -1)
				{
					EXPECT_TRUE(Evisited(i, 0) == 1) << "Expected to see edge " << i << " once at orient 0, but saw " << Evisited(i, 0);
				}
				else
				{
					EXPECT_TRUE(Evisited(i, 0) == 2) << "Expected to see edge " << i << " twice at orient 0, but saw " << Evisited(i, 0);
				}
				if (ED.EF(i, 1) == -1)
				{
					EXPECT_TRUE(Evisited(i, 1) == 1) << "Expected to see edge " << i << " once at orient 1, but saw " << Evisited(i, 1);
				}
				else
				{
					EXPECT_TRUE(Evisited(i, 1) == 2) << "Expected to see edge " << i << " twice at orient 1, but saw " << Evisited(i, 1);
				}

			}
			for (int i = 0; i < BoundaryVs.rows(); i++)
			{
				EXPECT_TRUE(BoundaryVs(i) == 0) << "Invalid boundary : at " << i << ": " << BoundaryVs(i);
			}
		}
	};
}

#endif