#ifndef DIRECTIONAL_SHM_EDGE_TOPOLOGY_H
#define DIRECTIONAL_SHM_EDGE_TOPOLOGY_H
#include <Eigen/Eigen>
#include <igl/sortrows.h>
namespace directional{
    void shm_edge_topology(
        const Eigen::MatrixXi& F,
        const int& vertexCount,
        Eigen::MatrixXi& E,
        Eigen::MatrixXi& EF,
        Eigen::MatrixXi& EI,
        Eigen::MatrixXi& SFE
    )
	{
		assert(F.cols() == 3, "Only triangle meshes supported");
		SFE.resize(F.rows(), 6);
		Eigen::MatrixXi rawEdges(3 * F.rows(), 2);
		for (int f = 0; f < F.rows(); f++)
		{
			SFE.row(f) = Eigen::RowVectorXi::Constant(6, -1);
			const auto& row = F.row(f);
			rawEdges.row(3 * f) = Eigen::RowVector2i(row(1), row(2));
			rawEdges.row(3 * f + 1) = Eigen::RowVector2i(row(2), row(0));
			rawEdges.row(3 * f + 2) = Eigen::RowVector2i(row(0), row(1));
			// Sort rows ascending. Record swap as opposite orientation.
			for (int i = 0; i < 3; i++) {
				if (rawEdges(3 * f + i, 0) > rawEdges(3 * f + i, 1)) {
					std::swap(rawEdges(3 * f + i, 0), rawEdges(3 * f + i, 1));
					SFE(f, 3 + i) = 1;
				}
				else
				{
					SFE(f, 3 + i) = 0;
				}
			}
		}
		// Sort total rows by first vertex index. Same vertex edges should now be adjacent rows in the matrix
		Eigen::VectorXi edgeMap;
		igl::sortrows(rawEdges, true, E, edgeMap);
		EF = Eigen::MatrixXi::Constant(edgeMap.rows(), 2,-1);
		EI = Eigen::MatrixXi::Constant(edgeMap.rows(), 2, -1);
		// Deduplicate edges and construct other matrices.
		int currId = 0;
		for (int i = 0; i < edgeMap.rows(); i++)
		{
			const int f = edgeMap(i) / 3;
			const int corner = edgeMap(i) - 3 * f;
			// Edge is not same as previous edge.
			if (E.row(i) != E.row(currId))
			{
				currId++;
				// Assign the row to the appropriate row in E.
				E.row(currId) = E.row(i);
			}
			// Update connectivity data.
			SFE(f, corner) = currId;
			EF(currId, SFE(f, 3 + corner)) = f;
			EI(currId, SFE(f, 3 + corner)) = corner;
		}
		currId++;

		// Resize the edge matrices to the appropriate size
		E.conservativeResize(currId, 2);
		EF.conservativeResize(currId, 2);
		EI.conservativeResize(currId, 2);

		//Acquire the boundary edge count
		// boundaryEdgeCount = 0;
		for (int e = 0; e < currId; e++)
		{
			// if (EF(e, 0) == -1 || EF(e, 1) == -1)boundaryEdgeCount++;
			// Fixup boundary to have outside at left of edge
			if (EF(e, 1) == -1)
			{
                // Flip edge direction
				std::swap(E(e, 0), E(e, 1));
				SFE(EF(e, 0), EI(e, 0) + 3) = 1 - SFE(EF(e, 0), EI(e, 0) + 3);
                std::swap(EF(e, 0), EF(e, 1));
                std::swap(EI(e, 0), EI(e, 1));
			}
		}
	}
}
#endif