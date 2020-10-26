#ifndef DIRECTIONAL_SHM_EDGE_TOPOLOGY_H
#define DIRECTIONAL_SHM_EDGE_TOPOLOGY_H
#include <Eigen/Eigen>
#include <igl/sortrows.h>
namespace directional{
	
    inline void shm_edge_topology(
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& EF,
        Eigen::MatrixXi& EI,
        Eigen::MatrixXi& SFE
    )
    {
        SFE = Eigen::MatrixXi::Constant(F.rows(), 6, -1);
        EI = Eigen::MatrixXi::Constant(EF.rows(), 2,-1);
        for (int e = 0; e < EF.rows(); ++e)
        {
            auto leftF = EF(e, 0);
            if(leftF != -1)
            {
                for (int j = 0; j < 3; ++j)
                {
                    auto v = F(leftF, j);
                    if (E(e, 0) != v && E(e, 1) != v)
                    {
                        SFE(leftF, j) = e;
                        SFE(leftF, 3 + j) = 0;
                        EI(e, 0) = j;
                    }
                }
            }
            
            auto rightF = EF(e, 1);
            if(rightF != -1)
            {
                for (int j = 0; j < 3; ++j)
                {
                    auto v = F(rightF, j);
                    if (E(e, 0) != v && E(e, 1) != v)
                    {
                        SFE(rightF, j) = e;
                        SFE(rightF, 3 + j) = 1;
                        EI(e, 1) = j;
                    }
                }
            }
        }
    }
    
    inline void shm_edge_topology(
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& E,
        Eigen::MatrixXi& EF,
        Eigen::MatrixXi& EI,
        Eigen::MatrixXi& SFE
    )
    {
        SFE = Eigen::MatrixXi::Constant(F.rows(), 6, -1);
        EI = Eigen::MatrixXi::Constant(E.rows(), 2, -1);
        EF = Eigen::MatrixXi::Constant(E.rows(), 2, -1);

        // Reconstruct EF
        {
            std::map<int,std::set<int>> VToF;
            for(auto f= 0; f < F.rows(); ++f)
            {
                for(int c = 0; c < 3; ++c)
                {
                    if (VToF.find(F(f, c)) == VToF.end()) VToF[F(f, c)] = {};
                    VToF[F(f, c)].insert(f);
                }
            }
            for(int e = 0; e < E.rows(); ++e)
            {
                const auto v0 = E(e, 0);
                const auto v1 = E(e, 1);
                std::vector<int> diff;
                auto end = std::set_intersection(VToF[v0].begin(), VToF[v0].end(), VToF[v1].begin(), VToF[v1].end(), diff.begin());
                assert(std::distance(diff.begin(), end) == 2);
                auto f0 = diff[0];
                auto f1 = diff[1];
                // Determine left and right
                for(int c = 0;c < 3; ++c)
                {
                    if(F(f0,c) == v0)
                    {
                        if(F(f0,(c+1)%3) == v1)
                        {
                            EF(e, 0) = f0;
                            EF(e, 1) = f1;
                        }
                        else
                        {
                            EF(e, 0) = f1;
                            EF(e, 1) = f0;
                        }
                        break;
                    }
                }
            }
        }

        for (int e = 0; e < EF.rows(); ++e)
        {
            auto leftF = EF(e, 0);
            if (leftF != -1)
            {
                for (int j = 0; j < 3; ++j)
                {
                    auto v = F(leftF, j);
                    if (E(e, 0) != v && E(e, 1) != v)
                    {
                        SFE(leftF, j) = e;
                        SFE(leftF, 3 + j) = 0;
                        EI(e, 0) = j;
                    }
                }
            }

            auto rightF = EF(e, 1);
            if (rightF != -1)
            {
                for (int j = 0; j < 3; ++j)
                {
                    auto v = F(rightF, j);
                    if (E(e, 0) != v && E(e, 1) != v)
                    {
                        SFE(rightF, j) = e;
                        SFE(rightF, 3 + j) = 1;
                        EI(e, 1) = j;
                    }
                }
            }
        }
    }

    inline void shm_edge_topology(
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& EF,
        const Eigen::MatrixXi& FE,
        bool isEdgeTopologyBased,
        Eigen::MatrixXi& EI,
        Eigen::MatrixXi& SFE
    )
    {
        SFE = Eigen::MatrixXi(EF.rows(), 6);
        for(int f = 0; f < FE.rows(); ++f)
        {
            // SHM follows edge_flaps based stuff, so compensate
            // Difference: edge FE(f,0) of edgetopology is between F(f,0) and F(f,1)
            // instead of opposite corner F(f,0), i.e. between vertices F(f,1) and F(f,2) (possibly flipped)
            if(isEdgeTopologyBased)
            {
                SFE(f, 2) = FE(f, 0);
                SFE(f, 0) = FE(f, 1);
                SFE(f, 1) = FE(f, 2);
            }
            else
            {
                SFE.block(f, 0, 1, 3) = FE.row(f);
            }
            for(int j = 0; j < 3; ++j)
            {
                SFE(f, 3 + j) = F(f, (j + 1) % 3) == E(SFE(f, j), 0) ? 0 : 1;
                // Update EI
                EI(SFE(f, j), SFE(f, 3 + j)) = j;
            }
        }
    }

    /**
     * TODO remove this in favour of centralized approach to edge topology
     * Converts the edge topology as used in the subdivision construction to igl::edge_topology compliant matrices.
     * The main difference is that the face-to-edge connectivity in igl::edge_topology makes it so that edge FE(f,0) 
     * is the edge between vertices F(f,0) and F(f,1), whereas the subdivision topology uses SFE, where edge SFE(f,0)
     * is the edge between vertices F(f,1,) and F(f,2), i.e. opposite corner (vertex) 0 of the face.
     * 
     * Input:
     * - F |F| x 3 matrix of face to vertex connectivity
     * - E |E| x 2 matrix of edge to vertex connectivity, with edge e directed from vertex E(e,0) to E(e,1).
     * - EF |E| x 2 matrix of edge to face connectivity. Given the direction as described by matrix E, EF(e,0) is the face
     * left of edge e and EF(e,1) the face to the right
     * - SFE |F| x 6 matrix of face to edge connectivity, such that edge SFE(f,c) is opposite corner c in face f. In addition,
     * contains orientations of the edges, such that edge SFE(f,c) is CCW with respect to f if SFE(f,c+3) = 0, and CW if this is 1.
     * Output:
     * - EI |E| x 2 matrix such that for edge e, EI(e,i) is the index of edge e in face EF(e,i). Note that this is different for
     * shm topology and igl::edge_topology, since this depends on FE/SFE.
     * - FE |F| x 3 matrix of face to edge connectivity such that edge FE(f,c) is between vertices F(f,c) and F(f,(c+1)%3).  
     */
    inline void shm_edge_topology_to_igledgetopology(
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& EF,
        const Eigen::MatrixXi& SFE,
        Eigen::MatrixXi& EI,
        Eigen::MatrixXi& FE
    )
    {
        FE = Eigen::MatrixXi(EF.rows(), 3);
        EI = Eigen::MatrixXi(EF.rows(), 2);
        for (int f = 0; f < FE.rows(); ++f)
        {
            for(int j = 0; j <3; ++j)
            {
                auto e = SFE(f, j);
                FE(f, (j + 1) % 3) = e;
                EI(e, EF(e, 0) == f ? 0 : 1) = (j + 1) % 3;
            }
        }
    }

    inline void shm_edge_topology(
        const Eigen::MatrixXi& F,
        const int& vertexCount,
        Eigen::MatrixXi& E,
        Eigen::MatrixXi& EF,
        Eigen::MatrixXi& EI,
        Eigen::MatrixXi& SFE
    )
	{
		assert((F.cols() == 3) &&  "Only triangle meshes supported");
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
