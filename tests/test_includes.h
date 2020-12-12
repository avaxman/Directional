#ifndef DIRECTIONAL_TEST_INCLUDES_H
#define DIRECTIONAL_TEST_INCLUDES_H
#include <Eigen/Eigen>
#include <igl/edge_topology.h>
#include <igl/read_triangle_mesh.h>

namespace directional_testing
{

    struct TriangleMesh;

    struct FEM_operators
    {
        Eigen::SparseMatrix<double>Gv, Ge, J, C, D;

        const TriangleMesh* mesh = nullptr;

        FEM_operators(){}
        FEM_operators(const TriangleMesh& mesh):mesh(&mesh){}

        inline void eliminateBoundary(const Eigen::SparseMatrix<double>& edgeOperator,
                                     Eigen::SparseMatrix<double>& newOperator);


        void adjacencyMatrix(const Eigen::MatrixXi& adjacency, Eigen::SparseMatrix<double>& adjacencyMatrix)
        {
            adjacencyMatrix = Eigen::SparseMatrix<double>(adjacency.rows(), adjacency.maxCoeff() + 1);
            std::vector<Eigen::Triplet<double>> trips;
            trips.reserve(adjacency.rows() * adjacency.cols());
            for (auto i = 0; i < adjacency.rows(); ++i)
            {
                for (auto j = 0; j < adjacency.cols(); ++j)
                {
                    if (adjacency(i, j) < 0) continue;
                    trips.emplace_back(i, adjacency(i, j), 1.0);
                }
            }
            adjacencyMatrix.setFromTriplets(trips.begin(), trips.end());
        }

        inline void shm_divergence(Eigen::SparseMatrix<double>& shmD);

        inline void eliminateBoundaryOp(Eigen::SparseMatrix<double>& elim);

        void dec_D1(const directional_testing::TriangleMesh& data, Eigen::SparseMatrix<double>& d1);
    };
    /**
     * \brief Creates a sparse block matrix out of a regular matrix, grouping a fixed number of rows/columns to form a block in the
     * matrix.
     * \param mat Matrix to convert to blocks.
     * \param dim Dimension to use as elements. 1 = column direction, element is a row. 2 = row direction, element is a column.
     * \param elementsPerBlock Number of elements per block. Suppose dim=1, then every `elementsPerBlock` rows will form a block.
     * Note that we have to have that if dim=1, mat.rows() % elementsPerBlock == 0, or when dim=2, mat.cols() % elementsPerBlock == 0.
     * \param output The output block matrix
     */
    inline void toBlocks(const Eigen::MatrixXd& mat, int dim, int elementsPerBlock, Eigen::SparseMatrix<double>& output)
    {
        std::vector<Eigen::Triplet<double>> trips;
        trips.reserve(mat.rows() * mat.cols());
        if (dim == 1)
        {
            assert(mat.rows() % elementsPerBlock == 0);
            output = Eigen::SparseMatrix<double>(mat.rows(), mat.cols() * mat.rows() / elementsPerBlock);
            for (auto i = 0; i < mat.rows() / elementsPerBlock; ++i)
            {
                for (auto j = 0; j < elementsPerBlock; ++j)
                {
                    Eigen::RowVectorXd row = mat.row(i*elementsPerBlock + j);
                    for (int c = 0; c < mat.cols(); ++c)
                    {
                        trips.emplace_back(i * elementsPerBlock + j, c + i * mat.cols(), row(c));
                    }
                }
            }
        }
        else
        {
            assert(mat.cols() % elementsPerBlock == 0);
            output = Eigen::SparseMatrix<double>(mat.cols() * mat.rows() / elementsPerBlock, mat.cols());
            for (auto i = 0; i < mat.cols() / elementsPerBlock; ++i)
            {
                for (auto j = 0; j < elementsPerBlock; ++j)
                {
                    Eigen::VectorXd col = mat.col(i*elementsPerBlock + j);
                    for (int r = 0; r < mat.rows(); ++r)
                    {
                        trips.emplace_back(i*mat.rows() + r, i * elementsPerBlock + j, col(r));
                    }
                }
            }
        }
        output.setFromTriplets(trips.begin(), trips.end());
    }

    struct TriangleMesh
    {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F, E, EF, FE;
        TriangleMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F) : V(V), F(F) {}
        TriangleMesh() {}

        void setMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
        {
            this->V = V;
            this->F = F;
        }

        /**
         * \brief Outputs an orthogonal basis for the mesh to define a PCVF on its surface.
         * Will be returned as a matrix of a vector per row, with two rows 2*i, 2*i+1 being E0 and the pi/2 rotated E0 of face i.
         * \param output The basis output
         */
        void basis(Eigen::MatrixXd& output)
        {
            output.setConstant(2 * F.rows(), 3, 0);
            Eigen::MatrixXd N;
            igl::per_face_normals(V, F, N);
            for (auto f = 0; f < F.rows(); ++f)
            {
                Eigen::RowVector3d E = V.row(F(f, 1)) - V.row(F(f, 1));
                Eigen::RowVector3d normal = N.row(f);
                Eigen::RowVector3d EPerp = normal.cross(E);
                output.row(2 * f) = E;
                output.row(2 * f + 1) = EPerp;
            }
        }

        /**
         * \brief Compute IGL edge topology, setting E, FE and EF, based on V and F.
         */
        void compute_edge_topology()
        {
            assert(V.rows() > 0 && F.rows() > 0);
            igl::edge_topology(V, F, E, FE, EF);
        }
        /**
         * \brief Compute EF and FE, based on F and E. Uses igl::edge_topology standard for EI/FE.
         */
        void compute_edge_topology_fixed_E()
        {
            assert(V.rows() > 0 && F.rows() > 0 && E.rows() > 0);
            Eigen::MatrixXi EI, SFE;
            directional::shm_edge_topology(F, std::cref(E), EF, EI, SFE);
            directional::shm_edge_topology_to_igledgetopology(F, E, EF, SFE, EI, FE);
        }
        /**
         * \brief Reads the file in the tests models directory, using igl::read_triangle_mesh, and sets V and F
         * \param fileName The path to the model file
         */
        void read(const std::string& fileName)
        {
            igl::read_triangle_mesh(DIRECTIONAL_TEST_MODELS_DIR + std::string("/") + fileName, V, F);
        }

        /**
         * \brief Computes the FEM suite via directional::FEM_suite using the mesh data.
         * \param Gv
         * \param Ge
         * \param J
         * \param C
         * \param D
         */
        void FEM_suite(Eigen::SparseMatrix<double>& Gv, Eigen::SparseMatrix<double>& Ge, Eigen::SparseMatrix<double>& J, Eigen::SparseMatrix<double>& C,
            Eigen::SparseMatrix<double>& D)
        {
            assert(V.rows() > 0 && F.rows() > 0 && E.rows() > 0 && EF.rows() > 0 && FE.rows() > 0);
            directional::FEM_suite(V, F, E, FE, EF, Gv, Ge, J, C, D);
        }
        void FEM_suite(FEM_operators& operators)
        {
            operators.mesh = this;
            FEM_suite(operators.Gv, operators.Ge, operators.J, operators.C, operators.D);
        }
        void boundary_vertices(std::set<int>& vertices)
        {
            for(auto e = 0; e < EF.rows(); ++e)
            {
                if(EF(e,0) == -1 || EF(e,1) == -1)
                {
                    vertices.insert(E(e, 0));
                    vertices.insert(E(e, 1));
                }
            }
        }
        void boundary_faces(std::set<int>& faces)
        {
            for (auto e = 0; e < EF.rows(); ++e)
            {
                if (EF(e, 0) == -1 || EF(e, 1) == -1)
                {
                    faces.insert(EF(e, 0) == -1 ? EF(e, 1) : EF(e, 0));
                }
            }
        }
        void boundary_edge(std::set<int>& edges)
        {
            for (auto e = 0; e < EF.rows(); ++e)
            {
                if (EF(e, 0) == -1 || EF(e, 1) == -1)
                {
                    edges.insert(e);
                }
            }
        }
        void boundary_vertex_edge_mapping(std::vector<std::vector<std::pair<int,int>>>& mapping)
        {
            mapping.resize(V.rows(), {});

            for (auto e = 0; e < EF.rows(); ++e)
            {
                if (EF(e, 0) == -1 || EF(e, 1) == -1)
                {
                    mapping[E(e, 0)].push_back(std::make_pair(e, 0));
                    mapping[E(e, 1)].push_back(std::make_pair(e, 1));
                }
            }
        }
    };


    inline void FEM_operators::eliminateBoundary(const Eigen::SparseMatrix<double>& edgeOperator,
        Eigen::SparseMatrix<double>& newOperator)
    {
        Eigen::SparseMatrix<double> elim(edgeOperator.rows(), edgeOperator.rows());
        Eigen::VectorXd diag;
        diag.setConstant(edgeOperator.rows(), 1.0);
        for (auto e = 0; e < mesh->EF.rows(); ++e)
        {
            // Boundary: one of the faces of the edge is missing
            if (mesh->EF(e, 0) == -1 || mesh->EF(e, 1) == -1)
            {
                diag(e) = 0;
            }
        }
        elim = diag.asDiagonal();
        newOperator = elim * edgeOperator;
    }
    inline void FEM_operators::eliminateBoundaryOp(Eigen::SparseMatrix<double>& elim)
    {
        elim = Eigen::SparseMatrix<double>(mesh->EF.rows(), mesh->EF.rows());
        Eigen::VectorXd diag;
        diag.setConstant(mesh->EF.rows(), 1.0);
        for (auto e = 0; e < mesh->EF.rows(); ++e)
        {
            // Boundary: one of the faces of the edge is missing
            if (mesh->EF(e, 0) == -1 || mesh->EF(e, 1) == -1)
            {
                diag(e) = 0;
            }
        }
        elim = diag.asDiagonal();
    }

    inline void FEM_operators::shm_divergence(Eigen::SparseMatrix<double>& shmD)
    {
        Eigen::SparseMatrix<double> add(D.rows(), D.cols());
        std::vector<Eigen::Triplet<double>> trips;
        Eigen::VectorXd DA;
        igl::doublearea(mesh->V, mesh->F, DA);

        Eigen::MatrixXd N;
        igl::per_face_normals(mesh->V, mesh->F, N);

        for (auto e = 0; e < mesh->EF.rows(); ++e)
        {
            // Only apply to boundary edges
            if (mesh->EF(e, 0) == -1 || mesh->EF(e, 1) == -1)
            {
                int f = mesh->EF(e, 0) == -1 ? mesh->EF(e, 1) : mesh->EF(e, 0);
                // Find vertex opposite edge
                int c = 0;
                for(; c<3; ++c)
                {
                    // Vertex is not in edge
                    if(mesh->F(f, c) != mesh->E(e,0) && mesh->F(f, c) != mesh->E(e, 1))
                    {
                        break;
                    }
                }
                double sign = mesh->E(e, 0) == mesh->F(f, (c + 1) % 2) ? 1.0 : -1.0;

                Eigen::RowVector3d edgeDir = mesh->V.row(mesh->F(f, (c + 2) % 3)) - mesh->V.row(mesh->F(f, (c + 1) % 3));
                Eigen::RowVector3d edgeDirPerp = Eigen::RowVector3d(N.row(f)).cross(edgeDir);
                // Now, project on negated gradient direction times face area for the face and sum to the vertices.
                for(int i = 0; i < 3; ++i)
                {
                    trips.emplace_back(mesh->E(e, 0), 3 * f + i, 1.0 * edgeDirPerp(i) / 2.0 * sign);
                    trips.emplace_back(mesh->E(e, 1), 3 * f + i, -1.0 * edgeDirPerp(i) / 2.0 * sign);
                }
            }
        }
        add.setFromTriplets(trips.begin(), trips.end());
        shmD = D + add;
    }

    inline void FEM_operators::dec_D1(const directional_testing::TriangleMesh& data, Eigen::SparseMatrix<double>& d1)
    {
        d1 = Eigen::SparseMatrix<double>(data.F.rows(), data.E.rows());
        std::vector<Eigen::Triplet<double>> trips;
        for (auto f = 0; f < data.FE.rows(); ++f)
        {
            for (auto j = 0; j < data.FE.cols(); ++j)
            {
                const auto e = data.FE(f, j);
                double v = 0;
                for (int c = 0; c < 3; ++c)
                {
                    if (data.F(f, c) == data.E(e, 0))
                    {
                        if (data.F(f, (c + 1) % 3) == data.E(e, 1))
                        {
                            v = 1.0;
                        }
                        else
                        {
                            v = -1.0;
                        }
                    }
                }
                trips.emplace_back(f, e, v);
            }
        }
        d1.setFromTriplets(trips.begin(), trips.end());
    }
}


#endif