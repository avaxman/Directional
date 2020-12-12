#ifndef DIRECTIONAL_TEST_INCLUDES_H
#define DIRECTIONAL_TEST_INCLUDES_H
#include <Eigen/Eigen>
#include <igl/edge_topology.h>
#include <igl/read_triangle_mesh.h>

struct FEM_operators
{
    Eigen::SparseMatrix<double>Gv, Ge, J, C, D;
};

/**
 * \brief 
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
    if(dim == 1)
    {
        assert(mat.rows() % elementsPerBlock == 0);
        output = Eigen::SparseMatrix<double>(mat.rows(), mat.cols() * mat.rows() / elementsPerBlock);
        for(auto i = 0; i < mat.rows() / elementsPerBlock; ++i)
        {
            for(auto j =0; j < elementsPerBlock; ++j)
            {
                Eigen::RowVectorXd row = mat.row(i*elementsPerBlock + j);
                for(int c = 0; c < mat.cols(); ++c)
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
                    trips.emplace_back(i*mat.rows() + r, i * elementsPerBlock + j,col(r));
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
    TriangleMesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F): V(V),F(F){}
    TriangleMesh(){}

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
        output.setConstant(2 * F.rows(), 3,0);
        Eigen::MatrixXd N;
        igl::per_face_normals(V, F, N);
        for(auto f =0; f < F.rows();++f)
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
        igl::read_triangle_mesh(DIRECTIONAL_TEST_MODELS_DIR + std::string("/")+ fileName, V, F);
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
        FEM_suite(operators.Gv, operators.Ge, operators.J, operators.C, operators.D);
    }
};

#endif