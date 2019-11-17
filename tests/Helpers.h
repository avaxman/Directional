#ifndef DIRECTIONAL_TEST_HELPERS_H
#define DIRECTIONAL_TEST_HELPERS_H
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <gtest/gtest.h>

// Gtest does not support testing outside of test class...
//#define DIR_ASSERT(x) ASSERT_TRUE(x) << "Failed condition " << #x  << " at " << __FILE__;
//#define DIR_ASSERT_M(x, m) ASSERT_TRUE(x) << (m);
// So, we will just disable it then...
#define DIR_ASSERT(x) 
#define DIR_ASSERT_M(x, m) 

#include <igl/read_triangle_mesh.h>
#include <directional/Gamma_suite.h>
#include <vector>
#include <tuple>
#include <directional/Subdivision/EdgeData.h>

namespace directional_fixtures
{
	struct Mesh
	{
		Eigen::MatrixXd V;
		Eigen::MatrixXi F;
		std::string fileName;

		EdgeData getED()
		{
			return EdgeData(F);
		}
		int faceCount() const
		{
			return F.rows();
		}
		int vertCount() const
		{
			return V.rows();
		}
	};
	/**
	 * \brief Mesh test with parameter the file with a mesh.
	 */
	class MeshTestFixture : public testing::TestWithParam<const char*>
	{
	protected:
		Mesh m_mesh;
		void SetUp() override
		{
			// Resolve against shared path
			m_mesh.fileName = std::string(TEST_SHARED_PATH) + "/";
			m_mesh.fileName = m_mesh.fileName + GetParam();
			bool success = igl::read_triangle_mesh(m_mesh.fileName, m_mesh.V, m_mesh.F);
			ASSERT_TRUE(success);
			ASSERT_GT(m_mesh.F.rows(), 0) << "No faces in mesh";
		}
		void TearDown() override
		{
			
		}
	};
}

namespace helpers{

	inline std::string description(const Eigen::SparseMatrix<double>& mat)
	{
		std::stringstream str;
		str << "Cols:" << mat.cols() << ", rows: " << mat.rows();
		return str.str();
	}

	struct FEM_Ops
	{
		using SparseMat = Eigen::SparseMatrix<double>;
		SparseMat Gv, Ge, J, C, D;

		FEM_Ops(const directional_fixtures::Mesh& m)
		{
			EdgeData ED;
			directional::construct_edge_topology(m.F, ED.E, ED.EF, ED.sFE, ED.EI);
			directional::FEM_suite(m.V, m.F, ED.E, ED.sFE, ED.EF, Gv, Ge, J, C, D);
		}
		FEM_Ops() {}
		void constructOperators(const EdgeData& ED, const directional_fixtures::Mesh& m)
		{
			directional::FEM_suite(m.V, m.F, ED.E, ED.sFE, ED.EF, Gv, Ge, J, C, D);
		}
	};

	/**
	 * \brief Helper struct that constructs all gamma2 related operators
	 */
	struct Gamma2_Ops
	{
		using SparseMat = Eigen::SparseMatrix<double>;
		SparseMat Gv, Ge, J, C, D, Gamma2_To_Oneform, G2_To_Decomp, Decomp_To_G2, Chi_To_Gamma2, Gamma2_To_Chi,
			G2ToG3, G3ToG2;

		Gamma2_Ops(const directional_fixtures::Mesh& m)
		{
			EdgeData ED(m.F);
			constructOperators(ED, m);
		}
		Gamma2_Ops(const EdgeData& ED, const Eigen::MatrixXd& V)
		{
			directional_fixtures::Mesh m;
			m.F = ED.F;
			m.V = V;
			constructOperators(ED, m);
		}
		Gamma2_Ops() {}
		void constructOperators(const EdgeData& ED, const directional_fixtures::Mesh& m)
		{
			directional::Gamma2_suite(m.V, m.F, ED.E, ED.EI, ED.sFE, ED.EF, Gv, Ge, J, C, D, Gamma2_To_Oneform, G2_To_Decomp,Decomp_To_G2);
			// Gamma2 to decmposition into oneform and curl component (average and half curl of the gamma's).
			directional::Gamma2_projector(m.V, m.F, ED.E, ED.sFE, ED.EF, Chi_To_Gamma2);
			directional::Gamma2_reprojector(m.V, m.F, ED.E, ED.sFE, ED.EF, Gamma2_To_Chi);

			directional::Gamma2_To_Gamma3(ED.sFE, G2ToG3);
			directional::Gamma3_To_Gamma2(ED.sFE, G3ToG2);
		}
	};

	/**
	 * \brief Helper struct that constructs all gamma2 related operators
	 */
	struct Gamma3_Ops
	{
		using SparseMat = Eigen::SparseMatrix<double>;
		SparseMat Gv, Ge, J, C, D, Gamma3_To_Oneform, G3_To_Decomp, Decomp_To_G3, Chi_To_Gamma3, Gamma3_To_Chi;

		Gamma3_Ops(const directional_fixtures::Mesh& m)
		{
			EdgeData ED;
			directional::construct_edge_topology(m.F, ED.E, ED.EF, ED.sFE, ED.EI);
			ED.F = m.F;
			constructOperators(ED, m);
		}
		Gamma3_Ops(const EdgeData& ED, const Eigen::MatrixXd& V)
		{
			directional_fixtures::Mesh m;
			m.F = ED.F;
			m.V = V;
			constructOperators(ED, m);
		}
		Gamma3_Ops() {}
		void constructOperators(const EdgeData& ED, const directional_fixtures::Mesh& m)
		{
			SparseMat G2ToG3, G3ToG2;
			directional::Gamma2_To_Gamma3(ED.sFE, G2ToG3);
			directional::Gamma3_To_Gamma2(ED.sFE, G3ToG2);
			directional::Gamma3_suite(m.V, m.F, ED.E,ED.EI, ED.sFE, ED.EF, Gv, Ge, J, C, D, Gamma3_To_Oneform, G3_To_Decomp, Decomp_To_G3);
			// Gamma2 to decmposition into oneform and curl component (average and half curl of the gamma's).
			directional::Gamma2_projector(m.V, m.F, ED.E, ED.sFE, ED.EF, Chi_To_Gamma3);
			Chi_To_Gamma3 = G2ToG3 * Chi_To_Gamma3;
			directional::Gamma2_reprojector(m.V, m.F, ED.E, ED.sFE, ED.EF, Gamma3_To_Chi);
			Gamma3_To_Chi = Gamma3_To_Chi * G3ToG2;

		}
	};

	std::ostream& operator<<(std::ostream& str, const directional_fixtures::Mesh& m)
	{
		str << m.fileName;
		return str;
	}
#ifdef WIN32
#define MD_SEP "/"
#else
#define MD_SEP "\\"
#endif

	inline Eigen::VectorXd random_edgefield()
	{
		
	}

	template<typename DataType>
	inline DataType maxAbs(const Eigen::SparseMatrix<DataType>& mat)
	{
		DataType currMax = -1;
		using SMat = Eigen::SparseMatrix<DataType>;
		for (int k = 0; k < mat.outerSize(); ++k) {
			for (typename SMat::InnerIterator it(mat, k); it; ++it) {
				currMax = std::max(std::abs(it.value()), currMax);
			}
		}
		return currMax;
	}


     template<typename DataType, int Opts>
     inline bool isApproximatelyDiagonal(const Eigen::SparseMatrix<DataType, Opts>& mat, DataType threshold){
		 // Non square = not diagonal.
		 if (mat.rows() != mat.cols()) return false;
		 using SMat = Eigen::SparseMatrix<DataType, Opts>;
		 for (int k = 0; k < mat.outerSize(); ++k) {
			 for (typename SMat::InnerIterator it(mat, k); it; ++it) {
				 if (it.row() != it.col() && std::abs(it.value()) > threshold) return false;
			 }
		 }
		 return true;
     }

	 template<typename DataType, int Opts>
	 inline bool isApproximatelyIdentity(const Eigen::SparseMatrix<DataType, Opts>& mat, DataType threshold) {
		 // Non square = not identity.
		 if (mat.rows() != mat.cols()) return false;
		 using SMat = Eigen::SparseMatrix<DataType, Opts>;
		 for (int k = 0; k < mat.outerSize(); ++k) {
			 for (typename SMat::InnerIterator it(mat, k); it; ++it) {
				 if (it.row() != it.col() && std::abs(it.value()) > threshold) return false;
				 else if (it.row() == it.col() && std::abs(it.value() - 1) > threshold) return false;
			 }
		 }
		 return true;
	 }
	 template<typename DataType, int Opts>
	 inline std::vector<std::tuple<int,int,double>> getDifferenceFromIdentity(const Eigen::SparseMatrix<DataType, Opts>& mat, DataType threshold) {
		 // Non square = not identity.
		 if (mat.rows() != mat.cols()) return {{-1,-1,-1}};
		 std::vector<std::tuple<int, int, double>> ret;
		 using SMat = Eigen::SparseMatrix<DataType, Opts>;
		 for (int k = 0; k < mat.outerSize(); ++k) {
			 for (typename SMat::InnerIterator it(mat, k); it; ++it) {
				 if (it.row() != it.col() && std::abs(it.value()) > threshold) ret.emplace_back(it.row(),it.col(), it.value());
				 else if (it.row() == it.col() && std::abs(it.value() - 1) > threshold) ret.emplace_back(it.row(), it.col(), it.value());
			 }
		 }
		 return ret;
	 }

	 template<typename DataType, int Opts>
	 inline std::vector<std::tuple<int, int, double>> getDifferenceFromZero(const Eigen::SparseMatrix<DataType, Opts>& mat, DataType threshold, int maxEls = -1) {
		 std::vector<std::tuple<int, int, double>> ret;
		 using SMat = Eigen::SparseMatrix<DataType, Opts>;
		 for (int k = 0; k < mat.outerSize(); ++k) {
			 for (typename SMat::InnerIterator it(mat, k); it; ++it) {
				 if (std::abs(it.value()) > threshold) ret.emplace_back(it.row(), it.col(), it.value());
				 if (maxEls != -1 && ret.size() == maxEls) return ret;
			 }
		 }
		 return ret;
	 }

	/**
	  * \brief Finds elements in two matrices that violate the maximum absolute difference per cell in the matrix.
	  * \tparam DataType Datatype of sparse matrix
	  * \tparam Opts Options for sparse matrix
	  * \param mat First matrix
	  * \param mat2 Second matrix
	  * \param threshold Max abs difference thresholds
	  * \param maxEls Maximum number of elements to report. -1 for all elements
	  * \return List of elements that violate the abs difference, given as a tuple with row, columne, first mat value, second mat value.
	  */
	 template<typename DataType, int Opts>
	 inline std::vector<std::tuple<int, int, double,double>> getAbsDifference(const Eigen::SparseMatrix<DataType, Opts>& mat, 
		 const Eigen::SparseMatrix<DataType, Opts>& mat2, DataType threshold, int maxEls = -1) 
	{
		 std::vector<std::tuple<int, int, double,double>> ret;
		 int cnt = 0;
		 using SMat = Eigen::SparseMatrix<DataType, Opts>;
		 for (int k = 0; k < mat.outerSize(); ++k) {
			 using IIt = typename SMat::InnerIterator;
			 IIt it1(mat, k);
			 IIt it2(mat2, k);
			 for (typename SMat::InnerIterator it(mat, k); it; ++it) {
				 if (std::abs(it.value() - mat2.coeff(it.row(), it.col())) > threshold) {
					 cnt++;
					 if(maxEls == -1 || ret.size() != maxEls)
					 {
						 ret.emplace_back(it.row(), it.col(), it.value(), mat2.coeff(it.row(), it.col()));
					 }
				 }
			 }
		 }
		 //std::cout << "Failures: " << cnt << "/" << mat.nonZeros();
		 return ret;
	 }

	inline std::string tripletsToString(const std::vector<std::tuple<int,int,double>>& trips)
	{
		std::stringstream str;
		for(auto el : trips)
		{
			str << "[" << std::get<0>(el) << "," << std::get<1>(el) << "]: " << std::get<2>(el) << ",";
		}
		return str.str();
	}
	inline std::string tripletsToString(const std::vector<std::tuple<int, int, double,double>>& trips)
	{
		std::stringstream str;
		for (auto el : trips)
		{
			str << "[" << std::get<0>(el) << "," << std::get<1>(el) << "]: " << std::get<2>(el) << " vs " << std::get<3>(el) << ",";
		}
		return str.str();
	}
 }

#endif