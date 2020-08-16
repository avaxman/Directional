#ifndef DIRECTIONAL_DIRECTIONALGAMMA_SUITE_H
#define DIRECTIONAL_DIRECTIONALGAMMA_SUITE_H
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <directional/block_diag.h>
#include "Gamma_suite.h"
namespace directional
{
	/**
	 * Helper struct for constructing a sparse double matrix.
	 */
	struct SparseHelper
	{
		std::vector<Eigen::Triplet<double>> trips;
		int rows, cols;
		SparseHelper(int rows, int cols, int expectedCoeffs):rows(rows),cols(cols)
		{
			trips.reserve(expectedCoeffs);
		}
		void addCoeff(int r, int c, double coeff)
		{
			trips.emplace_back(r, c, coeff);
		}
		Eigen::SparseMatrix<double> toMat()
		{
			Eigen::SparseMatrix<double> m{ rows,cols };
			m.setFromTriplets(trips.begin(), trips.end());
			return m;
		}
	};

	/**
	 * TODO
	 */
	void DirectionalGamma_Suite(
		int N,
		const Eigen::MatrixXi& matching,
		Eigen::SparseMatrix<double>& C
	)
	{

	}

    inline void columndirectional_to_gamma2_matrix(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXi& E, const Eigen::MatrixXi&SFE, const Eigen::MatrixXi&EF, int N,
        Eigen::SparseMatrix<double>& projector)
	{
        Eigen::SparseMatrix<double> P;
        directional::Gamma2_projector(V, F, E, SFE, EF, P);
        std::vector< Eigen::SparseMatrix<double>*> parts(N, &P);
        directional::block_diag(parts, projector);
	}

	/**
	 * Converts Gamma2 elements to Gamma3 elements for a degree N directional field.
	 * Input:
	 * - sFE: |F| x 6 matrix containing face to edge connectivity with orientation information
	 * - N: degree of the directional field
	 * Output:
	 * - G2ToG3_N: 3 N |F| x 2 N |F| matrix converting the gamma2 directional field to a gamma3 directional field.
	 */
	void Matched_G2_To_G3(
		const Eigen::MatrixXi& sFE,
		int N,
		Eigen::SparseMatrix<double>& G2ToG3_N
	)
	{
		Eigen::SparseMatrix<double> G2ToG3;
		directional::Gamma2_To_Gamma3(sFE, G2ToG3);

		// Convert gamma2 to gamma3, then apply operator. Note that gamma2 to gamma3 does not require matching.
		std::vector<Eigen::SparseMatrix<double>*> refs;
		for (int i = 0; i < N; i++) refs.push_back(&G2ToG3);
		directional::block_diag(refs, G2ToG3_N);
	}
	/**
	 * Converts Gamma3 elements to Gamma2 elements for a degree N directional field.
	 * Input:
	 * - sFE: |F| x 6 matrix containing face to edge connectivity with orientation information
	 * - N: degree of the directional field
	 * Output:
	 * - G3ToG2_N: 2 N |F| x 3 N |F| matrix converting the gamma2 directional field to a gamma3 directional field.
	 */
	void Matched_G3_To_G2(
		const Eigen::MatrixXi& sFE,
		int N,
		Eigen::SparseMatrix<double>& G2ToG3_N
	)
	{
		Eigen::SparseMatrix<double> G3ToG2;
		directional::Gamma3_To_Gamma2(sFE, G3ToG2);

		// Convert gamma2 to gamma3, then apply operator. Note that gamma2 to gamma3 does not require matching.
		std::vector<Eigen::SparseMatrix<double>*> refs;
		for (int i = 0; i < N; i++) refs.push_back(&G3ToG2);
		directional::block_diag(refs, G2ToG3_N);
	}

	/**
	 * Computes the curl for directional Gamma2 elements.
	 * Input:
	 *  - EF: |E| x 2 matrix. See EdgeData struct
	 *  - sFE: |F| x 6 matrix. See EdgeData struct
	 *  - EI: |E| x 2 matrix. See EdgeData struc
	 *  - Matching: |E| x 2 matrix. Level change when moving over edge
	 *  - faceCount : Number of aces
	 *  - N : Directional degree
	 * Output:
	 *  - Matched directional curl operator
	 */
	inline void Matched_Curl(
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EI,
		const Eigen::MatrixXi& Matching,
		int faceCount,
		int N,
		Eigen::SparseMatrix<double>& Curl
	)
	{
		std::cout << "Building matched curl" << std::endl;
		std::vector<Eigen::Triplet<double>> trips;
		trips.reserve(2 * N * EF.rows());
		for (int e = 0; e < EF.rows(); e++)
		{
			// Gamma3 IDs
			const int lGamma = 3 * EF(e, 0) + EI(e, 0);
			const int rGamma = 3 * EF(e, 1) + EI(e, 1);
			for (int n = 0; n < N; n++)
			{
				trips.emplace_back(e + n * EF.rows(), lGamma + n * faceCount * 3, -1);

				//Right gamma has to be compensated for.
				const int matchedLevel = (Matching(e, 0) + n) % N;
				trips.emplace_back(e + n * EF.rows(), rGamma + matchedLevel * faceCount * 3, 1);
			}
		}
		Eigen::SparseMatrix<double> G2ToG3;
		directional::Matched_G2_To_G3(sFE, N, G2ToG3);
		Curl = Eigen::SparseMatrix<double>(N * EF.rows(), N * 3 * sFE.rows());
		Curl.setFromTriplets(trips.begin(), trips.end());
		Curl = Curl * G2ToG3;
		std::cout << "Built matched curl" << std::endl;
	}

	void Matched_Curl_To_Gamma3(
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EI,
		const Eigen::MatrixXi& Matching,
		int N,
		Eigen::SparseMatrix<double>& Curl_To_Gamma3
	)
	{
		const int faceCount = sFE.rows();
		const int gammaCount = 3 * faceCount;
		const int edgeCount = EF.rows();
		SparseHelper sh(3 * faceCount * N, edgeCount * N, 2 * edgeCount * N);
		for(int e = 0; e < edgeCount; e++)
		{
			const int lGamma = 3 * EF(e, 0) + EI(e, 0);
			const int rGamma = 3 * EF(e, 1) + EI(e, 1);
			for(int n = 0; n < N; n++)
			{
				sh.addCoeff(lGamma + n * gammaCount, e + n * edgeCount, -1);
				const int level = (n + Matching(e, 0)) % N;
				sh.addCoeff(rGamma + level * gammaCount, e + n * edgeCount, 1);
			}
		}
		Curl_To_Gamma3 = sh.toMat();
	}
	void Matched_Curl_To_Gamma2(
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EI,
		const Eigen::MatrixXi& Matching,
		int N,
		Eigen::SparseMatrix<double>& Curl_To_Gamma2
	)
	{
		Matched_Curl_To_Gamma3(EF, sFE, EI, Matching, N, Curl_To_Gamma2);
		// Map gamma2 to gamma3 before applying operator.
		Eigen::SparseMatrix<double> G3ToG2_N;
		Matched_G3_To_G2(sFE, N, G3ToG2_N);
		Curl_To_Gamma2 = G3ToG2_N * Curl_To_Gamma2;
	}

	inline void Matched_A_To_Gamma3(
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EI,
		const Eigen::MatrixXi& Matching,
		int faceCount,
		int N,
		Eigen::SparseMatrix<double>& A_To_Gamma3
	)
	{
		const int edgeCount = EF.rows();
		const int gammaCount = 3 * faceCount;
		SparseHelper sh(N * gammaCount, edgeCount * N, N * edgeCount * 2);
		for(int f = 0; f < sFE.rows(); f++)
		{
			for(int j = 0; j < 3; j++)
			{
				const int e = sFE(f, j);
				const int s = sFE(f, j+3);
				const int gam = 3 * f + j;
				if(s == 0)
				{
					for(int n = 0; n < N; n++)
					{
						sh.addCoeff(gam + n * gammaCount, e + n * edgeCount, 1);
					}
				}
				else
				{
					for (int n = 0; n < N; n++)
					{
						const int level = (n+Matching(e, 0)) % N;
						sh.addCoeff(gam + level * gammaCount, e + n * edgeCount, 1);
					}
				}
			}
		}
		std::cout << "Matched_A_To_Gamma3 built" << std::endl;
		A_To_Gamma3 = sh.toMat();
	}


	void Matched_A_To_Gamma2(
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EI,
		const Eigen::MatrixXi& Matching,
		int faceCount,
		int N,
		Eigen::SparseMatrix<double>& A_To_Gamma2
	)
	{
		Matched_A_To_Gamma3(EF, sFE, EI, Matching, faceCount, N, A_To_Gamma2);
		Eigen::SparseMatrix<double> G3ToG2_N;
		directional::Matched_G3_To_G2(sFE,N, G3ToG2_N);

		// Convert gamma2 to gamma3, then apply operator. Note that gamma2 to gamma3 does not require matching.
		A_To_Gamma2 = G3ToG2_N * A_To_Gamma2;
	}

    /**
	 * \brief Trivially copies each edge field value to the left and right gamma nex to it. Follows the matching for selecting appropriate levels.
	 * \param EI 
	 * \param EF 
	 * \param sFE 
	 * \param Matching 
	 * \param faceCount 
	 * \param N 
	 * \param Gamma2_To_A 
	 */
	inline void Matched_Gamma2_To_E(
		const Eigen::MatrixXi& EI,
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& Matching,
		int faceCount,
		int N,
		Eigen::SparseMatrix<double>& Gamma2_To_A
	)
	{
		const int edgeCount = EI.rows();
		const int gammaCount = 3 * faceCount;
		SparseHelper sh(N * edgeCount, N * 3 * faceCount, 2 * 3 * N * faceCount);
		for(int e = 0; e < EI.rows(); e++)
		{
			const int gL = 3 * EF(e, 0) + EI(e, 0);
			const int gR = 3 * EF(e, 1) + EI(e, 1);
			for(int n = 0; n < N; n++)
			{
				sh.addCoeff(e + n * edgeCount, gL + n * gammaCount, 1);
				const int rLevel = (n + Matching(e, 0)) % N;
				sh.addCoeff(e + n  * edgeCount, gL + rLevel * gammaCount, 1);
			}
		}
		Eigen::SparseMatrix<double> G2ToG3_N;
		Matched_G2_To_G3(sFE, N, G2ToG3_N);
		Gamma2_To_A = sh.toMat() * G2ToG3_N;
	}

	void Matched_A_C_To_F(
		const Eigen::MatrixXi& sFE,
		int edgeCount,
		int N,
		const Eigen::MatrixXi& Matching,
		Eigen::SparseMatrix<double>& A_C_To_F
	)
	{
		const int faceCount = sFE.rows();
		SparseHelper sh(N * sFE.rows(), edgeCount, 3 * sFE.rows() * N);
		for(int f = 0; f < sFE.rows(); f++)
		{
			for(int j = 0; j < 3; j++)
			{
				const int e = sFE(f, j);
				const int s = sFE(f, j + 3);
				//CCW oriented
				if(s == 0)
				{
					for(int n = 0; n < N; n++)
					{
						// Level of curl corresponds to that of the face.
						sh.addCoeff(n * faceCount + f, e + n * edgeCount, 1);
					}
				}
				// CW oriented: matching is needed
				else
				{
					for (int n = 0; n < N; n++)
					{
						const int level = (n + N - Matching(e, 1)) % N;
						// Level of curl corresponds to that of the face.
						sh.addCoeff(n * faceCount + f, e + level * edgeCount, 1);
					}
				}
			}
		}
		A_C_To_F = sh.toMat();
	}

	void Matched_D1(
		const Eigen::MatrixXi& sFE,
		int edgeCount,
		int N,
		const Eigen::MatrixXi& Matching,
		Eigen::SparseMatrix<double>& D1
	)
	{
		const int faceCount = sFE.rows();
		SparseHelper sh(N * sFE.rows(), edgeCount, 3 * sFE.rows() * N);
		for (int f = 0; f < sFE.rows(); f++)
		{
			for (int j = 0; j < 3; j++)
			{
				const int e = sFE(f, j);
				const int s = sFE(f, j + 3);
				//CCW oriented
				if (s == 0)
				{
					for (int n = 0; n < N; n++)
					{
						// Level of curl corresponds to that of the face.
						sh.addCoeff(n * faceCount + f, e + n * edgeCount, 1);
					}
				}
				// CW oriented: matching is needed
				else
				{
					for (int n = 0; n < N; n++)
					{
						const int level = (n + N - Matching(e, 1)) % N;
						// Level of curl corresponds to that of the face.
						sh.addCoeff(n * faceCount + f, e + level * edgeCount, -1);
					}
				}
			}
		}
		D1 = sh.toMat();
	}


	void Matched_Gamma2_To_AC(
		const Eigen::MatrixXi& EI,
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& Matching,
		int N,
		Eigen::SparseMatrix<double>& G2ToAC
	)
	{
		const int faceCount = sFE.rows();
		Eigen::SparseMatrix<double> Gamma2_To_A, Curl;
		Matched_Gamma2_To_E(EI, EF, sFE, Matching, faceCount, N, Gamma2_To_A);
		Matched_Curl(EF, sFE, EI, Matching, faceCount, N, Curl);
		std::cout << "G2 To AC ops sizes: " <<
			Gamma2_To_A.rows() << "," << Gamma2_To_A.cols() << " && " <<
			Curl.rows() << "," << Curl.cols() << std::endl;
		G2ToAC = 0.5* igl::cat(1, Gamma2_To_A, Curl);
	}

	void Matched_AC_To_Gamma3(
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EI,
		const Eigen::MatrixXi& Matching,
		int N,
		Eigen::SparseMatrix<double>& AC_To_Gamma3
	)
	{
		const int faceCount = sFE.rows();
		Eigen::SparseMatrix<double> A_To_Gamma3, C_To_Gamma3;
		Matched_A_To_Gamma3(EF, sFE, EI, Matching, faceCount, N, A_To_Gamma3);
		Matched_Curl_To_Gamma3(EF, sFE, EI, Matching, N, C_To_Gamma3);
		AC_To_Gamma3 = igl::cat(2, A_To_Gamma3, C_To_Gamma3);
	}
	void Matched_AC_To_Gamma2(
		const Eigen::MatrixXi& EF,
		const Eigen::MatrixXi& sFE,
		const Eigen::MatrixXi& EI,
		const Eigen::MatrixXi& Matching,
		int N,
		Eigen::SparseMatrix<double>& AC_To_Gamma2
	)
	{
		const int faceCount = sFE.rows();
		Eigen::SparseMatrix<double> A_To_Gamma2, C_To_Gamma2;
		std::cout << "Building A to G2" << std::endl;
		Matched_A_To_Gamma2(EF, sFE, EI, Matching, faceCount, N, A_To_Gamma2);
		std::cout << "Building C to G2" << std::endl;
		Matched_Curl_To_Gamma2(EF, sFE, EI, Matching, N, C_To_Gamma2);
		std::cout << "Combining, sizes: " << 
			C_To_Gamma2.rows() << "," << C_To_Gamma2.cols() << " && " <<
			A_To_Gamma2.rows() << "," << A_To_Gamma2.cols() <<
			std::endl;
		AC_To_Gamma2 = igl::cat(2, A_To_Gamma2, C_To_Gamma2);
		std::cout << "Done" << std::endl;
	}
}
#endif