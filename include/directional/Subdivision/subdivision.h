#ifndef DIRECTIONAL_SUBDIVISION_H
#define DIRECTIONAL_SUBDIVISION_H
#include <directional/Subdivision/SubdivisionBuilder.h>

// Subdividers
#include <directional/Subdivision/VertexFieldSubdivider.h>
#include <directional/Subdivision/OneFormSubdivider.h>
#include <directional/Subdivision/FaceFieldSubdivider.h>
#include <directional/Subdivision/EdgeFieldSubdivider.h>

// Coefficient providers
#include <directional/Subdivision/LoopCoefficientProvider.h>
#include <directional/Subdivision/HalfBoxSplineCoefficientProvider.h>
#include <directional/Subdivision/HalfCurlCoefficientProvider.h>
#include <directional/Subdivision/OneFormCoefficientProvider.h>


#include <Eigen/Sparse>
#include <Eigen/Eigen>

namespace directional
{
	namespace subdivision
	{
		using Loop_S_v = VertexFieldSubdivider<LoopCoefficientProvider>;
		using HB_S_f = FaceFieldSubdivider<HalfboxSplineCoefficientProvider>;
		using DirFS_S_c = EdgeFieldSubdivider<HalfCurlCoefficientProvider>;
		using HDirFS_S_e = OneFormSubdivider<OneFormCoefficientProvider>;
	}
}

void subdivision_operators(const Eigen::MatrixXi& F, const Eigen::MatrixXd& V, int level,
Eigen::SparseMatrix<double>& S_V,
Eigen::SparseMatrix<double>& S_F,
Eigen::SparseMatrix<double>& S_E,
Eigen::SparseMatrix<double>& S_C,
EdgeData& ED0,
EdgeData& EDL)
{
	std::cout << "### Subdivision construction" << std::endl;
	SimpleTimer st;
	st.start();
	
	// Subdividers for different field types
	VertexFieldSubdivider<LoopCoefficientProvider> vSub;
	EdgeFieldSubdivider<HalfCurlCoefficientProvider> cSub;
	FaceFieldSubdivider<HalfboxSplineCoefficientProvider> fSub;
	OneFormSubdivider<OneFormCoefficientProvider> eSub;

	// Create the subdivision builder
	auto builder = create_sub_builder(F, vSub, cSub, eSub, fSub);
	// Record initial edge topology
	ED0 = builder.ED;
	st.stop();
	std::cout << "ED construction " << st.elapsed() << " ms" << std::endl;
	assert(ED0.isConsistent());
	std::cout << "Mesh: |F| = " << F.rows() << " , |V| = " << ED0.vertexCountFast() << ", |E| = " << ED0.E.rows() << std::endl;

	// Create the subdivision matrices incrementally.
	builder.constructUpToLevel(level);

	st.reset().start();
	// Retrieve subdivision matrices
	S_V = vSub.getMatrix();
	S_F = fSub.getMatrix();
	S_E = eSub.getMatrix();
	S_C = cSub.getMatrix();
	st.stop();
	std::cout << "Getting matrices in " << st.elapsed() << " ms" << std::endl;

	// Save the fine level topology
	EDL = builder.ED;
}
void subdivision_operators(const Eigen::MatrixXi& F, const Eigen::MatrixXd& V, int level,
	Eigen::SparseMatrix<double>& S_V,
	Eigen::SparseMatrix<double>& S_E,
	EdgeData& ED0,
	EdgeData& EDL)
{
	SimpleTimer st;
	st.start();


	assert(ED0.isConsistent());
	
	VertexFieldSubdivider<LoopCoefficientProvider> vSub;
	OneFormSubdivider<OneFormCoefficientProvider> eSub;

	auto builder = create_sub_builder(F, vSub, eSub);
	ED0 = builder.ED;
	std::cout << "Mesh: |F| = " << F.rows() << " , |V| = " << V.rows() << ", |E| = " << ED0.E.rows() << std::endl;
	// Create the subdivision matrices incrementally.
	builder.constructUpToLevel(level);

	st.reset().start();
	S_V = vSub.getMatrix();
	S_E = eSub.getMatrix();
	st.stop();
	std::cout << "Getting matrices in " << st.elapsed() << " ms" << std::endl;
	EDL = builder.ED;
}
#endif