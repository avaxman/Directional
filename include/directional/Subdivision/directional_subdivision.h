#ifndef DIRECTIONAL_DIRECTIONAL_SUBDIVISION_H
#define DIRECTIONAL_DIRECTIONAL_SUBDIVISION_H
#include <Eigen/Eigen>
#include <directional/Subdivision/EdgeData.h>
#include <directional/Subdivision/DirectionalSubdivider.h>
#include <directional/Subdivision/VertexFieldSubdivider.h>

#include <directional/Subdivision/HalfCurlCoefficientProvider.h>
#include <directional/Subdivision/OneFormCoefficientProvider.h>
#include <directional/Subdivision/LoopCoefficientProvider.h>

void directional_subdivision_operators(
	const Eigen::MatrixXi& F, 
	int level,
	int N,
	const Eigen::MatrixXi& Matching,
	const Eigen::MatrixXi& SingularityMatrix,
	Eigen::SparseMatrix<double>& S_V,
	Eigen::SparseMatrix<double>& S_DE,
	Eigen::SparseMatrix<double>& S_DC,
	Eigen::MatrixXi& MatchingK,
	EdgeData& ED0,
	EdgeData& EDL)
{
	std::cout << "### Subdivision construction" << std::endl;
	SimpleTimer st;
	st.start();

	MatchingK = Matching;

	// Subdividers for different field types
	DirectionalSubdivider<HalfCurlCoefficientProvider> cSub(N, &SingularityMatrix, &MatchingK, false);
	DirectionalSubdivider<OneFormCoefficientProvider> eSub(N, &SingularityMatrix, &MatchingK, true);
	VertexFieldSubdivider<LoopCoefficientProvider> vSub;

	// Create the subdivision builder
	auto builder = create_sub_builder(F, vSub, cSub, eSub);
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
	S_DE = eSub.getMatrix();
	S_DC = cSub.getMatrix();
	st.stop();
	std::cout << "Getting matrices in " << st.elapsed() << " ms" << std::endl;

	// Save the fine level topology
	EDL = builder.ED;
}

#endif