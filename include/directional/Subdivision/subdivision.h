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

void subdivision_operators(const Eigen::MatrixXi& F, const Eigen::MatrixXd& V, int level,
Eigen::SparseMatrix<double>& S_V,
Eigen::SparseMatrix<double>& S_F,
Eigen::SparseMatrix<double>& S_E,
Eigen::SparseMatrix<double>& S_C,
EdgeData& ED0,
EdgeData& EDL)
{
	SimpleTimer st;
	st.start();
	

	ED0.construct(F);
	st.stop();
	std::cout << "ED construction " << st.elapsed() << " ms" << std::endl;
	assert(ED0.isConsistent());
	std::cout << "Mesh: |F| = " << F.rows() << " , |V| = " << V.rows() << ", |E| = " << ED0.E.rows() << std::endl;
	
	VertexFieldSubdivider<LoopCoefficientProvider> vSub;
	EdgeFieldSubdivider<HalfCurlCoefficientProvider> cSub;
	FaceFieldSubdivider<HalfboxSplineCoefficientProvider> fSub;
	OneFormSubdivider<OneFormCoefficientProvider> eSub;

	auto builder = create_sub_builder(F, vSub, cSub, eSub, fSub);
	// Create the subdivision matrices incrementally.
	builder.constructUpToLevel(level);

	st.reset().start();
	S_V = vSub.getMatrix();
	S_F = fSub.getMatrix();
	S_E = eSub.getMatrix();
	S_C = cSub.getMatrix();
	st.stop();
	std::cout << "Getting matrices in " << st.elapsed() << " ms" << std::endl;
}