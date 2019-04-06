#include "SubdivisionBuilder.h"

// Subdividers
#include "VertexFieldSubdivider.h"
#include "OneFormSubdivider.h"
#include "FaceFieldSubdivider.h"
#include "EdgeFieldSubdivider.h"

// Coefficient providers
#include "LoopCoefficientProvider.h"
#include "HalfBoxSplineCoefficientProvider.h"
#include "HalfCurlCoefficientProvider.h"
#include "OneFormCoefficientProvider.h"

#include <Eigen/Sparse>
#include <Eigen/Eigen>

void subdivision_operator(const Eigen::MatrixXi& F, const Eigen::MatrixXd& V, int level,
Eigen::SparseMatrix<double>& S_V,
Eigen::SparseMatrix<double>& S_F,
Eigen::SparseMatrix<double>& S_E,
Eigen::SparseMatrix<double>& S_C,
EdgeData& ED0,
EdgeData& EDL)
{
	ED0.construct(F);
	VertexFieldSubdivider<LoopCoefficientProvider> vSub(V.rows(), V.rows());
	EdgeFieldSubdivider<HalfCurlCoefficientProvider> cSub(ED.edgeCount(), ED.edgeCount());
	FaceFieldSubdivider<HalfboxSplineCoefficientProvider> fSub(ED.faceCount(), ED.faceCount());
	OneFormSubdivider<OneFormCoefficientProvider> eSub(ED.edgeCount(), ED.edgeCount());

	auto builder = create_sub_builder(F, vSub, cSub, eSub, fSub);
	// Create the subdivision matrices incrementally.
	builder.constructUpToLevel(level);

	vSub.
}