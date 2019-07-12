#include <directional/DirectionalGamma_Suite.h>
//#include <directional/Subdivision/DirectionalSubdivider.h>
#include <directional/Subdivision/directional_subdivision.h>

#include <mex.h>
#include <mat.h>
#include "../mexHelpers.h"
#include "../FuncArgs.h"

struct IndicesConstructor
{
	Eigen::VectorXi& out;
	Eigen::VectorXi& loopOut;
	const EdgeData& ED;
	const Eigen::MatrixXi& Matching;
	int N;
	IndicesConstructor(const EdgeData& ED, const Eigen::MatrixXi& Matching, int N, int vCount, Eigen::VectorXi& out, Eigen::VectorXi& loopOut):
	out(out),ED(ED),Matching(Matching),
	N(N),loopOut(loopOut)
	{
		this->out.setConstant(vCount, 0);
		loopOut.setConstant(vCount, 0);
	}

	void handleRegularRing(const std::vector<int>& edges, const std::vector<int>& edgeSides)
	{
		const int centralV = ED.E(edges[0], edgeSides[0]);
		int level = 0;
		for(int i = 1; i < edges.size() / 2; i++)
		{
			const int e = 2 * i;
			level += Matching(edges[e], 1 - edgeSides[e]);
		}
		level += Matching(edges[0], 1 - edgeSides[0]);
		level = level % N;
		out(centralV) = level;

		Eigen::MatrixXi faceLevels(N, edges.size() / 2);
		faceLevels.col(0) = Eigen::VectorXi::LinSpaced(N, 0, N - 1);
		for (int i = 1; i < edges.size()/2; i++) {
			const int e = 2 * i;
			// Update new levels for target face by applying matching
			faceLevels.col(i) = circShift(faceLevels.col(i - 1), -Matching(edges[e], 1 - edgeSides[e]));
			// Next face. Same face will be seen multiple times due to unwrapping.
		}
		Eigen::VectorXi out = circShift(faceLevels.col(edges.size()/2 - 1), -Matching(edges[0], 1 - edgeSides[0]));
		loopOut(centralV) = out(0);
	}
	void handleBoundaryRing(const std::vector<int>& edges, const std::vector<int>& edgeSides)
	{

	}
};

void indicesFromMatching(EdgeData& ED, const Eigen::MatrixXi& matching, int N, int vCount, Eigen::VectorXi& sings, Eigen::VectorXi& singsLooped)
{
	IndicesConstructor ic(ED,matching, N, vCount,sings,singsLooped);
	DCEL dcel(ED);
	dcel.iterateRings(ic);
}

enum Input
{
	In_F,
	In_SubdivLevel,
	In_Matching,
	In_N,
	In_SingularityMat,
	In_Count
};
enum Output
{
	Out_ED0,
	Out_EDK,
	Out_SubdivisionStruct,
	Out_Matching,
	Out_Count	
};
#define ID_PREFIX "MATLAB:buildDirectionalSubdivisionOperators:"

mexhelpers::MexStructure EDToStruct(const EdgeData& ED)
{
	const char* fields[] = { "EF","sFE","EI","F","E" };
	mexhelpers::MexStructure str(fields, 5);
	str.setIndexValue("EF", ED.EF);
	str.setIndexValue("sFE", ED.sFE);
	str.setIndexValue("EI", ED.EI);
	str.setIndexValue("F", ED.F);
	str.setIndexValue("E", ED.E);
	dir_LOGLN("Added edge data struct");
	return str;
}

void mainFunc(mexhelpers::FuncArgs& args)
{
	// Verify number of inputs given
	args.verifyInputCount(In_Count, ID_PREFIX "inputmismatch", "Expected different number of input");
	args.verifyOutputCount(Out_Count, ID_PREFIX "outputmismatch", "Expected different number of outputs");
	
	// Get inputs
	Eigen::MatrixXi F = args.getIndexMat(In_F);
	Eigen::MatrixXi matching = args.getIntegerMat(In_Matching);
	int subdivisionLevel = args.getInt(In_SubdivLevel);
	int N = args.getInt(In_N);
	// Matrix containing 0/1 in column 0 signalling a singularity, and the value N * fracIndex in column 1 if it is a singularity.
	Eigen::MatrixXi singularityMatrix = args.getIntegerMat(In_SingularityMat);
	
	// Extra checking
	args.failIf(F.cols() != 3, "Expected F with 3 columns");
	args.failIf(matching.cols() != 2, "Expected matching with 2 columns");
	args.failIf(singularityMatrix.cols() != 2, "Expected singularityMatrix with 2 columns");

	// Output data
	Eigen::SparseMatrix<double> S_V, S_DE, S_DC, S_DGamma;
	EdgeData ED0(F), EDL;
	Eigen::MatrixXi matchingK;

	Eigen::VectorXi inds, indsLooped;
	indicesFromMatching(ED0, matching, N, ED0.vertexCount(), inds, indsLooped);
	
	for(int i = 0; i < inds.rows(); i++)
	{
		if(inds(i) != indsLooped(i))
		{
			std::stringstream sstr;
			sstr << "Failed singularity check: found " << inds(i) << ", but got from looped" << indsLooped(i) << std::endl;
			mexErrMsgTxt(sstr.str().c_str());
		}
		if (inds(i) == 0 && (singularityMatrix(i, 0) != 0 || singularityMatrix(i, 1) != 0)) {
			std::stringstream sstr;
			sstr << "Failed singularity check: found no, but got as input " << singularityMatrix(i, 0) << std::endl;
			mexErrMsgTxt(sstr.str().c_str());
		}
		if(inds(i) != 0 && (singularityMatrix(i,0) == 0 || singularityMatrix(i,1) != inds(i)))
		{
			std::stringstream sstr;
			sstr << "Failed singularity check: found ind " << inds(i) << ", but got as input " << singularityMatrix(i, 0) << std::endl;
			mexErrMsgTxt(sstr.str().c_str());
		}
	}

	// Construct subdivision operators
	directional_subdivision_operators(F, subdivisionLevel, N, matching, singularityMatrix,
		S_V, S_DE, S_DC, matchingK, ED0, EDL);

	Eigen::SparseMatrix<double> G2ToAC_0, ACToG2_K;
	directional::Matched_Gamma2_To_AC(ED0.EI, ED0.EF, ED0.sFE, matching, N, G2ToAC_0);
	
	directional::Matched_AC_To_Gamma2(EDL.EF, EDL.sFE, EDL.EI, matchingK, N, ACToG2_K);
	
	directional::block_diag({ &S_DE, &S_DC }, S_DGamma);

	S_DGamma = ACToG2_K * S_DGamma * G2ToAC_0;
	std::cout << "Full D Gamma done, starting output" << std::endl;
	const char* fields[] = { "V","DE","DC","DGamma" };
	mexhelpers::MexStructure subdivStruct(fields, 4);
	subdivStruct.setValue("V", S_V);
	subdivStruct.setValue("DE", S_DE);
	subdivStruct.setValue("DC", S_DC);
	subdivStruct.setValue("DGamma", S_DGamma);
	std::cout << "Setting output args" << std::endl;
	args.setOutput(Out_ED0, EDToStruct(ED0));
	args.setOutput(Out_EDK, EDToStruct(EDL));
	args.setOutput(Out_SubdivisionStruct, subdivStruct);
	args.setOutput(Out_Matching, matchingK);
	std::cout << "All done" << std::endl;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mexhelpers::FuncArgs args(plhs, nlhs, prhs, nrhs);
	try
	{
		mainFunc(args);
	}
	catch (const std::exception& e)
	{
		mexErrMsgTxt(e.what());
	}
}