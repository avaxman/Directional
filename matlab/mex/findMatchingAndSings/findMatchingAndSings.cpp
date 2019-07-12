#include <directional/DirectionalGamma_Suite.h>
//#include <directional/Subdivision/DirectionalSubdivider.h>
#include <directional/Subdivision/directional_subdivision.h>

#include <mex.h>
#include <mat.h>
#include "../mexHelpers.h"
#include "../FuncArgs.h"
#include <directional/effort_to_indices.h>
#include <directional/principal_matching.h>

enum Input
{
	In_F = 0,
	In_V,
	In_Rawfield,
	In_Count
};
enum Output
{
	Out_Matching = 0,
	Out_SingularityVertices,
	Out_SingularityIndices,
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

	Eigen::MatrixXd V = args.getMat(In_V);
	Eigen::MatrixXi F = args.getIndexMat(In_F);
	Eigen::MatrixXd rawField = args.getMat(In_Rawfield);

	EdgeData ED(F);

	Eigen::VectorXi matching;
	Eigen::VectorXd effort;
	Eigen::VectorXi singVertices;
	Eigen::VectorXi singIndices;
	const int N = rawField.cols() / 3;
	std::cout << "Staring matching " << std::endl;
	//combing and cutting
	directional::principal_matching(V, F, ED.E, ED.EF, ED.sFE, rawField, matching, effort);
	directional::effort_to_indices(V, F, ED.E, ED.EF, effort, matching, N, singVertices, singIndices);

	args.setOutput(Out_Matching, matching);
	args.setOutput(Out_SingularityVertices, mexhelpers::Indexed(singVertices));
	std::cout << "Sings" << singVertices << std::endl;
	args.setOutput(Out_SingularityIndices, singIndices);
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