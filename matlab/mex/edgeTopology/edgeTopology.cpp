/* system header */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef VERBOSE
#define dir_LOGLN(expr) 
#else
#define dir_LOGLN(expr) std::cout << expr << std::endl ;
#endif

/* MEX header */
#include <mex.h> 
#include "matrix.h"
#include "../mexHelpers.h"
#include "../FuncArgs.h"
/* Directional/IGL headers*/
#include <directional/Subdivision/subdivision.h>
#include <directional/curl_matching.h>
#include <directional/parameterize.h>
#include <directional/effort_to_indices.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/combing.h>
#include <directional/polycurl_reduction.h>

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

enum Input
{
	In_F = 0,
	In_V,
	In_Count
};
enum Output
{
	Out_ED =0,
	Out_Count
};
#define ID_PREFIX "MATLAB:curlFreeParametrization:"

void mainFunc(mexhelpers::FuncArgs& args)
{
	using namespace mexhelpers;
	
    /* argument check */
	args.verifyInputCount(In_Count, ID_PREFIX "inputmismatch", "Expected 2 input argument.");
	args.verifyOutputCount(Out_Count, ID_PREFIX "outputmismatch", "Expected 1 output argument.");

	// Get inputs
	Eigen::MatrixXi F = args.getIndexMat(Input::In_F);
	Eigen::MatrixXd V = args.getMat(In_V);

	args.failIf(F.cols() != 3, "Expected F with 3 cols");
	args.failIf(V.cols() != 3, "Expected V with 3 cols");
	// Setup topology
	EdgeData ed;
	ed.construct(F);



	// Setup output
	args.setOutput(Out_ED, EDToStruct(ed));
}
/* Entry for building subdivision operators
 * Function signature:
 *  [ED0_Struct, EDK_Struct, S_Struct] = buildSubdivisionOperators(F, level)
 *  Here the ED structs have keys EF< sFE, EI and F, corresponding to the topology data.
 *  The S struct contains the subdivision operators under the keys F, C, E, V and Gamma
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	mexhelpers::FuncArgs args(plhs, nlhs, prhs, nrhs);
	try
	{
		mainFunc(args);
	}
	catch(const std::exception& e)
	{
		mexErrMsgTxt(e.what());
	}
}