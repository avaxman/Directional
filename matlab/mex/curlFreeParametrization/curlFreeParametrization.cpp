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

void to_sparse_adjacency(const Eigen::MatrixXi& adjLists, int srcSize, Eigen::SparseMatrix<double>& output)
{
	std::vector<Eigen::Triplet<double>> trips;
	output = Eigen::SparseMatrix<double>(adjLists.rows(), srcSize);
	trips.reserve(adjLists.rows() * adjLists.cols());
	for(int i = 0; i < adjLists.rows(); i++)
	{
		for(int j =0; j < adjLists.cols(); j++)
		{
			trips.emplace_back(i, adjLists(i, j), 1);
		}
	}
	output.setFromTriplets(trips.begin(), trips.end());
}
void to_sparse_directed_adjacency(const Eigen::MatrixXi& adjLists, int srcSize, Eigen::SparseMatrix<double>& output)
{
	assert(adjLists.cols() == 2);
	std::vector<Eigen::Triplet<double>> trips;
	output = Eigen::SparseMatrix<double>(adjLists.rows(), srcSize);
	trips.reserve(adjLists.rows() * adjLists.cols());
	for (int i = 0; i < adjLists.rows(); i++)
	{
		trips.emplace_back(i, adjLists(i,0), -1);
		trips.emplace_back(i, adjLists(i, 1), -1);
	}
	output.setFromTriplets(trips.begin(), trips.end());
}

enum Input
{
	In_F = 0,
	In_V,
	In_RawField,
	In_Count
};
enum Output
{
	Out_CutV =0,
	Out_CutF = 1,
	Out_UV,
	Out_CombedField,
	Out_Matching,
	Out_EF,
	Out_SingVerts,
	Out_SingInds,
	Out_PoissonError,
	Out_Count
};
#define ID_PREFIX "MATLAB:curlFreeParametrization:"

void mainFunc(mexhelpers::FuncArgs& args)
{
	using namespace mexhelpers;
	
    /* argument check */
	args.verifyInputCount(In_Count);
	args.verifyOutputCount(Out_Count);

	// Get inputs
	Eigen::MatrixXi F = args.getIndexMat(Input::In_F);
	Eigen::MatrixXd V = args.getMat(In_V);
	//Eigen::VectorXi singVertices = args.getIndexVec(Input::In_SingularVertices);
	//Eigen::VectorXi matching = args.getIntegerVec(Input::In_Matching);
	Eigen::MatrixXd rawField = args.getMat(In_RawField);
	int N = rawField.cols() / 3;

	args.failIf(rawField.cols() % 3 != 0, "Expected rawfield as N directionals per face in a row");
	args.failIf(F.cols() != 3, "Expected F with 3 cols");
	args.failIf(V.cols() != 3, "Expected V with 3 cols");
	//args.failIf(singVertices.cols() != 1, "Expected SingVertices with 1 col");
	//args.failIf(singVertices.minCoeff() < 0 || singVertices.maxCoeff() >= V.rows(), "Singularity indices out of bounds");
	//args.failIf(matching.cols() != 1, "Singularity indices out of bounds");
	// Setup topology
	Eigen::MatrixXi EV, FE, EF;
	igl::edge_topology(V, F, EV, FE, EF);
	args.failIf(F.rows() != rawField.rows(), "Expected F and rawfield with same number of rows");
	//EdgeData ed;
	//ed.construct(F);
	Eigen::VectorXi matching, sings, singInds;
	Eigen::VectorXd effortInit, curlNormInit;
	std::cout << "Redo curl matching" << std::endl;
	directional::curl_matching(V, F, EV, EF, FE, rawField, matching, effortInit, curlNormInit);
	std::cout << "Cal indices" << std::endl;
	directional::effort_to_indices(V, F, EV, EF, effortInit, matching, N, sings, singInds);

	// Start parameterization
	directional::ParameterizationData pd;
	std::cout << "Cutting mesh" << std::endl;
	directional::cut_mesh_with_singularities(V, F, sings, pd.face2cut);
	Eigen::VectorXi combedMatching;
	Eigen::VectorXd combedEffort;

	Eigen::MatrixXd combedField;

	std::cout << "Setting up matching" << std::endl;

	//Eigen::MatrixXi FE;
	//ed.getEdgeTopologyCompatibleFE(FE);

	/*directional::effort_to_indices(V, F, ed.E, ed.EF, effortCF, matchingCF, N, singVerticesCF, singIndicesCF);
	directional::combing(V, F, ed.E, ed.EF, FE, rawField, matchingCF, combedFieldCF);
	directional::curl_matching(V, F, ed.E, ed.EF, FE, combedFieldCF, combedMatchingCF, combedEffortCF, combedCurlCF);*/

	directional::combing(V, F, EV, EF, FE, pd.face2cut, rawField, matching, combedField);
	std::cout << "Performing principal matching" << std::endl;
	Eigen::VectorXd curlNorm;
	directional::curl_matching(V, F, EV, EF, FE, combedField, combedMatching, combedEffort, curlNorm);
	Eigen::VectorXi singsCF, singIndsCF;
	directional::effort_to_indices(V, F, EV, EF, combedEffort, combedMatching, N, singsCF, singIndsCF);
	//directional::principal_matching(V, F, ed.E, ed.EF, FE, combedField, combedMatching, combedEffort);

	std::cout << "Setting up parameterization" << std::endl;

	Eigen::MatrixXd CutV, cutUV;
	Eigen::MatrixXi CutF;
	directional::setup_parameterization(N, V, F, EV, EF, FE, combedMatching, singsCF, pd, CutV, CutF);

	double lengthRatio = 0.01;
	bool isInteger = false;  //do not do translational seamless.
	std::cout << "Solving parameterization" << std::endl;
	directional::parameterize(V, F, FE, combedField, lengthRatio, pd, CutV, CutF, isInteger, cutUV);
	std::cout << "Done!" << std::endl;

	// Setup output
	args.setOutput(Out_CutF, Indexed(CutF));
	args.setOutput(Out_CutV, CutV);
	args.setOutput(Out_UV, cutUV);
	args.setOutput(Out_CombedField, combedField);
	args.setOutput(Out_Matching, combedMatching);
	args.setOutput(Out_SingInds, singIndsCF);
	args.setOutput(Out_SingVerts, Indexed(singsCF));
	args.setOutput(Out_EF, Indexed(EF));
	args.setOutput(Out_PoissonError, pd.poissonError);
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