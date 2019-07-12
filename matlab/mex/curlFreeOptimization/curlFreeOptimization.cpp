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
	In_N,
	In_Batches,
	In_Count
};
enum Output
{
	Out_CF_Field =0,
	Out_CF_Combed_Field = 1,
	Out_CF_Matching,
	Out_CF_Combed_Matching,
	Out_SingVertices,
	Out_SingIndices,
	Out_Max_Curl,
	//Out_EF,
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
	Eigen::MatrixXd rawField = args.getMat(In_RawField);
	int N = args.getInt(In_N);
	int batches = args.getInt(In_Batches);

	args.failIf(rawField.cols() % 3 != 0, "Expected rawfield as N directionals per face per column in matrix");
	args.failIf(N != 4, "Optimization only supported for N = 4");

	// Setup topology
	EdgeData ed;
	ed.construct(F);

	// The set of parameters for calculating the curl-free fields
	directional::polycurl_reduction_parameters params;

	// Solver data (needed for precomputation)
	directional::PolyCurlReductionSolverData pcrdata;

	//trivial constraints
	Eigen::VectorXi b; b.resize(1); b << 0;
	Eigen::MatrixXd bc; bc.resize(1, 6); bc << rawField.row(0).head(6);
	Eigen::VectorXi blevel; blevel.resize(1); b << 1;
	directional::polycurl_reduction_precompute(V, F, b, bc, blevel, rawField, pcrdata);


	//do a batch of iterations
	dir_LOGLN("--Improving Curl--\n");
	SimpleTimer timer;
	for (int bi = 0; bi < batches; ++bi)
	{
		timer.reset();
		timer.start();
		dir_LOGLN(std::endl << std::endl << "**** Batch " << bi << " ****");
		directional::polycurl_reduction_solve(pcrdata, params, rawField, bi == 0);
		params.wSmooth *= params.redFactor_wsmooth;
		timer.stop();
		std::cout << "Elapsed time:" << timer.elapsed() << std::endl;
		std::cout << "Estimated remaining: " << (batches - bi) * timer.elapsed() << std::endl;
	}

	Eigen::VectorXi matchingCF, singVerticesCF, singIndicesCF, combedMatchingCF;
	VectorXd effortCF, curlCF, combedEffortCF, combedCurlCF;
	MatrixXd combedFieldCF;

	directional::curl_matching(V, F, ed.E, ed.EF, ed.sFE, rawField, matchingCF, effortCF, curlCF);
	directional::effort_to_indices(V, F, ed.E, ed.EF, effortCF, matchingCF, N, singVerticesCF, singIndicesCF);
	directional::combing(V, F, ed.E, ed.EF, ed.sFE, rawField, matchingCF, combedFieldCF);
	directional::curl_matching(V, F, ed.E, ed.EF, ed.sFE, combedFieldCF, combedMatchingCF, combedEffortCF, combedCurlCF);


	const double cfMax =curlCF.maxCoeff();
	const double cfCombMax = combedCurlCF.maxCoeff();
	//Set outputs
	dir_LOGLN("curlMax optimized: " << cfMax);
	dir_LOGLN("curlMax combed optimized: " << cfCombMax);
	args.setOutput(Out_CF_Field, rawField);
	args.setOutput(Out_CF_Matching, matchingCF);
	args.setOutput(Out_CF_Combed_Matching, combedMatchingCF);
	args.setOutput(Out_CF_Combed_Field, combedFieldCF);
	args.setOutput(Out_SingIndices, singIndicesCF);
	args.setOutput(Out_SingVertices, Indexed(singVerticesCF));
	args.setOutput(Out_Max_Curl, cfMax);
	//args.setOutput(Out_EF, Indexed(EF));
}

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