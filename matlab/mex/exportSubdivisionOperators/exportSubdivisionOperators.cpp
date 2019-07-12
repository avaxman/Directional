/* system header */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* MEX header */
#include <mex.h> 
#include "matrix.h"
#include <igl/matlab/MatlabWorkspace.h>
#include "../mexHelpers.h"
#include "../FuncArgs.h"

#include <directional/Subdivision/subdivision.h>
#include <directional/block_diag.h>
#include <directional/Gamma_suite.h>

/**
 * \brief Converts the EdgeData object to a MEX structure
 * \param target Pointer to mxArray pointer, to be set up as mex structure (array).
 * \param ED Edge data object
 * \return MexStructure object that wraps the MEX structure (of size 1x1).
 */
mexhelpers::MexStructure EDToStruct(mxArray** target, const EdgeData& ED)
{
	const char* fields[] = { "EF","sFE","EI","F","E" };
	mexhelpers::MexStructure str(target, fields, 5);
	str.setIndexValue("EF", ED.EF);
	str.setIndexValue("sFE", ED.sFE);
	str.setIndexValue("EI", ED.EI);
	str.setIndexValue("F", ED.F);
	str.setIndexValue("E", ED.E);
	std::cout << "Added edge data struct" << std::endl;
	return str;
}

template<typename T>
std::ostream& writeDims(std::ostream& str, const Eigen::SparseMatrix<T>& t)
{
	str << t.rows() << "," << t.cols() << std::endl;
	return str;
}

/**
 * Export subdivision functions per level
 */
void mainFunc(mexhelpers::FuncArgs& args)
{
	using namespace mexhelpers;

    /* argument check */
	args.verifyInputCount(3, "MATLAB:exportSubdivisionOperators:inputmismatch", "Expected 3 input argument.");
	args.verifyOutputCount(0, "MATLAB:exportSubdivisionOperators:outputmismatch", "Expected 0 output arguments.");

	//args.

	// Acquire face indices
	Eigen::MatrixXi F = args.getIndexMat(0);

	// Acquire level
	int level = 0;
	level = args.getInt(1);

	if(F.cols() != 3)
	{
		mexErrMsgTxt("Expected F to have 3 columns");
	}
	Eigen::MatrixXd dum;
	using SMat = Eigen::SparseMatrix<double>;
	SMat S_v, S_f, S_e, S_c, S_Gamma;
	//std::cout << "Starting subdiv" << std::endl;
	EdgeData ED0, EDL;
	subdivision_operators(F, dum, level, S_v, S_f, S_e, S_c, ED0, EDL);
	//std::cout << "Subdiv done" << std::endl;

	if(EDL.faceCount() == 0)
	{
		mexErrMsgTxt("Invalid face count on new edge data");
	}

	//std::cout << "Setting up output " << std::endl;
	//Export coarse level edge data
	EDToStruct(args.getRawOut(0), ED0);
	//Export fine level edge data
	EDToStruct(args.getRawOut(1), EDL);

	//std::cout << "Setting up subdiv output " << std::endl;
	// Export subdivision operators
	const char* fields[] = { "E","V","C","F" , "Gamma"};
	MexStructure S_Struct(args.getRawOut(2), fields, 5);
	//std::cout << "Setting up subdivs in struct" << std::endl;
	S_Struct.setValue("E", MexSparse::fromEigen(S_e));
	S_Struct.setValue("V", MexSparse::fromEigen(S_v));
	S_Struct.setValue("C", MexSparse::fromEigen(S_c));
	S_Struct.setValue("F", MexSparse::fromEigen(S_f));
	//std::cout << "Setting up gamma output" << std::endl;
	// Build gamma subdivider
	SMat temp, G2ToG30, G3ToG2L, G3ToDC0, DCToG3L;
	//std::cout << "Setting up compound subdiv" << std::endl;
	directional::block_diag({ &S_e, &S_c }, temp);
	std::cout << "S_E ";  writeDims(std::cout, S_e);
	std::cout << "S_C ";  writeDims(std::cout, S_c);
	std::cout << "temp "; writeDims(std::cout, temp);
	if(temp.rows() != S_e.rows()+ S_c.rows())
	{
		mexErrMsgTxt("block_diag produced invalid matrix");
	}

	std::cout << "Building decomposition" << std::endl;
	directional::Gamma3_To_Decomp(ED0.EF, ED0.EI, ED0.faceCount(), G3ToDC0);
	std::cout << "G3ToDC0 "; writeDims(std::cout, G3ToDC0);
	directional::Decomp_To_Gamma3(EDL.EF, EDL.EI,  EDL.faceCount(), DCToG3L);
	std::cout << "DCToG3L "; writeDims(std::cout, DCToG3L);
	directional::Gamma2_To_Gamma3(ED0.sFE, G2ToG30);
	std::cout << "G2ToG30 "; writeDims(std::cout, G2ToG30);
	directional::Gamma3_To_Gamma2(EDL.sFE, G3ToG2L);
	std::cout << "G3ToG2L "; writeDims(std::cout, G3ToG2L);

	std::cout << "Setting output with full matrix" << std::endl;
	// Decompose and apply S_E and S_C subdivision.
	S_Gamma = G3ToG2L * DCToG3L * temp * G3ToDC0 * G2ToG30;
	std::cout << "Assigning to output" << std::endl;
	S_Struct.setValue("Gamma", MexSparse::fromEigen(S_Gamma));
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