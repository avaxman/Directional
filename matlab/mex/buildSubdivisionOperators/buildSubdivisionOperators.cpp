/* system header */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* MEX header */
#include <mex.h> 
#include "matrix.h"
#include "../mexHelpers.h"
#include "../FuncArgs.h"
#include <directional/Subdivision/subdivision.h>
#include <directional/block_diag.h>
#include <directional/Gamma_suite.h>

mexhelpers::MexStructure EDToStruct(mxArray** target, const EdgeData& ED)
{
	const char* fields[] = { "EF","sFE","EI","F","E" };
	mexhelpers::MexStructure str(target, fields, 5);
	str.setIndexValue("EF", ED.EF);
	str.setIndexValue("sFE", ED.sFE);
	str.setIndexValue("EI", ED.EI);
	str.setIndexValue("F", ED.F);
	str.setIndexValue("E", ED.E);
	dir_LOGLN("Added edge data struct");
	return str;
}
template<typename EigenType>
struct Dims
{
	const EigenType& type;
	Dims(const EigenType& type): type(type){}
	int rows()const
	{
		return type.rows();
	}
	int cols()const
	{
		return type.rows();
	}
};
template<typename EigenType>
std::ostream& operator<<(std::ostream& str, const Dims<EigenType>& type)
{
	str << type.rows() << "," << type.cols();
	return str;
}
template<typename EigenType>
Dims<EigenType> dims(const EigenType& type)
{
	return Dims<EigenType>(type);
}

void mainFunc(mexhelpers::FuncArgs& args)
{
	using namespace mexhelpers;
	

    /* argument check */
	args.verifyInputCount(2, "MATLAB:buildSubdivisionOperators:inputmismatch", "Expected 2 input argument.");
	args.verifyOutputCount(3, "MATLAB:buildSubdivisionOperators:outputmismatch", "Expected 3 output argument.");

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

	dir_LOGLN("Setting up subdiv output ");
	// Export subdivision operators
	const char* fields[] = { "E","V","C","F" , "Gamma"};
	MexStructure S_Struct(args.getRawOut(2), fields, 5);
	//std::cout << "Setting up subdivs in struct" << std::endl;
	S_Struct.setValue("E", MexSparse::fromEigen(S_e));
	S_Struct.setValue("V", MexSparse::fromEigen(S_v));
	S_Struct.setValue("C", MexSparse::fromEigen(S_c));
	S_Struct.setValue("F", MexSparse::fromEigen(S_f));
	dir_LOGLN("Setting up gamma output");
	// Build gamma subdivider
	SMat temp, G2ToG30, G3ToG2L, G3ToDC0, DCToG3L;
	dir_LOGLN("Setting up compound subdiv");
	directional::block_diag({ &S_e, &S_c }, temp);
	/*std::cout << "S_E ";  writeDims(std::cout, S_e);
	std::cout << "S_C ";  writeDims(std::cout, S_c);
	std::cout << "temp "; writeDims(std::cout, temp);*/
	if(temp.rows() != S_e.rows()+ S_c.rows())
	{
		mexErrMsgTxt("block_diag produced invalid matrix");
	}

	dir_LOGLN("Building decomposition");
	directional::Gamma3_To_Decomp(ED0.EF, ED0.EI, ED0.faceCount(), G3ToDC0);
	dir_LOGLN("G3ToDC0 " << dims(G3ToDC0));
	directional::Decomp_To_Gamma3(EDL.EF, EDL.EI,  EDL.faceCount(), DCToG3L);
	dir_LOGLN("DCToG3L " << dims(DCToG3L));
	directional::Gamma2_To_Gamma3(ED0.sFE, G2ToG30);
	dir_LOGLN("G2ToG30 " << dims(G2ToG30));
	directional::Gamma3_To_Gamma2(EDL.sFE, G3ToG2L);
	dir_LOGLN("G3ToG2L " << dims(G3ToG2L));

	dir_LOGLN("Setting output with full matrix");
	// Decompose and apply S_E and S_C subdivision.
	S_Gamma = G3ToG2L * DCToG3L * temp * G3ToDC0 * G2ToG30;
	dir_LOGLN("Assigning to output");
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