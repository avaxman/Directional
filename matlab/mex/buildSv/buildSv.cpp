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

enum Input
{
	F,
	Level,
	InputCount
};
enum Output
{
	ED0Struct,
	EDKStruct,
	SV,
	OutputCount
};

void computeSv(const Eigen::MatrixXi& F, const Eigen::MatrixXd& V, int level,
	Eigen::SparseMatrix<double>& S_V,
	EdgeData& ED0,
	EdgeData& EDL)
{
	std::cout << "### Subdivision construction" << std::endl;
	SimpleTimer totalTime;
	totalTime.start();
	SimpleTimer st;
	st.start();

	// Subdividers for different field types
	VertexFieldSubdivider<LoopCoefficientProvider> vSub;
	// Create the subdivision builder
	auto builder = create_sub_builder(F, vSub);
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
	st.stop();
	std::cout << "Getting matrices in " << st.elapsed() << " ms" << std::endl;

	totalTime.stop();
	std::cout << "Total time : " << totalTime.elapsed() << " ms" << std::endl;
	// Save the fine level topology
	EDL = builder.ED;
}

void mainFunc(mexhelpers::FuncArgs& args)
{
	using namespace mexhelpers;
	

    /* argument check */
	args.verifyInputCount(Input::InputCount, "MATLAB:buildSubdivisionOperators:inputmismatch", "Expected 2 input argument.");
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
	computeSv(F, dum, level, S_v, ED0, EDL);
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

	args.setOutput(2, MexSparse::fromEigen(S_v));
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