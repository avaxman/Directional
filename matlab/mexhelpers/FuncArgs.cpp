#include "FuncArgs.h"
#include "mexHelpers.h"

mexhelpers::FuncArgs::
FuncArgs(mxArray** plhs, int lhsCount, const mxArray** prhs, int rhsCount): prhs(prhs), plhs(plhs), lhsCount(lhsCount),
                                                                            rhsCount(rhsCount)
{
}

void mexhelpers::FuncArgs::verifyInputCount(int num, const char* errId, const char* errmMsg)
{
    if (rhsCount != num)
    {
        mexErrMsgIdAndTxt(errId, errmMsg);
    }
}

void mexhelpers::FuncArgs::verifyInputCount(int num)
{
    if (rhsCount != num)
    {
        std::stringstream msg;
        msg << "Input mismatch: got " << rhsCount << " args, expected " << num << std::endl;
        mexErrMsgTxt(msg.str().c_str());
    }
}

void mexhelpers::FuncArgs::verifyOutputCount(int num)
{
    if (lhsCount != num)
    {
        std::stringstream msg;
        msg << "Output mismatch: got " << rhsCount << " args, expected " << num << std::endl;
        mexErrMsgTxt(msg.str().c_str());
    }
}

void mexhelpers::FuncArgs::failIf(bool condition, const std::string& msg)
{
    if (condition)
    {
        throw std::runtime_error(msg.c_str());
    }
}

void mexhelpers::FuncArgs::verifyOutputCount(int num, const char* errId, const char* errmMsg)
{
    if (lhsCount != num)
    {
        mexErrMsgIdAndTxt(errId, errmMsg);
    }
}

void mexhelpers::FuncArgs::checkIn(int i)
{
    if (i < 0 || i >= rhsCount) throw std::out_of_range(
        std::string("Requested input el") + std::to_string(i) + " out of range");
}

void mexhelpers::FuncArgs::checkOut(int i)
{
    if (i < 0 || i >= lhsCount) throw std::out_of_range(
        std::string("Requested output el") + std::to_string(i) + " out of range");
}



mxArray** mexhelpers::FuncArgs::getRawOut(int i)
{
    checkOut(i);
    return &plhs[i];
}

int mexhelpers::FuncArgs::getInt(int i)
{
    checkIn(i);
    return (int)mxGetScalar(prhs[i]);
}

double mexhelpers::FuncArgs::getDouble(int i)
{
    checkIn(i);
    return mxGetScalar(prhs[i]);
}

Eigen::MatrixXi mexhelpers::FuncArgs::getIndexMat(int i)
{
    checkIn(i);
    Eigen::MatrixXi mat;
    MexArrayReadonly(prhs[i]).mapIndicesTo(mat);
    return mat;
}

Eigen::VectorXi mexhelpers::FuncArgs::getIndexVec(int i)
{
    checkIn(i);
    Eigen::VectorXi mat;
    MexArrayReadonly(prhs[i]).mapIndicesTo(mat);
    return mat;
}

Eigen::MatrixXd mexhelpers::FuncArgs::getMat(int i)
{
    checkIn(i);
    Eigen::MatrixXd mat;
    MexArrayReadonly(prhs[i]).mapTo(mat);
    return mat;
}

Eigen::VectorXd mexhelpers::FuncArgs::getVec(int i)
{
    checkIn(i);
    Eigen::VectorXd mat;
    MexArrayReadonly(prhs[i]).mapTo(mat);
    return mat;
}

Eigen::MatrixXi mexhelpers::FuncArgs::getIntegerMat(int i)
{
    checkIn(i);
    Eigen::MatrixXi mat;
    MexArrayReadonly(prhs[i]).mapTo(mat);
    return mat;
}

Eigen::VectorXi mexhelpers::FuncArgs::getIntegerVec(int i)
{
    checkIn(i);
    Eigen::VectorXi mat;
    MexArrayReadonly(prhs[i]).mapTo(mat);
    return mat;
}
