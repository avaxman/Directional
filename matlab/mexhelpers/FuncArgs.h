#ifndef MEXHELPERS_FUNCARGS_H
#define MEXHELPERS_FUNCARGS_H
#include "mex.h"
#include "matrix.h"
#include <sstream>
#include <Eigen/Eigen>

namespace mexhelpers
{
	/**
	 * Structure for handling input and output for a MEX function
	 */
	struct FuncArgs
	{
		// Pointers to mxArray pointers of the input arguments
		const mxArray ** prhs;

		// Pointers to mxArray pointers of the output arguments
		mxArray** plhs;

		// Number of outputs and inputs
		int lhsCount, rhsCount;

		/**
		 * \brief Constructs the function arguments helper structure.
		 * \param plhs Pointer to pointer of mxArray objects, representing the left hand side (output) of the function
		 * \param lhsCount Number of output arguments
		 * \param prhs Pointer to pointer of mxArray objects, representing the right hand side (input) of the function
		 * \param rhsCount Number of input arguments
		 */
        FuncArgs(mxArray** plhs, int lhsCount, const mxArray** prhs, int rhsCount);

		/**
		 * \brief Verifies the number of inputs. Returns a mex error with the given ID and message if the input count is incorrect.
		 * \param num Number of expected inputs
		 * \param errId Err ID
		 * \param errmMsg Error message to display
		 */
        void verifyInputCount(int num, const char* errId, const char* errmMsg);

		/**
		 * \brief Verifies the number of inputs. Generates a MEX error and quits if the number is incorrect.
		 * \param num Number of expected inputs
		 */
        void verifyInputCount(int num);

		/**
		 * \brief Verifies the number of outputs. Generates a MEX error and quits if the number is incorrect.
		 * \param num Number of expected outputs
		 */
        void verifyOutputCount(int num);

		/**
		 * \brief Throws an error with the given message if the given condition does not hold
		 * \param num Number of expected outputs
		 */
        void failIf(bool condition, const std::string& msg);

		/**
		 * \brief Verifies the number of outputs. Returns a mex error with the given ID and message if the output count is incorrect.
		 * \param num Number of expected outputs
		 * \param errId Err ID
		 * \param errmMsg Error message to display
		 */
        void verifyOutputCount(int num, const char* errId, const char* errmMsg);

		/**
		 * \brief Checks whether the index is not out of bounds for inputs. Throws otherwise
		 * \param i The index
		 */
        void checkIn(int i);

		/**
		 * \brief Checks whether the index is not out of bounds for outputs. Throws otherwise
		 * \param i The index
		 */
        void checkOut(int i);

		/**
		 * Sets the output argument at the specified index to the given value.
		 * \param i The output argument index
		 */
		template<typename T>
        void setOutput(int i, const T& value);
		{
			checkOut(i);
			std::cout << "Setting output " << i << std::endl;
			mxArray* arr = MexConvert<T>::toMex(value);
			if (arr == nullptr) std::cout << "NULLPTR!" << std::endl;
			plhs[i] = arr;
		}


		/**
		 * \brief Returns a pointer to the mxArray element pointer, if the index is in bounds. Throws otherwise
		 * \param i The index of the output element to return
		 * \return The pointer to mxArray pointer.
		 */
        mxArray** getRawOut(int i);

		/**
		 * \brief Retrieves the input element at the given index as an integer
		 * \return The input integer
		 */
        int getInt(int i);

		/**
		 * \brief Retrieves the input element at the given index as a double
		 * \return The input double
		 */
        double getDouble(int i);

		/**
		 * \brief Retrieves the input element at the given index as a double
		 * \return The input double
		 */
        Eigen::MatrixXi getIndexMat(int i);

        Eigen::VectorXi getIndexVec(int i);

        Eigen::MatrixXd getMat(int i);

        Eigen::VectorXd getVec(int i);

        Eigen::MatrixXi getIntegerMat(int i);

        Eigen::VectorXi getIntegerVec(int i);
	};
}
#endif