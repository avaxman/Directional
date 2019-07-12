#ifndef MEXHELPERS_FUNCARGS_H
#define MEXHELPERS_FUNCARGS_H
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
		FuncArgs(mxArray** plhs, int lhsCount, const mxArray** prhs, int rhsCount) :prhs(prhs), plhs(plhs), lhsCount(lhsCount), rhsCount(rhsCount) {}

		/**
		 * \brief Verifies the number of inputs. Returns a mex error with the given ID and message if the input count is incorrect.
		 * \param num Number of expected inputs
		 * \param errId Err ID
		 * \param errmMsg Error message to display
		 */
		void verifyInputCount(int num, const char* errId, const char* errmMsg)
		{
			if (rhsCount != num)
			{
				mexErrMsgIdAndTxt(errId, errmMsg);
			}
		}

		/**
		 * \brief Verifies the number of inputs. Generates a MEX error and quits if the number is incorrect.
		 * \param num Number of expected inputs
		 */
		void verifyInputCount(int num)
		{
			if(rhsCount != num)
			{
				std::stringstream msg;
				msg << "Input mismatch: got " << rhsCount << " args, expected " << num << std::endl;
				mexErrMsgTxt(msg.str().c_str());
			}
		}

		/**
		 * \brief Verifies the number of outputs. Generates a MEX error and quits if the number is incorrect.
		 * \param num Number of expected outputs
		 */
		void verifyOutputCount(int num)
		{
			if (lhsCount != num)
			{
				std::stringstream msg;
				msg << "Output mismatch: got " << rhsCount << " args, expected " << num << std::endl;
				mexErrMsgTxt(msg.str().c_str());
			}
		}

		/**
		 * \brief Throws an error with the given message if the given condition does not hold
		 * \param num Number of expected outputs
		 */
		void failIf(bool condition, const std::string& msg)
		{
			if(condition)
			{
				throw std::runtime_error(msg.c_str());
			}
		}

		/**
		 * \brief Verifies the number of outputs. Returns a mex error with the given ID and message if the output count is incorrect.
		 * \param num Number of expected outputs
		 * \param errId Err ID
		 * \param errmMsg Error message to display
		 */
		void verifyOutputCount(int num, const char* errId, const char* errmMsg)
		{
			if (lhsCount != num)
			{
				mexErrMsgIdAndTxt(errId, errmMsg);
			}
		}

		/**
		 * \brief Checks whether the index is not out of bounds for inputs. Throws otherwise
		 * \param i The index
		 */
		void checkIn(int i)
		{
			if (i < 0 || i >= rhsCount) throw std::out_of_range(std::string("Requested input el") + std::to_string(i) + " out of range");
		}

		/**
		 * \brief Checks whether the index is not out of bounds for outputs. Throws otherwise
		 * \param i The index
		 */
		void checkOut(int i)
		{
			if (i < 0 || i >= lhsCount) throw std::out_of_range(std::string("Requested output el") + std::to_string(i) + " out of range");
		}

		/**
		 * Sets the output argument at the specified index to the given value.
		 * \param i The output argument index
		 */
		template<typename T>
		void setOutput(int i, const T& value)
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
		mxArray** getRawOut(int i)
		{
			checkOut(i);
			return &plhs[i];
		}

		/**
		 * \brief Retrieves the input element at the given index as an integer
		 * \return The input integer
		 */
		int getInt(int i)
		{
			checkIn(i);
			return (int)mxGetScalar(prhs[i]);
		}

		/**
		 * \brief Retrieves the input element at the given index as a double
		 * \return The input double
		 */
		double getDouble(int i)
		{
			checkIn(i);
			return mxGetScalar(prhs[i]);
		}

		/**
		 * \brief Retrieves the input element at the given index as a double
		 * \return The input double
		 */
		Eigen::MatrixXi getIndexMat(int i)
		{
			checkIn(i);
			Eigen::MatrixXi mat;
			MexArrayReadonly(prhs[i]).mapIndicesTo(mat);
			return mat;
		}
		Eigen::VectorXi getIndexVec(int i)
		{
			checkIn(i);
			Eigen::VectorXi mat;
			MexArrayReadonly(prhs[i]).mapIndicesTo(mat);
			return mat;
		}
		Eigen::MatrixXd getMat(int i)
		{
			checkIn(i);
			Eigen::MatrixXd mat;
			MexArrayReadonly(prhs[i]).mapTo(mat);
			return mat;
		}
		Eigen::VectorXd getVec(int i)
		{
			checkIn(i);
			Eigen::VectorXd mat;
			MexArrayReadonly(prhs[i]).mapTo(mat);
			return mat;
		}
		Eigen::MatrixXi getIntegerMat(int i)
		{
			checkIn(i);
			Eigen::MatrixXi mat;
			MexArrayReadonly(prhs[i]).mapTo(mat);
			return mat;
		}
		Eigen::VectorXi getIntegerVec(int i)
		{
			checkIn(i);
			Eigen::VectorXi mat;
			MexArrayReadonly(prhs[i]).mapTo(mat);
			return mat;
		}
	};
}
#endif