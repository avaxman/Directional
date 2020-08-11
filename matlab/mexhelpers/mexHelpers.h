#ifndef MEXHELPERS_H
#define MEXHELPERS_H
#include "mex.h"
#include "matrix.h"
#include <Eigen/Eigen>
#include <Eigen/SparseCore>
#include <type_traits>
#include <limits>
#include <iostream>
#include <sstream>
#include "FuncArgs.h"

#ifndef VERBOSE
#define dir_LOGLN(expr) 
#else
#define dir_LOGLN(expr) std::cout << expr << std::endl ;
#endif

#define DIR_THROW(clsType, text) throw clsType(text ## "at " ## __FILE__ ":" ## __LINE__)

#define DIR_MEX_ENTRY(mainFuncName) void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {\\
mexhelpers::FuncArgs args(plhs, nlhs, prhs, nrhs);\\
try{ mainFuncName##(args); } \\
catch (const std::exception& e) { mexErrMsgTxt(e.what()); }\\
}

// Derived from https://stackoverflow.com/questions/49952275/passing-sparse-arrays-from-matlab-to-eigen-c-and-back-to-matlab
namespace mexhelpers {

	class Logger
	{
		std::ostream* m_str;
		Logger(): activeLevel(Info)
		{
			m_str = &std::cout;
			globalLevel = None; // Default off
		}
	public:
		enum Level
		{
			Info = 0,
			Warning,
			None // = errors only
		};
		static char nl()
		{
			return '\n';
		}
		static void setLevel(Level level)
		{
			log().globalLevel = level;
		}
		static Logger& log(Level level = Info)
		{
			static Logger l;
			return l;
		}
		template<typename T>
		Logger& operator<<(T t)
		{
#ifndef MEX_DIR_NO_LOG
			if(activeLevel >= globalLevel)
				*m_str << t;
#endif
			return *this;
		}
	private:
		Level activeLevel;
		Level globalLevel;
	};

	/**
	 * \brief Wrapper for read-only mex array inputs
	 */
	struct MexArrayReadonly
	{
		const mxArray* arr;
		MexArrayReadonly(const mxArray* arr) : arr(arr) {}

		/**
		 * \brief Number of rows
		 * \return Returns the number of rows
		 */
		size_t rows() const
		{
			return mxGetM(arr);
		}


		double* data()
		{
			return mxGetPr(arr);
		}
		
		/**
		 * \brief Number of cols
		 * \return Returns the number of columns
		 */
		size_t cols() const
		{
			return mxGetN(arr);
		}

		/**
		 * \brief Maps the values of the mxarray as indices to an Eigen vector.
		 * Compensates for the 0-1 indexed problem between Matlab and C/C++.
		 * \param vec The output Eigen vector
		 */
		void mapIndicesTo(Eigen::VectorXi& vec)
		{
			vec = Eigen::VectorXi(std::max(rows(), cols()));
			const int rs = rows();
			const double* dat = data();
			// Column major.
			for (int c = 0; c < vec.size(); c++)
			{
				vec(c) = (int)dat[c] - 1;
			}
		}
		/**
		 * \brief Maps the values of the mxarray as indices to an Eigen matrix.
		 * Compensates for the 0-1 indexed problem between Matlab and C/C++.
		 * \param mat The output Eigen matrix
		 */
		void mapIndicesTo(Eigen::MatrixXi& mat)
		{
			mat = Eigen::MatrixXi(rows(), cols());
			const int rs = rows();
			const double* dat = data();
			// Column major.
			for (int c = 0; c < cols(); c++)
			{
				for (int r = 0; r < rs; r++)
				{
					mat(r, c) = (int)dat[r + c * rs] - 1;
				}
			}
		}
		void mapTo(Eigen::MatrixXd& mat)
		{
			mat = Eigen::MatrixXd(rows(), cols());
			const int rs = rows();
			const double* dat = data();
			// Column major.
			for (int c = 0; c < cols(); c++)
			{
				for (int r = 0; r < rs; r++)
				{
					mat(r, c) = dat[r + c * rs];
				}
			}
		}
		void mapTo(Eigen::VectorXd& mat)
		{
			const int size = rows() * cols();
			mat = Eigen::VectorXd(size);
			const int rs = rows();
			const double* dat = data();
			// Column major.
			for (int i = 0; i < size; i++)
			{
				mat(i) = dat[i];
			}
		}
		void mapTo(Eigen::VectorXi& mat)
		{
			const int size = rows() * cols();
			mat = Eigen::VectorXi(size);
			const int rs = rows();
			const double* dat = data();
			// Column major.
			for (int i = 0; i < size; i++)
			{
				mat(i) = dat[i];
			}
		}
		void mapTo(Eigen::MatrixXi& mat)
		{
			mat = Eigen::MatrixXi(rows(), cols());
			const int rs = rows();
			const double* dat = data();
			// Column major.
			for (int c = 0; c < cols(); c++)
			{
				for (int r = 0; r < rs; r++)
				{
					mat(r, c) = dat[r + c * rs];
				}
			}
		}
	};

	/**
	 * Wrapper to signal that the MEX element represents an index (possibly requires conversion between 0- and 1-based indexing).
	 */
	template<typename T>
	struct Indexed_type
	{
		T& t;
		Indexed_type(T& t) : t(t) {}
	};

	template<typename T>
	Indexed_type<T> Indexed(T& t)
	{
		return { t };
	}

	/**
	 * Wrapper for mutable MEX array
	 */
	struct MexArray
	{
		mxArray* arr;
		// Wraps an mx array
		MexArray(mxArray* arr) : arr(arr) {}

		//Explicitly creates an mx array
		MexArray(int rows, int cols)
		{
			arr = mxCreateDoubleMatrix(rows, cols, mxREAL);
			Logger::log() << "Mex arr created: " << rows << ',' << cols << Logger::nl();
		}

		double* data()
		{
			return mxGetPr(arr);
		}

		size_t rows() const
		{
			return mxGetM(arr);
		}
		size_t cols() const
		{
			return mxGetN(arr);
		}

		void mapIndicesTo(Eigen::VectorXi& vec)
		{
			vec = Eigen::VectorXi(std::max(rows(), cols()));
			const int rs = rows();
			const double* dat = data();
			// Column major.
			for (int c = 0; c < vec.size(); c++)
			{
				vec(c) = (int)dat[c] - 1;
			}
		}
		void mapIndicesTo(Eigen::MatrixXi& mat)
		{
			mat = Eigen::MatrixXi(rows(), cols());
			const int rs = rows();
			const double* dat = data();
			// Column major.
			for (int c = 0; c < cols(); c++)
			{
				for (int r = 0; r < rs; r++)
				{
					mat(r, c) = (int)dat[r + c * rs] - 1;
				}
			}
		}
		template<typename Derived>
		static MexArray fromEigenIndex(const Eigen::MatrixBase<Derived>& vec)
		{
			Logger::log() << "Setting up regular Eigen index to matlab" << Logger::nl();
			MexArray ma(vec.rows(), vec.cols());
			double* dat = ma.data();
			for (int i = 0; i < vec.cols(); i++)
			{
				for (int j = 0; j < vec.rows(); j++)
				{
					dat[i*vec.rows() + j] = vec(j, i) + 1;
				}
			}
			Logger::log() << "Converted Eigen vec to mex" << Logger::nl();
			return ma;
		}

		template<typename Derived>
		static MexArray fromEigen(const Eigen::MatrixBase<Derived>& vec)
		{
			Logger::log() << "Setting up regular Eigen to matlab" << Logger::nl();
			MexArray ma(vec.rows(), vec.cols());
			double* dat = ma.data();
			for (int i = 0; i < vec.cols(); i++)
			{
				for (int j = 0; j < vec.rows(); j++)
				{
					dat[i*vec.rows() + j] = vec(j, i);
				}
			}
			std::cout << "Converted Eigen vec to mex" << std::endl;
			return ma;
		}
	};


	using namespace Eigen;
	// Matlab sparse type in Eigen representation
	using MatlabSparse = Eigen::SparseMatrix<double, Eigen::ColMajor, std::make_signed<mwIndex>::type>;

	/**
	 * Wrapper for mutable sparse MEX array.
	 */
	struct MexSparse
	{
		mxArray* arr;
		MexSparse(mxArray* arr) :arr(arr) {}

		MexSparse(int rows, int cols, int nnz)
		{
			arr = mxCreateSparse(rows, cols, nnz, mxREAL);
			std::cout << "Created spare" << std::endl;
		}

		size_t rows() const
		{
			return mxGetM(arr);
		}

		size_t cols() const
		{
			return mxGetN(arr);
		}

		double* data()
		{
			return mxGetPr(arr);
		}

		Eigen::MappedSparseMatrix<double, ColMajor, long long> eigenView()
		{
			mxAssert(mxGetClassID(arr) == mxDOUBLE_CLASS,
				"Type of the input matrix isn't double");
			mwSize     m = rows();
			mwSize     n = cols();
			mwSize    nz = mxGetNzmax(arr);
			/*Theoretically fails in very very large matrices*/
			mxAssert(nz <= std::numeric_limits< std::make_signed<mwIndex>::type>::max(),
				"Unsupported Data size."
			);
			double  * pr = data();
			auto* ir = reinterpret_cast<MatlabSparse::StorageIndex*>(mxGetIr(arr));
			auto* jc = reinterpret_cast<MatlabSparse::StorageIndex*>(mxGetJc(arr));
			Eigen::MappedSparseMatrix<double, ColMajor, long long> result(m, n, nz, jc, ir, pr);
			return result;
		}
		static MexSparse fromEigen(const Eigen::SparseMatrix<double, ColMajor, long long>& mat)
		{
			std::cout << "Setting up a sparse mat" << std::endl;
			MexSparse ms(mat.rows(), mat.cols(), mat.nonZeros());
			const MatlabSparse::StorageIndex* ir = mat.innerIndexPtr();
			const MatlabSparse::StorageIndex* jc = mat.outerIndexPtr();
			const double* pr = mat.valuePtr();

			mwIndex * ir2 = mxGetIr(ms.arr);
			mwIndex * jc2 = mxGetJc(ms.arr);
			double  * pr2 = mxGetPr(ms.arr);
			std::cout << "Copying data" << std::endl;
			for (mwIndex i = 0; i < mat.nonZeros(); i++) {
				pr2[i] = pr[i];
				ir2[i] = ir[i];
			}
			for (mwIndex i = 0; i < mat.cols() + 1; i++) {
				jc2[i] = jc[i];
			}
			return ms;
		}
	};

	Eigen::MappedSparseMatrix<double, ColMajor, long long>
		matlab_to_eigen_sparse(const mxArray * mat)
	{
		mxAssert(mxGetClassID(mat) == mxDOUBLE_CLASS,
			"Type of the input matrix isn't double");
		mwSize     m = mxGetM(mat);
		mwSize     n = mxGetN(mat);
		mwSize    nz = mxGetNzmax(mat);
		/*Theoretically fails in very very large matrices*/
		mxAssert(nz <= std::numeric_limits< std::make_signed<mwIndex>::type>::max(),
			"Unsupported Data size."
		);
		double  * pr = mxGetPr(mat);
		MatlabSparse::StorageIndex* ir = reinterpret_cast<MatlabSparse::StorageIndex*>(mxGetIr(mat));
		MatlabSparse::StorageIndex* jc = reinterpret_cast<MatlabSparse::StorageIndex*>(mxGetJc(mat));
		Eigen::MappedSparseMatrix<double, ColMajor, long long> result(m, n, nz, jc, ir, pr);
		return result;
	}

	mxArray*
		eigen_to_matlab_sparse(const MatlabSparse& mat)
	{
		mxArray * result = mxCreateSparse(mat.rows(), mat.cols(), mat.nonZeros(), mxREAL);
		const MatlabSparse::StorageIndex* ir = mat.innerIndexPtr();
		const MatlabSparse::StorageIndex* jc = mat.outerIndexPtr();
		const double* pr = mat.valuePtr();

		mwIndex * ir2 = mxGetIr(result);
		mwIndex * jc2 = mxGetJc(result);
		double  * pr2 = mxGetPr(result);

		for (mwIndex i = 0; i < mat.nonZeros(); i++) {
			pr2[i] = pr[i];
			ir2[i] = ir[i];
		}
		for (mwIndex i = 0; i < mat.cols() + 1; i++) {
			jc2[i] = jc[i];
		}
		return result;
	}



	/**
	 * \brief Creates an mxArray object that Matlab can work with for the given type
	 * \tparam T The type
	 * \param t The value
	 * \return The matlab object
	 */
	template<typename T>
	struct MexConvert
	{
		static mxArray* toMex(const T& t)
		{
			std::stringstream ss;
			ss << "MexConvert::toMex not implemented! " << typeid(T).name();
			//static_assert(false, "MexConvert::toMax should be implemented");
			throw std::runtime_error(ss.str().c_str());
			
		}
	};

	template<> inline mxArray* MexConvert<MexArray>::toMex(const MexArray& str)
	{
		return str.arr;
	}
	template<> inline mxArray* MexConvert<MexSparse>::toMex(const MexSparse& str)
	{
		return str.arr;
	}
	template<> inline mxArray* MexConvert<double>::toMex(const double& t)
	{
		return mxCreateDoubleScalar(t);
	}
	template<> inline mxArray* MexConvert<int>::toMex(const int& t)
	{
		// Save it as double anyway...
		return mxCreateDoubleScalar(t);
	}

	template<typename T>
	struct MexConvert<Indexed_type<T>>
	{
		static mxArray* toMex(const Indexed_type<T>& i)
		{
			//static_assert(false, "Matched Indexed");
			return MexArray::fromEigenIndex(i.t).arr;
		}
	};
	template<typename Scalar, int RowsAtCompile, int ColsAtCompile, int Opts, int MaxR, int MaxC >
	struct MexConvert<Eigen::Matrix<Scalar, RowsAtCompile, ColsAtCompile, Opts, MaxR, MaxC>>
	{
		static mxArray* toMex(const Eigen::Matrix<Scalar, RowsAtCompile, ColsAtCompile, Opts, MaxR, MaxC>& i)
		{
			//static_assert(false, "Matched MatrixBase");
			return MexArray::fromEigen(i).arr;
		}
	};

	template<typename Scalar, int Opts, typename Index>
	struct MexConvert<Eigen::SparseMatrix<Scalar, Opts, Index>>
	{
		static mxArray* toMex(const Eigen::SparseMatrix<Scalar, Opts, Index>& i)
		{
			//static_assert(false, "Matched MatrixBase");
			return MexSparse::fromEigen(i).arr;
		}
	};
	
	template<>
	struct MexConvert<mxArray>
	{
		static mxArray* toMex(mxArray* arr)
		{
			return arr;
		}
	};


	struct MexStructure
	{
		mxArray* target;
		MexStructure(mxArray** target, const char** fields, int fieldNum)
		{
			mwSize sz[] = { 1,1 };
			this->target = mxCreateStructArray(2, sz, fieldNum, fields);
			*target = this->target;
			std::cout << "Struct array set up" << std::endl;
			mxAssert(mxGetNumberOfFields(*target) == fieldNum, "Invalid fiedlnum");
		}
		MexStructure(const char** fields, int fieldNum)
		{
			mwSize sz[] = { 1,1 };
			this->target = mxCreateStructArray(2, sz, fieldNum, fields);
			std::cout << "Struct array set up" << std::endl;
			mxAssert(mxGetNumberOfFields(target) == fieldNum, "Invalid fiedlnum");
		}
		template<typename T>
		void setValue(const char* name, const T& t)
		{
			mxSetField(target, 0, name, MexConvert<T>::toMex(t));
		}
		template<typename T>
		void setIndexValue(const char* name, const T& t)
		{
			// To be specialized.
		}
	};


	template<> inline void MexStructure::setIndexValue<Eigen::MatrixXi>(const char* name, const Eigen::MatrixXi& mat)
	{
		mxSetField(target, 0, name, MexArray::fromEigenIndex(mat).arr);
	}
	template<> inline void MexStructure::setIndexValue<Eigen::VectorXi>(const char* name, const Eigen::VectorXi& mat)
	{
		mxSetField(target, 0, name, MexArray::fromEigenIndex(mat).arr);
	}

	template<> inline mxArray* MexConvert<MexStructure>::toMex(const MexStructure& str)
	{
		return str.target;
	}

    class MexFunction
	{
	protected:
        mexhelpers::FuncArgs m_args;
        // To be overridden
        virtual void apply(const mexhelpers::FuncArgs& funcArgs) = 0;

        /**
         * \brief Throws an error with the given message if the given condition does not hold
         * \param num Number of expected outputs
         */
        void failIf(bool condition, const std::string& msg)
        {
            if (condition)
            {
                throw std::runtime_error(msg.c_str());
            }
        }
        void verifyInputCount(int num, const char* errId, const char* errmMsg)
        {
            if (m_args.rhsCount != num)
            {
                mexErrMsgIdAndTxt(errId, errmMsg);
            }
        }
        void verifyOutputCount(int num, const char* errId, const char* errmMsg)
        {
            if (rhsCount != num)
            {
                mexErrMsgIdAndTxt(errId, errmMsg);
            }
        }
	public:
        virtual ~MexFunction() = default;

        void execute(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
        {
            try
            {
                mexhelpers::FuncArgs args(plhs, nlhs, prhs, nrhs);
                apply(args);
            }
            catch (const std::exception& e)
            {
                mexErrMsgTxt(e.what());
            }
            
        }
	};
#define MEX_ENTRY(clsName) void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])\
    {\
        runner.execute(nlhs, plhs, nrhs, prhs);\
    }
}
#endif