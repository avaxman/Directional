#ifndef READ_TRIVIAL_FIELD
#define READ_TRIVIAL_FIELD
#include <cmath>
#include <igl/read_triangle_mesh.h>
#include <igl/writeOFF.h>
#include <igl/boundary_loop.h>
#include <directional/dual_cycles.h>
#include <Eigen/Core>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>


namespace directional
{
	// Writes a list of singularities to a file.
	// Inputs:
	//   fileName: The file name to which the singularities should be saved.
	//   singularities: The vector containing the singularities
	//   N: The degree of the field
	//   globalRotation: The angle of rotation between the vector on the first face and its basis in radians
	// Return:
	//   Whether or not the file was written successfully
	bool IGL_INLINE write_singularities(const std::string &fileName, const Eigen::VectorXi &singularities, const int N, const double globalRotation)
	{
		try
		{
			std::ofstream f(fileName, std::ios::trunc);
			f << N << " " << singularities.size() << " " << globalRotation << std::endl;
			for (int i = 0; i < singularities.rows(); i++)
			{
				if (singularities(i))
					f << i << " " << singularities(i) << std::endl;
			}
			f.close();
			return !f.fail();
		}
		catch(std::exception e)
		{
			return false;
		}
	}

	// Reads a list of singularities from a file
	// Inputs:
	//   fileName: The to be loaded file.
	// Outputs:
	//   singularities: The vector containing the singularities
	//   N: The degree of the field
	//   globalRotation: The angle of rotation between the vector on the first face and its basis in radians
	// Return:
	//   Whether or not the file was written successfully
	bool IGL_INLINE read_singularities(const std::string &fileName, Eigen::VectorXi &singularities, int& N, double& globalRotation)
	{
		try
		{
			std::ifstream f(fileName);
			int s;
			f >> N;
			f >> s;
			f >> globalRotation;

			singularities = Eigen::VectorXi::Zero(s);

			while(f)
			{
				int i;
				f >> i;
				f >> singularities.coeffRef(i);
			}
			f.close();
			return f.fail();
		}
		catch (std::exception e)
		{
			return false;
		}
	}

	// Writes a list of singularities and the mesh to files. The method will create a folder containing a "mesh.off" file with the mesh and "singularities.sing" file with the singularities.
	// Inputs:
	//   filename: The name used for the mesh and singularity file, without extention
	//   V: List of Vertices
	//   F: List of faces
	//   singularities: The vector containing the singularities
	//   N: The degree of the field
	//   globalRotation: The angle of rotation between the vector on the first face and its basis in radians
	// Return:
	//   Whether or not the file was written successfully
	bool IGL_INLINE write_trivial_field(const std::string filename, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXi &singularities, const int N, const double globalRotation)
	{
		return igl::writeOFF(filename + ".off", V, F) &&
			directional::write_singularities(filename + ".sing", singularities, N, globalRotation);
	}

	// Will search for a mesh file "mesh.off" and a singularity file "singularities.sing" in the given folder and loads their data.
	// Inputs:
	//   filename: The name used for the mesh and singularity file, without extention
	// Return:
	//   V: List of Vertices
	//   F: List of faces
	//   singularities: The vector containing the singularities
	//   N: The degree of the field
	//   globalRotation: The angle of rotation between the vector on the first face and its basis in radians
	//   Whether or not the file was written successfully
	bool IGL_INLINE read_trivial_field(const std::string& filename, Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::VectorXi &singularities, int N, double globalRotation)
	{
		return igl::readOFF(filename + ".off", V, F) &&
			directional::read_singularities(filename + ".sing", singularities, N, globalRotation);
	}


}

#endif