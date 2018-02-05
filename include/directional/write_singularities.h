#ifndef DIRECTIONAL_WRITE_SINGULARITIES_H
#define DIRECTIONAL_WRITE_SINGULARITIES_H
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
	bool IGL_INLINE write_singularities(const std::string &fileName, const int N, const Eigen::VectorXi &singIndices, const Eigen::VectorXi &singPositions)
	{
		try
		{
			std::ofstream f(fileName, std::ios::trunc);
			f<<N<<" "<<singIndices.size()<<std::endl;
      
      for (int i=0;i<singIndices.rows();i++)
        f<<singPositions(i)<<" "<<singIndices(i)<<std::endl;
      
			f.close();
			return !f.fail();
		}
		catch(std::exception e)
		{
			return false;
		}
	}
}

#endif
