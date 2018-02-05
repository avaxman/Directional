#ifndef DIRECTIONAL_READ_SINGULARITIES_H
#define DIRECTIONAL_READ_SINGULARITIES_H
#include <cmath>
#include <Eigen/Core>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>


namespace directional
{

	// Reads a list of singularities from a file
	// Inputs:
	//   fileName: The to be loaded file.
	// Outputs:
	//   singularities: The vector containing the singularities
	//   N: The degree of the field
	//   globalRotation: The angle of rotation between the vector on the first face and its basis in radians
	// Return:
	//   Whether or not the file was written successfully
  bool IGL_INLINE read_singularities(const std::string &fileName, int& N, Eigen::VectorXi &singPositions, Eigen::VectorXi& singIndices)
	{
		try
		{
			std::ifstream f(fileName);
			int numSings;
			f >> N;
			f >> numSings;

			singPositions = Eigen::VectorXi::Zero(numSings);
      singIndices = Eigen::VectorXi::Zero(numSings);
      
      for (int i=0;i<numSings;i++)
        f >> singPositions.coeffRef(i)>> singIndices.coeffRef(i);
      
			f.close();
			return f.fail();
		}
		catch (std::exception e)
		{
			return false;
		}
	}
}

#endif
