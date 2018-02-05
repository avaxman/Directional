#ifndef DIRECTIONAL_WRITE_RAW_FIELD_H
#define DIRECTIONAL_WRITE_RAW_FIELD_H

#include <Eigen/Core>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fstream>

namespace directional
{

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
	bool IGL_INLINE write_raw_field(const std::string fileName, Eigen::MatrixXd& rawField)
	{
    
  try
  {
    std::ofstream f(fileName);
    
    int N=rawField.cols()/3;
    f<<N<<" "<<rawField.rows()<<std::endl;
    for (int i=0;i<rawField.rows();i++){
      for (int j=0;j<rawField.cols();j++)
        f<<rawField(i,j)<<" ";
      f<<std::endl;
    }
    f.close();
    return !f.fail();
  }
  catch (std::exception e)
  {
    return false;
  }
}
  
}

#endif
