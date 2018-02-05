#ifndef DIRECTIONAL_READ_RAW_FIELD_H
#define DIRECTIONAL_READ_RAW_FIELD_H
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
  bool IGL_INLINE read_raw_field(const std::string &fileName, int& N, Eigen::MatrixXd& rawField)
  {
    try
    {
      std::ifstream f(fileName);
      int numF;
      f>>N;
      f>> numF;
      rawField.conservativeResize(numF, 3*N);
      
      //Can we do better than element-wise reading?
      for (int i=0;i<rawField.rows();i++)
        for (int j=0;j<rawField.cols();j++)
          f>>rawField(i,j);
      
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
