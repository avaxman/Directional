#include <igl/readOFF.h>
#include <igl/edge_topology.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/euler_characteristic.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/index_prescription.h>
#include <directional/rotation_to_representative.h>
#include <directional/representative_to_raw.h>
#include <directional/power_to_representative.h>
#include <directional/power_field.h>
#include <directional/directional_viewer.h>


Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::VectorXi singVertices,singIndices;
Eigen::VectorXi prinSingVertices, prinSingIndices;

Eigen::SparseMatrix<double> basisCycles;
Eigen::MatrixXi EV, FE, EF;
Eigen::MatrixXd barycenters, faceNormals;
Eigen::VectorXi matching;
Eigen::MatrixXd rawField;
Eigen::VectorXd effort;
Eigen::VectorXi b;
Eigen::MatrixXd bc;
Eigen::VectorXd cycleCurvature;
Eigen::VectorXi vertex2cycle;
Eigen::VectorXi innerEdges;

int N=2;
double globalRotation=0.0;

typedef enum {TRIVIAL_ONE_SING, TRIVIAL_PRINCIPAL_MATCHING, IMPLICIT_FIELD} ViewingModes;
ViewingModes viewingMode=TRIVIAL_ONE_SING;

directional::DirectionalViewer viewer;


void update_directional_field()
{
  
  using namespace Eigen;
  using namespace std;
  VectorXd rotationAngles;
  prinSingIndices=VectorXi::Zero(basisCycles.rows());
  for (int i=0;i<singVertices.size();i++)
    prinSingIndices(singVertices[i])=singIndices[i];
  
  double IPError;
  Eigen::VectorXi currIndices;
  directional::index_prescription(V,F,innerEdges, basisCycles,cycleCurvature,prinSingIndices,N,rotationAngles, IPError);
  
  Eigen::MatrixXd representative;
  directional::rotation_to_representative(V, F,EV,EF,rotationAngles,N,globalRotation, representative);
  directional::representative_to_raw(V,F,representative,N, rawField);
  
  if (viewingMode==TRIVIAL_PRINCIPAL_MATCHING){
    Eigen::VectorXd effort;
    directional::principal_matching(V, F,EV, EF, FE, rawField, matching, effort,prinSingVertices, prinSingIndices);
  }
  
  if (viewingMode==IMPLICIT_FIELD){
    bc.conservativeResize(b.rows(),3);
    for (int i=0;i<b.size();i++)
      bc.row(i)<<rawField.block(b(i),0,1,3).normalized();
    
    Eigen::VectorXd effort;
    Eigen::MatrixXcd powerField;
    directional::power_field(V, F, b, bc, N, powerField);
    directional::power_to_representative(V,F, powerField,N,representative);
    representative.rowwise().normalize();
    directional::representative_to_raw(V,F,representative,N, rawField);
    directional::principal_matching(V, F,EV,EF,FE,rawField, matching, effort,prinSingVertices, prinSingIndices);
  }
  
  viewer.set_field(rawField);
  
  if (viewingMode==TRIVIAL_ONE_SING)
    viewer.set_singularities(singVertices, singIndices);
   
  if ((viewingMode==TRIVIAL_PRINCIPAL_MATCHING)||(viewingMode==IMPLICIT_FIELD))
    viewer.set_singularities(prinSingVertices, prinSingIndices);
  
}



bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifiers)
{
  using namespace std;
  switch(key)
  {
    case '1': viewingMode=TRIVIAL_ONE_SING; cout<<"Showing prescribed singularity "<<std::endl;
      break;
    case '2': viewingMode=TRIVIAL_PRINCIPAL_MATCHING; cout<<"Principal-matching singularities "<<std::endl;
      break;
    case '3': viewingMode=IMPLICIT_FIELD; cout<<"Field interpolated from constraints with principal singularities "<<std::endl;
      break;
      
    case '4':{
      singIndices[0]--;
      singIndices[1]++;
      cout<<"Prescribed singularity index: "<<singIndices[0]<<"/"<<N<<std::endl;
      break;
    }
    case '5':{
      singIndices[0]++;
      singIndices[1]--;
      cout<<"Prescribed singularity index: "<<singIndices[0]<<"/"<<N<<std::endl;
      break;
    }
      
    case '6':{
      globalRotation+=igl::PI/16;
      std::cout<<"globalRotation: " <<globalRotation<<std::endl;
      break;
    }
    default: break;
  }
  update_directional_field();
  return true;
}

int main()
{
  std::cout <<
  "1-3      Toggle between singularity modes" << std::endl <<
  "4-5      Decrease/increase prescribed singularity index" << std::endl <<
  "6        Change global rotation" << std::endl;
  using namespace Eigen;
  using namespace std;
  igl::readOBJ(TUTORIAL_SHARED_PATH "/spherers.obj", V, F);
  igl::edge_topology(V, F, EV,FE,EF);
  igl::barycenter(V,F,barycenters);
  igl::per_face_normals(V,F,faceNormals);
  
  directional::dual_cycles(V, F,EV, EF, basisCycles, cycleCurvature, vertex2cycle, innerEdges);
  
  igl::readDMAT(TUTORIAL_SHARED_PATH "/spheres_constFaces.dmat",b);
  
  singVertices.resize(2);
  singIndices.resize(2);
  singVertices(0)=35;
  singVertices(1)=36;
  singIndices(0)=N;
  singIndices(1)=N;
  
  //viewing mesh
  viewer.set_mesh(V, F);
  viewer.set_selected_faces(b);
  update_directional_field();
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}
