#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOFF.h>
#include <igl/edge_topology.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/euler_characteristic.h>
#include <igl/readDMAT.h>
#include <igl/writeDMAT.h>
#include <directional/visualization_schemes.h>
#include <directional/dual_cycles.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/index_prescription.h>
#include <directional/rotation_to_representative.h>
#include <directional/representative_to_raw.h>
#include <directional/power_to_representative.h>
#include <directional/power_field.h>
#include <directional/singularity_spheres.h>
#include <directional/glyph_lines_raw.h>

Eigen::VectorXi singVertices,singIndices;
Eigen::VectorXi prinSingVertices,prinSingIndices;

Eigen::SparseMatrix<double> basisCycles;
Eigen::MatrixXi FMesh, FField, FSings;
Eigen::MatrixXd VMesh, VField, VSings;
Eigen::MatrixXi EV, FE, EF;
Eigen::MatrixXd barycenters, faceNormals;
Eigen::VectorXi matching;
Eigen::MatrixXd rawField;
Eigen::VectorXd effort;
Eigen::VectorXi constFaces;
Eigen::MatrixXd constVecMat;
Eigen::VectorXd cycleCurvature;
Eigen::VectorXi vertex2cycle;
Eigen::VectorXi innerEdges;

int N=2;
double globalRotation=0.0;

typedef enum {TRIVIAL_ONE_SING, TRIVIAL_PRINCIPAL_MATCHING, IMPLICIT_FIELD} ViewingModes;
ViewingModes viewingMode=TRIVIAL_ONE_SING;

igl::opengl::glfw::Viewer viewer;


void update_directional_field()
{
  
  using namespace Eigen;
  using namespace std;
  VectorXd rotationAngles;
  VectorXi prinSingIndices=VectorXi::Zero(basisCycles.rows());
  for (int i=0;i<singVertices.size();i++)
    prinSingIndices(singVertices[i])=singIndices[i];
  
  double TCError;
  Eigen::VectorXi currIndices;
  directional::index_prescription(VMesh,FMesh,innerEdges, basisCycles,cycleCurvature,prinSingIndices,N,rotationAngles, TCError);
  
  Eigen::MatrixXd representative;
  directional::rotation_to_representative(VMesh, FMesh,EV,EF,rotationAngles,N,globalRotation, representative);
  directional::representative_to_raw(VMesh,FMesh,representative,N, rawField);
  
  if (viewingMode==TRIVIAL_PRINCIPAL_MATCHING){
    Eigen::VectorXd effort;
    directional::principal_matching(VMesh, FMesh,EV, EF, FE, rawField, matching, effort);
    directional::effort_to_indices(VMesh,FMesh,EV, EF, effort,matching,N,prinSingVertices, prinSingIndices);
  }
  
  if (viewingMode==IMPLICIT_FIELD){
    constVecMat.conservativeResize(constFaces.rows(),3);
    for (int i=0;i<constFaces.size();i++)
      constVecMat.row(i)<<rawField.block(constFaces(i),0,1,3).normalized();
    
    Eigen::VectorXd effort;
    Eigen::MatrixXcd powerField;
    directional::power_field(VMesh, FMesh, constFaces, constVecMat, N, powerField);
    directional::power_to_representative(VMesh,FMesh, powerField,N,representative);
    representative.rowwise().normalize();
    directional::representative_to_raw(VMesh,FMesh,representative,N, rawField);
    directional::principal_matching(VMesh, FMesh,EV,EF,FE,rawField, matching, effort);
    directional::effort_to_indices(VMesh,FMesh,EV,EF,effort,matching,N,prinSingVertices, prinSingIndices);
  }
  
  Eigen::MatrixXd CField, CSings;
  directional::glyph_lines_raw(VMesh, FMesh, rawField, directional::default_glyph_color(), VField, FField, CField);
  
  if (viewingMode==TRIVIAL_ONE_SING)
    directional::singularity_spheres(VMesh, FMesh, N, singVertices, singIndices,  VSings, FSings, CSings);
  
  if ((viewingMode==TRIVIAL_PRINCIPAL_MATCHING)||(viewingMode==IMPLICIT_FIELD))
    directional::singularity_spheres(VMesh, FMesh, N, prinSingVertices, prinSingIndices,  VSings, FSings, CSings);
                                     
  
  viewer.data_list[1].clear();
  viewer.data_list[1].set_mesh(VField, FField);
  viewer.data_list[1].set_colors(CField);
  
  viewer.data_list[2].clear();
  viewer.data_list[2].set_mesh(VSings, FSings);
  viewer.data_list[2].set_colors(CSings);
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
      
    default: break;  //dunno why this is needed but it helps...
      
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
  igl::readOBJ(TUTORIAL_SHARED_PATH "/spherers.obj", VMesh, FMesh);
  igl::edge_topology(VMesh, FMesh, EV,FE,EF);
  igl::barycenter(VMesh,FMesh,barycenters);
  igl::per_face_normals(VMesh,FMesh,faceNormals);
  
  directional::dual_cycles(VMesh, FMesh,EV, EF, basisCycles, cycleCurvature, vertex2cycle, innerEdges);
  
  igl::readDMAT(TUTORIAL_SHARED_PATH "/spheres_constFaces.dmat",constFaces);
  
  singVertices.resize(2);
  singIndices.resize(2);
  singVertices(0)=35;
  singVertices(1)=36;
  singIndices(0)=N;
  singIndices(1)=N;
  
  //triangle mesh
  Eigen::MatrixXd CMesh=directional::default_mesh_color().replicate(FMesh.rows(),1);

  for (int i = 0; i < constFaces.rows(); i++)
    CMesh.row(constFaces(i)) = directional::selected_face_color();
  
  viewer.data().set_mesh(VMesh, FMesh);
  viewer.data().set_colors(CMesh);

  //directional & singularities meshes
  viewer.append_mesh();
  viewer.data().show_lines=false;
  viewer.append_mesh();
  viewer.data().show_lines=false;
  update_directional_field();

  viewer.callback_key_down = &key_down;
  viewer.launch();
}
