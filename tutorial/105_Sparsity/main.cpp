#include <directional/FaceField.h>
#include <directional/TriMesh.h>
#include <directional/directional_viewer.h>
#include <directional/read_raw_field.h>
#include <directional/read_singularities.h>
#include <directional/principal_matching.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/readOBJ.h>

int N=2;
directional::FaceField field, powerField;
directional::TriMesh mesh;
int sparsity=0;

directional::DirectionalViewer viewer;

bool key_down(igl::opengl::glfw::Viewer& iglViewer, int key, int modifiers)
{
  igl::opengl::glfw::Viewer* iglViewerPointer=&iglViewer;
  directional::DirectionalViewer* directionalViewer = static_cast<directional::DirectionalViewer*>(iglViewerPointer);
  switch (key)
  {
    case '1': sparsity++; break;
    case '2': sparsity--; break;
    default: return false;
  }
  if (sparsity<0) sparsity = 0;
  std::cout<<"Sparsity: "<<sparsity<<std::endl;
  directionalViewer->set_field(field,Eigen::MatrixXd(),0,0.9*((double)sparsity+1.0),sparsity);
  return true;
}

int main()
{
  std::cout <<"1  Sparser field view" << std::endl;
  std::cout <<"2  Denser field view" << std::endl;
 
  directional::readOBJ(TUTORIAL_SHARED_PATH "/armadillo.obj",mesh);
  Eigen::VectorXi bc(1); bc(0)=0;
  Eigen::MatrixXd b(1,3); b.row(0)=mesh.V.row(mesh.F(0,1))-mesh.V.row(mesh.F(0,2));
  directional::power_field(mesh, bc,b, Eigen::VectorXd::Constant(bc.size(),-1), N, powerField);
  directional::power_to_raw(powerField,N,field, true);
  
  viewer.set_mesh(mesh);
  viewer.set_field(field);
  viewer.toggle_mesh_edges(false);
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

