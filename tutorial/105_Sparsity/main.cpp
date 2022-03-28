#include <math.h>
#include <igl/edge_topology.h>
#include <directional/directional_viewer.h>
#include <directional/read_raw_field.h>
#include <directional/read_singularities.h>
#include <directional/principal_matching.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>

int N=2;
Eigen::MatrixXi F, EV, FE, EF;
Eigen::MatrixXd V, rawField, normals;
Eigen::VectorXi singVertices, singIndices, matching;
Eigen::MatrixXcd powerField;
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
  directionalViewer->set_field(rawField,Eigen::MatrixXd(),0,0.9*((double)sparsity+1.0),sparsity);
  return true;
}

int main()
{
  std::cout <<"1  Sparser field view" << std::endl;
  std::cout <<"2  Denser field view" << std::endl;
 
  igl::readOBJ(TUTORIAL_SHARED_PATH "/armadillo.obj", V, F);
  Eigen::VectorXi bc(1); bc(0)=0;
  Eigen::MatrixXd b(1,3); b.row(0)=V.row(F(0,1))-V.row(F(0,2));
  directional::power_field(V, F, bc,b, Eigen::VectorXd::Constant(bc.size(),-1), N, powerField);
  directional::power_to_raw(V,F,powerField,N,rawField, true);
  
  viewer.set_mesh(V,F);
  viewer.set_field(rawField);
  viewer.set_singularities(singVertices, singIndices);
  viewer.toggle_mesh_edges(false);
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

