#include <igl/opengl/glfw/Viewer.h>
#include <igl/edge_topology.h>
#include <igl/diag.h>
#include <directional/non_conforming_mesh.h>
#include <directional/edge_diamond_mesh.h>
#include <directional/vertex_area_mesh.h>
#include <directional/FEM_suite.h>
#include <directional/FEM_masses.h>
#include <igl/min_quad_with_fixed.h>

#include <directional/visualization_schemes.h>
#include <directional/read_raw_field.h>
#include <directional/glyph_lines_raw.h>
#include <directional/harmonic_basis.h>


Eigen::MatrixXi FMesh,  FField;
Eigen::MatrixXi EV, FE, EF;
Eigen::MatrixXd VMesh,  VField;
Eigen::MatrixXd CMesh, CField;
std::vector<Eigen::MatrixXd> harmFields;
Eigen::VectorXd exactFunc, coexactFunc;
igl::opengl::glfw::Viewer viewer;

Eigen::SparseMatrix<double> Gv, Ge, J, Mv, Mchi, Mf, Me, C, D;
int currField=0;



void update_field()
{

  //mesh
  directional::glyph_lines_raw(VMesh, FMesh,harmFields[currField], directional::default_glyph_color(),VField, FField, CField, 25.0);
  viewer.data_list[1].clear();
  viewer.data_list[1].set_mesh(VField, FField);
  viewer.data_list[1].set_colors(CField);
  
}


bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  
  switch (key)
  {
    case '1': currField = (currField+1)%harmFields.size();
  }
  update_field();
  return true;
}



int main()
{
  using namespace Eigen;
  std::cout <<"1    Switch between harmonic basis fields " << std::endl;
  int N;
  igl::readOFF(TUTORIAL_SHARED_PATH "/cup_input_simple_10000.off", VMesh, FMesh);
  igl::edge_topology(VMesh, FMesh, EV, FE, EF);
  
  directional::harmonic_basis(VMesh, FMesh, EV, FE, EF, harmFields);
  
  //Triangle mesh
  viewer.data_list[0].clear();
  viewer.data_list[0].set_mesh(VMesh, FMesh);
  viewer.data_list[0].set_colors(directional::default_mesh_color());
  viewer.data().show_lines = false;
  
  //Raw field
  viewer.append_mesh();
  viewer.data().show_lines = false;
  
  update_field();
  viewer.selected_data_index=0;
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

