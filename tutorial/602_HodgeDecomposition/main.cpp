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
#include <directional/hodge_decomposition.h>


Eigen::MatrixXi FMesh,  FField, FNCMesh, FVAMesh, FEDMesh;
Eigen::MatrixXi EV, FE, EF;
Eigen::MatrixXd VMesh,  VField, VNCMesh, VVAMesh, VEDMesh;
Eigen::MatrixXd CMesh, CField, CNCMesh, CVAMesh, CEDMesh;
Eigen::MatrixXd rawField, gradField, rotCogradField, harmField;
Eigen::VectorXd exactFunc, coexactFunc;
igl::opengl::glfw::Viewer viewer;

Eigen::SparseMatrix<double> Gv, Ge, J, Mv, Mchi, Mf, Me, C, D;


typedef enum {ORIGINAL_MESH, GRAD_MESH, COGRAD_MESH, HARMONIC_MESH} ViewingModes;
ViewingModes viewingMode=ORIGINAL_MESH;

void update_mesh()
{
  viewer.data_list[0].clear();
  
  //mesh
  switch(viewingMode){
    case ORIGINAL_MESH: viewer.data_list[0].set_mesh(VMesh, FMesh); viewer.data_list[0].set_colors(directional::default_mesh_color()); break;
    case GRAD_MESH: viewer.data_list[0].set_mesh(VMesh, FMesh); viewer.data_list[0].set_colors(exactFunc); break;
    case COGRAD_MESH: viewer.data_list[0].set_mesh(VNCMesh, FNCMesh); viewer.data_list[0].set_colors(CNCMesh); break;
    case HARMONIC_MESH: viewer.data_list[0].set_mesh(VMesh, FMesh); viewer.data_list[0].set_colors(directional::default_mesh_color()); break;
  }
  
  //field
  switch(viewingMode){
    case ORIGINAL_MESH: directional::glyph_lines_raw(VMesh, FMesh,rawField, directional::default_glyph_color(),VField, FField, CField, 3.0); break;
    case GRAD_MESH: directional::glyph_lines_raw(VMesh, FMesh,gradField, directional::default_glyph_color(),VField, FField, CField, 3.0); break;
    case COGRAD_MESH: directional::glyph_lines_raw(VMesh, FMesh,rotCogradField, directional::default_glyph_color(),VField, FField, CField, 3.0); break;
    case HARMONIC_MESH: directional::glyph_lines_raw(VMesh, FMesh,harmField, directional::default_glyph_color(),VField, FField, CField, 3.0); break;
  }

  viewer.data_list[1].clear();
  viewer.data_list[1].set_mesh(VField, FField);
  viewer.data_list[1].set_colors(CField);
  
}


bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  
  switch (key)
  {
    case '1': viewingMode=ORIGINAL_MESH; break;
    case '2': viewingMode=GRAD_MESH; break;
    case '3': viewingMode=COGRAD_MESH; break;
    case '4': viewingMode=HARMONIC_MESH; break;
  }
  update_mesh();
  return true;
}



int main()
{
  using namespace Eigen;
  std::cout <<"1    Original field " << std::endl;
  std::cout <<"2    Exact function + Gradient component " << std::endl;
  std::cout <<"3    Coexact function + Rotated-cogradient component " << std::endl;
  std::cout <<"4    Harmonic component " << std::endl;
  
  int N;
  igl::readOFF(TUTORIAL_SHARED_PATH "/cup_input_simple_10000.off", VMesh, FMesh);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/cup_input_simple_10000.rawfield", N, rawField);
  igl::edge_topology(VMesh, FMesh, EV, FE, EF);
  
  Eigen::Vector3d minV = VMesh.colwise().minCoeff();
  Eigen::Vector3d maxV = VMesh.colwise().maxCoeff();
  Eigen::RowVector3d spanV = maxV-minV;
  
  directional::hodge_decomposition(VMesh, FMesh, EV, FE, EF, rawField, exactFunc, coexactFunc, harmField);
  
  directional::FEM_suite(VMesh, FMesh, EV, FE, EF, Gv, Ge, J, C, D);

  Eigen::VectorXd gradFieldVec = Gv*exactFunc;
  Eigen::VectorXd rotCogradFieldVec = J*Ge*coexactFunc;

  gradField.resize(FMesh.rows(),3);
  rotCogradField.resize(FMesh.rows(),3);
  for (int i=0;i<FMesh.rows();i++)
    for (int j=0;j<3;j++){
      gradField(i,j)=gradFieldVec(3*i+j);
      rotCogradField(i,j)=rotCogradFieldVec(3*i+j);
    }
  
  //visualization meshes
  directional::non_conforming_mesh(VMesh, FMesh, EV, FE, EF, coexactFunc, VNCMesh, FNCMesh, CNCMesh);

  //Triangle mesh
  viewer.data().show_lines = false;
  
  //Raw field
  viewer.append_mesh();
  viewer.data().show_lines = false;
  
  update_mesh();
  viewer.selected_data_index=0;
  
  std::cout<<"Structure-preserving sanity check: "<<std::endl;
  std::cout<<"(C*gradField).lpNorm<Infinity>(): "<<(C*gradFieldVec).lpNorm<Eigen::Infinity>()<<std::endl;
  //std::cout<<"(C*harmField).lpNorm<Infinity>(): "<<(C*harmFieldVec).lpNorm<Eigen::Infinity>()<<std::endl;
  std::cout<<"(D*rotCogradField).lpNorm<Infinity>(): "<<(D*rotCogradFieldVec).lpNorm<Eigen::Infinity>()<<std::endl;
  //std::cout<<"(D*harmField).lpNorm<Infinity>(): "<<(D*harmFieldVec).lpNorm<Eigen::Infinity>()<<std::endl;
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

