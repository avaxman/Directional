#include <igl/opengl/glfw/Viewer.h>
#include <directional/visualization_schemes.h>
#include <directional/read_raw_field.h>
#include <directional/read_singularities.h>
#include <directional/singularity_spheres.h>
#include <directional/glyph_lines_raw.h>

int N;
Eigen::MatrixXi FMesh, FField, FSings;
Eigen::MatrixXi EV, FE, EF;
Eigen::MatrixXd VMesh, VField, VSings;
Eigen::MatrixXd CMesh, CField, CSings;
Eigen::MatrixXd rawField;
Eigen::VectorXi singVertices, singIndices;
igl::opengl::glfw::Viewer viewer;

Eigen::MatrixXd positiveIndexColors, negativeIndexColors;
Eigen::RowVector3d rawGlyphColor;
bool drawSingularities=false;


void create_meshes()
{
  directional::glyph_lines_raw(VMesh, FMesh, rawField, directional::default_glyph_color(), VField, FField, CField);
  directional::singularity_spheres(VMesh, FMesh, N, singVertices, singIndices, VSings, FSings, CSings);
  
  //triangle mesh
  viewer.data().set_mesh(VMesh, FMesh);
  viewer.data().set_colors(directional::default_mesh_color());
  viewer.data().set_face_based(true);
  viewer.data().show_lines=true;
  
  //Raw field mesh
  viewer.append_mesh();
  viewer.data().set_mesh(VField, FField);
  viewer.data().set_colors(CField);
  viewer.data().set_face_based(true);
  viewer.data().show_lines=false;
  
  //Singularities mesh
  viewer.append_mesh();
  viewer.data().set_mesh(VSings, FSings);
  viewer.data().set_colors(CSings);
  viewer.data().set_face_based(true);
  viewer.data().show_lines=false;
  viewer.data().show_faces=false;
  
  viewer.selected_data_index=0;  //for all generic libigl commands.
  
}


bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  
  switch (key)
  {
    case GLFW_KEY_SPACE: drawSingularities=!drawSingularities; viewer.data_list[2].show_faces=!viewer.data_list[2].show_faces; /*update_mesh();*/ break;
    default: return false;
  }
  return true;
}



int main()
{
  using namespace Eigen;
  std::cout <<"<space bar>  Show/hide Singularities" << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/rocker-arm2500.obj", VMesh, FMesh);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/rocker-arm2500.rawfield", N, rawField);
  
 
  directional::FEM_suite(VMesh, FMesh, EV, FE, EF, Gv, Ge, J, C, D);
  directional::FEM_masses(VMesh, FMesh, EV, FE, EF, MvVec, MeVec, MfVec, MchiVec);
  
  igl::diag(MvVec,Mv);
  igl::diag(MeVec,Me);
  igl::diag(MfVec,Mf);
  igl::diag(MchiVec,Mchi);
  
  Lv = D*Gv;   //Gv^T * Mchi * Gv
  Le = C*J*Ge; //(JGe)^T * Mchi * JGe
  
  //solving for exact part
   VectorXd scalarFunc;
  igl::min_quad_with_fixed_data<double> mqwfExact;
  // Linear term is 0
  VectorXd B = VectorXd::Zero(V.rows(),1);
  VectorXd Beq;
  SparseMatrix<double> Aeq;
  Eigen::VectorXi b(1); b(0)=0;
  Eigen::VectorXd bc(1); bc(0)=0;
  
  igl::min_quad_with_fixed_precompute(Lv,b,Aeq,true,mqwfExact);
  igl::min_quad_with_fixed_solve(mqwfExact,B,bc,Beq,scalarFunc);
  
  VectorXd exactFunc(V.rows());
  exactFunc<<0.0,scalarFunc;
  
  //solving for coexact part
  igl::min_quad_with_fixed_data<double> mqwfCoexact;
  // Linear term is 0
  VectorXd B = VectorXd::Zero(V.rows(),1);
  VectorXd Beq;
  SparseMatrix<double> Aeq;
  Eigen::VectorXi b(1); b(0)=0;
  Eigen::VectorXd bc(1); bc(0)=0;
  
  igl::min_quad_with_fixed_precompute(Le,b,Aeq,true,mqwfCoexact);
  igl::min_quad_with_fixed_solve(mqwfCoexact,B,bc,Beq,scalarFunc);
  
  VectorXd exactFunc(V.rows());
  exactFunc<<0.0,scalarFunc;

  
  Eigen::VectorXd gradFieldVec = Gv*vScalar;
  Eigen::VectorXd rotCogradFieldVec = J*Ge*eScalar;
  
  gradField.resize(FMesh.rows(),3);
  rotCogradField.resize(FMesh.rows(),3);
  for (int i=0;i<FMesh.rows();i++)
    for (int j=0;j<3;j++){
      gradField(i,j)=gradFieldVec(3*i+j);
      rotCogradField(i,j)=rotCogradFieldVec(3*i+j);
    }
  
  fieldDiv = D*gradFieldVec;
  fieldCurl = C*rotCogradFieldVec;
  
  //visualization meshes
  directional::non_conforming_mesh(VMesh, FMesh, EV, FE, EF, eScalar, VNCMesh, FNCMesh, CNCMesh);
  directional::edge_diamond_mesh(VMesh, FMesh, EV, FE, EF, fieldCurl.cwiseAbs(), VEDMesh, FEDMesh, CEDMesh);
  directional::vertex_area_mesh(VMesh, FMesh, EV, FE, EF, fieldDiv.cwiseAbs(), VVAMesh, FVAMesh, CVAMesh);
  
  std::cout<<"Structure-preserving sanity check: "<<std::endl;
  std::cout<<"(C*gradField).lpNorm<Infinity>(): "<<(C*gradFieldVec).lpNorm<Eigen::Infinity>()<<std::endl;
  std::cout<<"(D*rotCogradField).lpNorm<Infinity>(): "<<(D*rotCogradFieldVec).lpNorm<Eigen::Infinity>()<<std::endl;
  
  create_meshes();
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

