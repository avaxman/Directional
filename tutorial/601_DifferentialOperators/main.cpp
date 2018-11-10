#include <igl/opengl/glfw/Viewer.h>
#include <igl/edge_topology.h>
#include <igl/diag.h>
#include <directional/midedge_mesh.h>
#include <directional/FEM_suite.h>
#include <directional/FEM_masses.h>

#include <directional/visualization_schemes.h>
#include <directional/read_raw_field.h>
#include <directional/glyph_lines_raw.h>


Eigen::MatrixXi FMesh, FMidEdge, FField;
Eigen::MatrixXi EV, FE, EF;
Eigen::MatrixXd VMesh, VMidEdge, VField;
Eigen::MatrixXd CMesh, CField;
Eigen::MatrixXd gradField, rotCogradField;
Eigen::VectorXd fieldDiv, fieldCurl, vScalar, eScalar;
igl::opengl::glfw::Viewer viewer;

Eigen::SparseMatrix<double> Gv, Ge, J, Mv, Mchi, Mf, Me, C, D;


typedef enum {GRAD_MESH, DIV_MESH, COGRAD_MESH, CURL_MESH} ViewingModes;
ViewingModes viewingMode=GRAD_MESH;

void update_mesh()
{
  viewer.data_list[0].clear();
  
  if ((viewingMode==GRAD_MESH)||(viewingMode==DIV_MESH)){
    viewer.data_list[0].set_mesh(VMesh, FMesh);
    viewer.data_list[0].set_colors(viewingMode==GRAD_MESH ? vScalar : fieldDiv);
    directional::glyph_lines_raw(VMesh, FMesh, gradField, directional::default_glyph_color(),VField, FField, CField,2.0);
  }else{  //non-conforming mesh
    viewer.data_list[0].set_mesh(VMesh, FMesh);
    //viewer.data_list[0].set_mesh(VMidEdge, FMidEdge);
    //viewer.data_list[0].set_colors(viewingMode==COGRAD_MESH ? eScalar : fieldCurl);
    directional::glyph_lines_raw(VMesh, FMesh, rotCogradField, directional::default_glyph_color(),VField, FField, CField,2.0);
    
  }
  viewer.data_list[1].clear();
  viewer.data_list[1].set_mesh(VField, FField);
  viewer.data_list[1].set_colors(CField);

}


bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  
  switch (key)
  {
    case '1': viewingMode=GRAD_MESH; break;
    case '2': viewingMode=DIV_MESH; break;
    case '3': viewingMode=COGRAD_MESH; break;
    case '4': viewingMode=CURL_MESH; break;
  }
  update_mesh();
  return true;
}

int main()
{
  std::cout <<"1    Gradient field " << std::endl;
  std::cout <<"2    The divergence of the gradient field " << std::endl;
  std::cout <<"3    Rotated cogradient field " << std::endl;
  std::cout <<"4    The curl of the cogradient field " << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/dragon.off", VMesh, FMesh);
  igl::edge_topology(VMesh, FMesh, EV, FE, EF);
  
  directional::midedge_mesh(VMesh, FMesh, EV, FE, EF, VMidEdge, FMidEdge);
 
  //Procedurally created scalar vertex-based function
  Eigen::Vector3d minV = VMesh.colwise().minCoeff();
  Eigen::Vector3d maxV = VMesh.colwise().maxCoeff();
  Eigen::RowVector3d spanV = maxV-minV;
  
  Eigen::VectorXd MvVec, MeVec, MfVec, MchiVec;
  
  directional::FEM_suite(VMesh, FMesh, EV, FE, EF, Gv, Ge, J, C, D);
  directional::FEM_masses(VMesh, FMesh, EV, FE, EF, MvVec, MeVec, MfVec, MchiVec);
  
  igl::diag(MvVec,Mv);
  igl::diag(MeVec,Me);
  igl::diag(MfVec,Mf);
  igl::diag(MchiVec,Mchi);
  
  vScalar = VMesh.col(0);//.cwiseProduct(VMesh.col(1));
  eScalar = VMidEdge.col(0);//.cwiseProduct(VMidEdge.col(1));
  
  Eigen::VectorXd gradFieldVec = Gv*vScalar;
  Eigen::VectorXd rotCogradFieldVec = J*Ge*eScalar;
  
  /*std::cout<<"Gv*VectorXd::Ones(V.rows()): "<<Gv*Eigen::VectorXd::Ones(VMesh.rows())<<std::endl;
  std::cout<<"Gv*V: "<<Gv*VMesh<<std::endl;*/
  
  gradField.resize(FMesh.rows(),3);
  rotCogradField.resize(FMesh.rows(),3);
  for (int i=0;i<FMesh.rows();i++)
    for (int j=0;j<3;j++){
      gradField(i,j)=gradFieldVec(3*i+j);
      rotCogradField(i,j)=rotCogradFieldVec(3*i+j);
    }
  
  //std::cout<<"gradField: "<<gradField<<std::endl;
  //std::cout<<"rotCogradField: "<<rotCogradField<<std::endl;
  
  fieldDiv = D*gradFieldVec;
  fieldCurl = C*rotCogradFieldVec;
  
  std::cout<<"(C*gradField).lpNorm<Infinity>(): "<<(C*gradFieldVec).lpNorm<Eigen::Infinity>()<<std::endl;
  std::cout<<"(D*rotCogradField).lpNorm<Infinity>(): "<<(D*rotCogradFieldVec).lpNorm<Eigen::Infinity>()<<std::endl;
  
 
  //Triangle mesh
  viewer.data().show_lines = false;

  //Raw field
  viewer.append_mesh();
  viewer.data().show_lines = false;
  
  update_mesh();
  viewer.selected_data_index=0;
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

