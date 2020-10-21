#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cut_mesh.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/curl_matching.h>
#include <directional/combing.h>
#include <directional/setup_parameterization.h>
#include <directional/parameterize.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/directional_viewer.h>


int N;
Eigen::MatrixXi FWhole, FCut;
Eigen::MatrixXd VWhole, VCut;
Eigen::MatrixXd CField;
Eigen::MatrixXd rawField, combedField;
Eigen::VectorXd effort, combedEffort;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi singIndices, singVertices;
Eigen::MatrixXd cutReducedUV, cutFullUV, cornerWholeUV;

directional::DirectionalViewer viewer;

typedef enum {FIELD, PARAMETERIZATION} ViewingModes;
ViewingModes viewingMode=FIELD;

void update_triangle_mesh()
{
  viewer.set_active(viewingMode==FIELD,0);
  viewer.set_active(!(viewingMode==FIELD),1);
}

// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = FIELD; break;
    case '2': viewingMode = PARAMETERIZATION; break;
    case 'W':
      Eigen::MatrixXd emptyMat;
      igl::writeOBJ(TUTORIAL_SHARED_PATH "/botijo-param.obj", VCut, FCut, emptyMat, emptyMat, cutReducedUV, FCut);
      break;
  }
  update_triangle_mesh();
  return true;
}


int main()
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show textured mesh" << std::endl <<
  "  W  Save parameterized OBJ file "<< std::endl;
  
  igl::readOBJ(TUTORIAL_SHARED_PATH "/botijo.obj", VWhole, FWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/botijo.rawfield", N, rawField);
  igl::edge_topology(VWhole, FWhole, EV, FE, EF);
  
  //combing and cutting
  directional::principal_matching(VWhole, FWhole,EV, EF, FE, rawField, matching, effort,singVertices, singIndices);
  directional::ParameterizationData pd;
  directional::cut_mesh_with_singularities(VWhole, FWhole, singVertices, pd.face2cut);
  directional::combing(VWhole,FWhole, EV, EF, FE, pd.face2cut, rawField, matching, combedField, combedMatching);
   //directional::principal_matching(VWhole, FWhole,EV, EF, FE, combedField,  combedMatching, combedEffort);
 
  std::cout<<"Setting up parameterization"<<std::endl;
  directional::setup_parameterization(directional::sign_symmetry(N), VWhole, FWhole,  EV, EF, FE, combedMatching, singVertices, pd, VCut, FCut);
  
  double lengthRatio=0.01;
  bool isInteger = false;  //do not do translational seamless.
  std::cout<<"Solving parameterization"<<std::endl;
  directional::parameterize(VWhole, FWhole, FE, combedField, lengthRatio, pd, VCut, FCut, isInteger, cutReducedUV,  cutFullUV,cornerWholeUV);
  cutFullUV=cutFullUV.block(0,0,cutFullUV.rows(),2);
  std::cout<<"Done!"<<std::endl;
  
  viewer.set_mesh(VWhole, FWhole,directional::default_mesh_color(),0);
  viewer.set_field(combedField, directional::indexed_glyph_colors(combedField),0);
  viewer.set_singularities(N, singVertices, singIndices,0);
  viewer.set_seams(EV,combedMatching,0);
  viewer.toggle_mesh_edges(false,0);
  viewer.set_active(true,0);
  
  viewer.set_mesh(VCut, FCut,directional::default_mesh_color(),1);
  viewer.set_texture(directional::DirectionalViewer::default_texture(),1);
  viewer.set_uv(cutFullUV,1);
  viewer.toggle_mesh_edges(false,1);
  viewer.set_active(false,1);
    
  update_triangle_mesh();
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


