#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cut_mesh.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
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
Eigen::MatrixXd cutUVFull, cutUVRot, cornerWholeUV, cutReducedUV;
directional::DirectionalViewer viewer;

typedef enum {FIELD, ROT_PARAMETERIZATION, FULL_PARAMETERIZATION, UV_COORDS} ViewingModes;
ViewingModes viewingMode=FIELD;

void update_triangle_mesh()
{
  viewer.set_active(viewingMode==FIELD,0);
  viewer.set_active(!(viewingMode==FIELD),1);
  viewer.set_uv(viewingMode==ROT_PARAMETERIZATION ? cutUVRot : cutUVFull,1);
}

// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = FIELD; break;
    case '2': viewingMode = ROT_PARAMETERIZATION; break;
    case '3': viewingMode = FULL_PARAMETERIZATION; break;
    case 'W':
      Eigen::MatrixXd emptyMat;
      igl::writeOBJ(TUTORIAL_SHARED_PATH "/horsers-param-rot-seamless.obj", VCut, FCut, emptyMat, emptyMat, cutUVRot, FCut);
      igl::writeOBJ(TUTORIAL_SHARED_PATH "/horsers-param-full-seamless.obj", VCut, FCut, emptyMat, emptyMat, cutUVFull, FCut);
      break;
  }
  update_triangle_mesh();
  return true;
}


int main()
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show textured rotationally-seamless parameterization mesh" << std::endl <<
  "  3  Show textured fully-seamless parameterization mesh" << std::endl <<
  "  W  Save parameterized OBJ file "<< std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/horsers.off", VWhole, FWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/horsers-cf.rawfield", N, rawField);
  igl::edge_topology(VWhole, FWhole, EV, FE, EF);
  
  //combing and cutting
  Eigen::VectorXd curlNorm;
  directional::curl_matching(VWhole, FWhole,EV, EF, FE, rawField, matching, effort, curlNorm,singVertices, singIndices);
  
  directional::ParameterizationData pd;
  directional::cut_mesh_with_singularities(VWhole, FWhole, singVertices, pd.face2cut);
  directional::combing(VWhole,FWhole, EV, EF, FE, pd.face2cut, rawField, matching, combedField, combedMatching);
  std::cout<<"curlNorm max: "<<curlNorm.maxCoeff()<<std::endl;
  
  std::cout<<"Setting up parameterization"<<std::endl;
  
  directional::setup_parameterization(directional::sign_symmetry(N), VWhole, FWhole, EV, EF, FE, combedMatching, singVertices, pd, VCut, FCut);
  
  double lengthRatio=0.01;
  bool isInteger = false;  //do not do translational seamless.
  std::cout<<"Solving rotationally-seamless parameterization"<<std::endl;
  directional::parameterize(VWhole, FWhole, FE, combedField, lengthRatio, pd, VCut, FCut, isInteger, cutReducedUV, cutUVRot, cornerWholeUV);
  
  cutUVRot=cutUVRot.block(0,0,cutUVRot.rows(),2);
  std::cout<<"Done!"<<std::endl;
  
  isInteger = true;  //do not do translational seamless.
  std::cout<<"Solving fully-seamless parameterization"<<std::endl;
  directional::parameterize(VWhole, FWhole, FE, combedField, lengthRatio, pd, VCut, FCut, isInteger,   cutReducedUV, cutUVFull, cornerWholeUV);
  cutUVFull=cutUVFull.block(0,0,cutUVFull.rows(),2);
  std::cout<<"Done!"<<std::endl;
  
  viewer.set_mesh(VWhole, FWhole,directional::default_mesh_color(),0);
  viewer.set_field(combedField, directional::indexed_glyph_colors(combedField),0);
  viewer.set_singularities(N, singVertices, singIndices,0);
  viewer.set_seams(EV,combedMatching,0);
  viewer.toggle_mesh_edges(false,0);
  viewer.set_active(true,0);
  
  viewer.set_mesh(VCut, FCut,directional::default_mesh_color(),1);
  viewer.set_texture(directional::DirectionalViewer::default_texture(),1);
  viewer.set_uv(cutUVRot,1);
  viewer.toggle_mesh_edges(false,1);
  viewer.set_active(false,1);
  
  /*Eigen::VectorXi isSeam=Eigen::VectorXi::Zero(EV.rows());
   for (int i=0;i<FE.rows();i++)
   for (int j=0;j<3;j++)
   if (pd.face2cut(i,j))
   isSeam(FE(i,j))=1;
   directional::seam_lines(VWhole, FWhole, EV, combedMatching, VSeams, FSeams, CSeams,2.5);*/
  
  update_triangle_mesh();
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


