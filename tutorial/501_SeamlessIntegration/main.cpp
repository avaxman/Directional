#include <iostream>
#include <Eigen/Core>
#include <igl/unproject_onto_mesh.h>
#include <igl/cut_mesh.h>
#include <directional/TriMesh.h>
#include <directional/FaceField.h>
#include <directional/readOFF.h>
#include <directional/writeOBJ.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/combing.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/directional_viewer.h>


int N;
directional::TriMesh meshWhole, meshCut;
//Eigen::MatrixXi FMeshWhole, FMeshCut;
//Eigen::MatrixXd VMeshWhole, VMeshCut;
directional::FaceField rawField, combedField;
Eigen::MatrixXd cutUVFull, cutUVRot, cornerWholeUV;
directional::DirectionalViewer viewer;

typedef enum {FIELD, ROT_INTEGRATION, FULL_INTEGRATION} ViewingModes;
ViewingModes viewingMode=FIELD;

//texture image
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> texture_R, texture_G, texture_B;

// Create a texture that hides the integer translation in the parametrization
void setup_line_texture()
{
  unsigned size = 128;
  unsigned size2 = size/2;
  unsigned lineWidth = 5;
  texture_B.setConstant(size, size, 0);
  texture_G.setConstant(size, size, 0);
  texture_R.setConstant(size, size, 0);
  for (unsigned i=0; i<size; ++i)
    for (unsigned j=size2-lineWidth; j<=size2+lineWidth; ++j)
      texture_B(i,j) = texture_G(i,j) = texture_R(i,j) = 255;
  for (unsigned i=size2-lineWidth; i<=size2+lineWidth; ++i)
    for (unsigned j=0; j<size; ++j)
      texture_B(i,j) = texture_G(i,j) = texture_R(i,j) = 255;
}

void update_viewer()
{
  if (viewingMode==FIELD){
    viewer.set_active(true,0);
    viewer.set_active(false,1);
    viewer.toggle_texture(false,1);
  } else if ((viewingMode==ROT_INTEGRATION) || (viewingMode==FULL_INTEGRATION)){
    viewer.set_uv(viewingMode==ROT_INTEGRATION ? cutUVRot : cutUVFull,1);
    viewer.set_active(true,1);
    viewer.set_active(false,0);
    viewer.toggle_texture(true,1);
  }
}

// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = FIELD; break;
    case '2': viewingMode = ROT_INTEGRATION; break;
    case '3': viewingMode = FULL_INTEGRATION; break;
    case 'W':
      Eigen::MatrixXd emptyMat;
      directional::writeOBJ(TUTORIAL_SHARED_PATH "/horsers-param-rot-seamless.obj", meshCut, cutUVRot, meshCut.F);
      directional::writeOBJ(TUTORIAL_SHARED_PATH "/horsers-param-full-seamless.obj", meshCut, cutUVFull, meshCut.F);
      break;
  }
  update_viewer();
  return true;
}


int main()
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show textured rotationally-seamless parameterization mesh" << std::endl <<
  "  3  Show textured fully-seamless parameterization mesh" << std::endl <<
  "  W  Save parameterized OBJ file "<< std::endl;
  
  setup_line_texture();
  directional::readOFF(TUTORIAL_SHARED_PATH "/horsers.off", meshWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/horsers-cf.rawfield", meshWhole, N, rawField);
  
  //combing and cutting
  Eigen::VectorXd curlNorm;
  directional::curl_matching(rawField, curlNorm);
  std::cout<<"curlNorm max: "<<curlNorm.maxCoeff()<<std::endl;

  directional::IntegrationData intData(N);
  std::cout<<"Setting up Integration"<<std::endl;
  directional::setup_integration(meshWhole, rawField, intData, meshCut, combedField);
  
  intData.verbose=true;
  intData.integralSeamless=false;
  
  std::cout<<"Solving for permutationally-seamless integration"<<std::endl;
  directional::integrate(meshWhole, combedField, intData, meshCut, cutUVRot ,cornerWholeUV);
  //Extracting the UV from [U,V,-U, -V];
  cutUVRot=cutUVRot.block(0,0,cutUVRot.rows(),2);
  std::cout<<"Done!"<<std::endl;
  
  intData.integralSeamless = true;  //do not do translational seamless.
  std::cout<<"Solving for integrally-seamless integration"<<std::endl;
  directional::integrate(meshWhole, combedField,  intData, meshCut, cutUVFull,cornerWholeUV);
  cutUVFull=cutUVFull.block(0,0,cutUVFull.rows(),2);
  std::cout<<"Done!"<<std::endl;
  
  //viewer cut (texture) and whole (field) meshes
  viewer.set_mesh(meshWhole,0);
  viewer.set_mesh(meshCut,1);
  viewer.set_field(combedField,directional::DirectionalViewer::indexed_glyph_colors(combedField.extField, false));
  viewer.set_seams(combedField.matching);
  viewer.set_texture(texture_R,texture_G,texture_B,1);
  
  viewer.toggle_texture(false,0);
  viewer.toggle_field(true,0);
  viewer.toggle_seams(true,0);
  
  viewer.toggle_texture(true,1);
  viewer.toggle_field(false,1);
  viewer.toggle_seams(false,1);
  
  update_viewer();
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


