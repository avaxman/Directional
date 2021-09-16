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
#include <directional/effort_to_indices.h>
#include <directional/combing.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/directional_viewer.h>


int N;
Eigen::MatrixXi FMeshWhole, FMeshCut;
Eigen::MatrixXd VMeshWhole, VMeshCut;
Eigen::MatrixXd rawField, combedField;
Eigen::VectorXd effort, combedEffort;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi singIndices, singVertices;
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
  } else if ((viewingMode==ROT_INTEGRATION) || (viewingMode==FULL_INTEGRATION)){
    viewer.set_uv(viewingMode==ROT_INTEGRATION ? cutUVRot : cutUVFull,1);
    viewer.set_active(true,1);
    viewer.set_active(false,0);
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
      igl::writeOBJ(TUTORIAL_SHARED_PATH "/horsers-param-rot-seamless.obj", VMeshCut, FMeshCut, emptyMat, emptyMat, cutUVRot, FMeshCut);
      igl::writeOBJ(TUTORIAL_SHARED_PATH "/horsers-param-full-seamless.obj", VMeshCut, FMeshCut, emptyMat, emptyMat, cutUVFull, FMeshCut);
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
  igl::readOFF(TUTORIAL_SHARED_PATH "/horsers.off", VMeshWhole, FMeshWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/horsers-cf.rawfield", N, rawField);
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);

  //combing and cutting
  Eigen::VectorXd curlNorm;
  directional::curl_matching(VMeshWhole, FMeshWhole,EV, EF, FE, rawField, matching, effort, curlNorm,singVertices, singIndices);
  std::cout<<"curlNorm max: "<<curlNorm.maxCoeff()<<std::endl;

  directional::IntegrationData intData(N);
  std::cout<<"Setting up Integration"<<std::endl;
  directional::setup_integration(VMeshWhole, FMeshWhole,  EV, EF, FE, rawField, matching, singVertices, intData, VMeshCut, FMeshCut, combedField, combedMatching);
  
  intData.verbose=true;
  intData.integralSeamless=false;
  
  std::cout<<"Solving for permutationally-seamless integration"<<std::endl;
  directional::integrate(VMeshWhole, FMeshWhole, FE, combedField, intData, VMeshCut, FMeshCut, cutUVRot ,cornerWholeUV);
  //Extracting the UV from [U,V,-U, -V];
  cutUVRot=cutUVRot.block(0,0,cutUVRot.rows(),2);
  std::cout<<"Done!"<<std::endl;
  
  intData.integralSeamless = true;  //do not do translational seamless.
  std::cout<<"Solving for integrally-seamless integration"<<std::endl;
  directional::integrate(VMeshWhole, FMeshWhole, FE, combedField,  intData, VMeshCut, FMeshCut, cutUVFull,cornerWholeUV);
  cutUVFull=cutUVFull.block(0,0,cutUVFull.rows(),2);
  std::cout<<"Done!"<<std::endl;
  
  //viewer cut (texture) and whole (field) meshes
  viewer.set_mesh(VMeshWhole, FMeshWhole,0);
  viewer.set_mesh(VMeshCut, FMeshCut,1);
  viewer.set_field(rawField);
  viewer.set_singularities(singVertices, singIndices);
  viewer.set_seams(combedMatching);
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


