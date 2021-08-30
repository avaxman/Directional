#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <igl/cut_mesh.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/combing.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/cut_mesh_with_singularities.h>
#include <directional/branched_isolines.h>
#include <directional/directional_viewer.h>

int N;
Eigen::MatrixXi FMeshWhole, FMeshCut;
Eigen::MatrixXd VMeshWhole, VMeshCut;
Eigen::MatrixXd rawField, combedField;
Eigen::VectorXd effort, combedEffort;
Eigen::VectorXi matching, combedMatching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi singIndices, singVertices;
Eigen::MatrixXd NFunctionSign, NFunctionTri, NCornerFunc;
directional::DirectionalViewer viewer;

Eigen::MatrixXd P1Sign, P2Sign, P1Tri, P2Tri;
Eigen::VectorXi funcNumSign, funcNumTri;

typedef enum {FIELD, SIGN_SYMMETRY, TRI_SYMMETRY} ViewingModes;
ViewingModes viewingMode=FIELD;


void update_viewer()
{
  if (viewingMode==FIELD){
    viewer.toggle_seams(true);
    viewer.toggle_field(true);
    viewer.toggle_singularities(true);
    viewer.toggle_isolines(false);
  } else if ((viewingMode==SIGN_SYMMETRY) || (viewingMode==TRI_SYMMETRY)){
    viewer.toggle_seams(false);
    viewer.toggle_field(false);
    viewer.toggle_singularities(true);
    viewer.toggle_isolines(true);
  }
  
  if (viewingMode==SIGN_SYMMETRY)
    viewer.set_isolines(VMeshCut, FMeshCut, NFunctionSign);
  if (viewingMode==TRI_SYMMETRY)
    viewer.set_isolines(VMeshCut, FMeshCut, NFunctionTri);
}

// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = FIELD; break;
    case '2': viewingMode = SIGN_SYMMETRY; break;
    case '3': viewingMode = TRI_SYMMETRY; break;
  }
  update_viewer();
  return true;
}


int main()
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show only sign-symmetric integrated functions" << std::endl <<
  "  3  Show triangular-symmetric integrated functions" << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/dome.off", VMeshWhole, FMeshWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/dome-6.rawfield", N, rawField);
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);

  //combing and cutting
  directional::principal_matching(VMeshWhole, FMeshWhole,EV, EF, FE, rawField, matching, effort,singVertices, singIndices);

  directional::IntegrationData intData(N);
  std::cout<<"Setting up Integration"<<std::endl;
  directional::setup_integration(VMeshWhole, FMeshWhole,  EV, EF, FE, rawField, matching, singVertices, intData, VMeshCut, FMeshCut, combedField, combedMatching);
  
  intData.verbose=false;
  intData.integralSeamless=true;

  std::cout<<"Free (sign-symmetric) Integrating..."<<std::endl;
  directional::integrate(VMeshWhole, FMeshWhole, FE, combedField, intData, VMeshCut, FMeshCut, NFunctionSign, NCornerFunc);
  std::cout<<"Done!"<<std::endl;
  
  
  std::cout<<"Solving triangular-constrained integration..."<<std::endl;
  intData.set_triangular_symmetry(N);
  directional::setup_integration(VMeshWhole, FMeshWhole,  EV, EF, FE, rawField, matching, singVertices, intData, VMeshCut, FMeshCut, combedField, combedMatching);
  directional::integrate(VMeshWhole, FMeshWhole, FE, combedField,  intData, VMeshCut, FMeshCut,  NFunctionTri, NCornerFunc);
  std::cout<<"Done!"<<std::endl;
  
  viewer.set_mesh(VMeshWhole, FMeshWhole,directional::DirectionalViewer::default_mesh_color(), 0);
  viewer.set_field(rawField);
  viewer.set_singularities(singVertices, singIndices);
  viewer.set_seams(EV, combedMatching);
  
  update_viewer();
 
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


