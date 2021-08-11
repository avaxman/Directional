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
Eigen::MatrixXd NFunctionSeams, NFunctionSings, NCornerFunc;

Eigen::MatrixXd P1Seams, P2Seams, P1Sings, P2Sings;
Eigen::VectorXi funcNumSeams, funcNumSings;

directional::DirectionalViewer viewer;

typedef enum {FIELD, SEAMS_ROUNDING, SINGS_ROUNDING} ViewingModes;
ViewingModes viewingMode=FIELD;

void update_viewer()
{
  if (viewingMode==FIELD){
    viewer.toggle_seams(true);
    viewer.toggle_field(true);
    viewer.toggle_singularities(true);
    viewer.toggle_isolines(false);
  } else if ((viewingMode==SEAMS_ROUNDING) || (viewingMode==SINGS_ROUNDING)){
    viewer.toggle_seams(false);
    viewer.toggle_field(false);
    viewer.toggle_singularities(true);
    viewer.toggle_isolines(true);
  }
  
  if (viewingMode==SEAMS_ROUNDING)
    viewer.set_isolines(P1Seams, P2Seams, funcNumSeams);
  if (viewingMode==SINGS_ROUNDING)
    viewer.set_isolines(P1Sings, P2Sings, funcNumSings);
}


bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = FIELD; break;
    case '2': viewingMode = SEAMS_ROUNDING; break;
    case '3': viewingMode = SINGS_ROUNDING; break;
  }
  update_viewer();
  return true;
}


int main()
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show integral-seams function" << std::endl <<
  "  3  Show integral-singularity function" << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/train-station.off", VMeshWhole, FMeshWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/train-station-5.rawfield", N, rawField);
  igl::edge_topology(VMeshWhole, FMeshWhole, EV, FE, EF);

  //combing and cutting
  directional::principal_matching(VMeshWhole, FMeshWhole,EV, EF, FE, rawField, matching, effort,singVertices, singIndices);

  directional::IntegrationData intData(N);
  std::cout<<"Setting up Integration"<<std::endl;
  directional::setup_integration(VMeshWhole, FMeshWhole,  EV, EF, FE, rawField, matching, singVertices, intData, VMeshCut, FMeshCut, combedField, combedMatching);
  
  intData.verbose=false;
  intData.integralSeamless=true;
  intData.roundSeams=true;
    
  std::cout<<"Seams-rounding Integrating..."<<std::endl;
  directional::integrate(VMeshWhole, FMeshWhole, FE, combedField, intData, VMeshCut, FMeshCut,  NFunctionSeams,NCornerFunc);
  std::cout<<"Done!"<<std::endl;
  
  intData.roundSeams=false;
  directional::setup_integration(VMeshWhole, FMeshWhole,  EV, EF, FE, rawField, matching, singVertices, intData, VMeshCut, FMeshCut, combedField, combedMatching);
  std::cout<<"Singularity-rounding integration..."<<std::endl;
  directional::integrate(VMeshWhole, FMeshWhole, FE, combedField,  intData, VMeshCut, FMeshCut, NFunctionSings,NCornerFunc);
  std::cout<<"Done!"<<std::endl;
  
  viewer.set_mesh(VMeshWhole, FMeshWhole,directional::DirectionalViewer::default_mesh_color(), 0);
  viewer.set_field(rawField);
  viewer.set_singularities(singVertices, singIndices);
  viewer.set_seams(EV, combedMatching);
  
  directional::branched_isolines(VMeshCut, FMeshCut, NFunctionSeams, P1Seams, P2Seams, funcNumSeams);
  directional::branched_isolines(VMeshCut, FMeshCut, NFunctionSings, P1Sings, P2Sings, funcNumSings);

  viewer.callback_key_down = &key_down;
  viewer.launch();
}


