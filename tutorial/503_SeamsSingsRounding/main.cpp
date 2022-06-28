#include <iostream>
#include <Eigen/Core>
#include <igl/unproject_onto_mesh.h>
#include <directional/readOFF.h>
#include <directional/TriMesh.h>
#include <directional/FaceField.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/combing.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/branched_isolines.h>
#include <directional/directional_viewer.h>


int N;
directional::TriMesh meshWhole, meshCut;
directional::FaceField rawField, combedField;
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
    viewer.set_isolines(meshCut, NFunctionSeams);
  if (viewingMode==SINGS_ROUNDING)
    viewer.set_isolines(meshCut, NFunctionSings);
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
  
  directional::readOFF(TUTORIAL_SHARED_PATH "/train-station.off", meshWhole);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/train-station-5.rawfield", meshWhole, N, rawField);
  
  //combing and cutting
  directional::principal_matching(rawField);

  directional::IntegrationData intData(N);
  std::cout<<"Setting up Integration"<<std::endl;
  directional::setup_integration(meshWhole, rawField, intData, meshCut, combedField);
  
  intData.verbose=false;
  intData.integralSeamless=true;
  intData.roundSeams=true;
    
  std::cout<<"Seams-rounding Integrating..."<<std::endl;
  directional::integrate(meshWhole, combedField, intData, meshCut,  NFunctionSeams,NCornerFunc);
  std::cout<<"Done!"<<std::endl;
  
  intData.roundSeams=false;
  directional::setup_integration(meshWhole, rawField, intData, meshCut,combedField);
  std::cout<<"Singularity-rounding integration..."<<std::endl;
  directional::integrate(meshWhole, combedField,  intData, meshCut, NFunctionSings,NCornerFunc);
  std::cout<<"Done!"<<std::endl;
  
  viewer.set_mesh(meshWhole,0);
  viewer.set_field(rawField);
  viewer.set_seams(combedField.matching);

  viewer.callback_key_down = &key_down;
  viewer.launch();
}


