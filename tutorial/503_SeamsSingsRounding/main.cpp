#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
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
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField rawField, combedField;
Eigen::MatrixXd NFunctionSeams, NFunctionSings, NCornerFunc;

Eigen::MatrixXd P1Seams, P2Seams, P1Sings, P2Sings;
Eigen::VectorXi funcNumSeams, funcNumSings;

directional::DirectionalViewer viewer;

typedef enum {FIELD, SEAMS_ROUNDING, SINGS_ROUNDING} ViewingModes;
ViewingModes viewingMode=FIELD;

/*void update_viewer()
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
}*/

void callbackFunc(){
    ImGui::PushItemWidth(100);

    const char* items[] = { "Original field", "Rounding seams", "Rounding Singularities"};
    static const char* current_item = NULL;

    if (ImGui::BeginCombo("##combo", current_item)) // The second parameter is the label previewed before opening the combo.
    {
        for (int n = 0; n < IM_ARRAYSIZE(items); n++)
        {
            bool is_selected = (current_item == items[n]); // You can store your selection however you want, outside or inside your objects
            if (ImGui::Selectable(items[n], is_selected)) {
                switch (n) {
                    case 0:
                        viewingMode = FIELD;
                        break;
                    case 1:
                        viewingMode = SEAMS_ROUNDING;
                        viewer.set_isolines(meshCut, NFunctionSeams);
                        break;
                    case 2:
                        viewingMode = SINGS_ROUNDING;
                        viewer.set_isolines(meshCut, NFunctionSings);
                        break;
                }

                current_item = items[n];
            }// You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
            if (is_selected)
                ImGui::SetItemDefaultFocus();
        }
        ImGui::EndCombo();
    }

    ImGui::PopItemWidth();
}


/*bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
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
}*/


int main()
{
  /*std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show integral-seams function" << std::endl <<
  "  3  Show integral-singularity function" << std::endl;*/
  
  directional::readOFF(TUTORIAL_DATA_PATH "/train-station.off", meshWhole);
  ftb.init(meshWhole);
  directional::read_raw_field(TUTORIAL_DATA_PATH "/train-station-5.rawfield", ftb, N, rawField);
  
  //combing and cutting
  directional::principal_matching(rawField);

  directional::IntegrationData intData(N);
  std::cout<<"Setting up Integration"<<std::endl;
  directional::setup_integration(rawField, intData, meshCut, combedField);
  
  intData.verbose=false;
  intData.integralSeamless=true;
  intData.roundSeams=true;
    
  std::cout<<"Seams-rounding Integrating..."<<std::endl;
  directional::integrate(combedField, intData, meshCut,  NFunctionSeams,NCornerFunc);
  std::cout<<"Done!"<<std::endl;
  
  intData.roundSeams=false;
  directional::setup_integration(rawField, intData, meshCut,combedField);
  std::cout<<"Singularity-rounding integration..."<<std::endl;
  directional::integrate(combedField,  intData, meshCut, NFunctionSings,NCornerFunc);
  std::cout<<"Done!"<<std::endl;

  viewer.init();
  viewer.set_mesh(meshWhole,0);
  viewer.set_field(rawField);
  viewer.set_seams(combedField.matching);
  viewer.set_isolines(meshCut, NFunctionSings);
  viewer.set_callback(callbackFunc);
  viewer.launch();
}


