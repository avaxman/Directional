#include <iostream>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/readOFF.h>
#include <directional/read_raw_field.h>
#include <directional/write_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/setup_integration.h>
#include <directional/integrate.h>
#include <directional/branched_isolines.h>
#include <directional/directional_viewer.h>

int N;
directional::TriMesh meshWhole, meshCut;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField rawField, combedField;
Eigen::MatrixXd NFunctionSign, NFunctionTri, NCornerFunc;
directional::DirectionalViewer viewer;

Eigen::MatrixXd P1Sign, P2Sign, P1Tri, P2Tri;
Eigen::VectorXi funcNumSign, funcNumTri;

typedef enum {FIELD, SIGN_SYMMETRY, TRI_SYMMETRY} ViewingModes;
ViewingModes viewingMode=FIELD;


/*void update_viewer()
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
    viewer.set_isolines(meshCut, NFunctionSign);
  if (viewingMode==TRI_SYMMETRY)
    viewer.set_isolines(meshCut, NFunctionTri);
}*/

void callbackFunc(){
    ImGui::PushItemWidth(100);

    const char* items[] = { "Original field", "Sign Symmetry", "Triangular Symmetry"};
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
                        viewingMode = SIGN_SYMMETRY;
                        viewer.set_isolines(meshCut, NFunctionSign);
                        break;
                    case 2:
                        viewingMode = TRI_SYMMETRY;
                        viewer.set_isolines(meshCut, NFunctionTri);
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

// Handle keyboard input
/*bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
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
}*/


int main()
{

    directional::readOFF(TUTORIAL_DATA_PATH "/dome.off", meshWhole);
    ftb.init(meshWhole);
    directional::read_raw_field(TUTORIAL_DATA_PATH "/dome-6.rawfield", ftb, N, rawField);

    //combing and cutting
    directional::principal_matching(rawField);

    directional::IntegrationData intData(N);
    std::cout<<"Setting up Integration"<<std::endl;
    directional::setup_integration(rawField, intData,meshCut, combedField);

    intData.verbose=false;
    intData.integralSeamless=true;

    std::cout<<"Free (sign-symmetric) Integrating..."<<std::endl;
    directional::integrate(combedField, intData, meshCut, NFunctionSign, NCornerFunc);
    std::cout<<"Done!"<<std::endl;


    std::cout<<"Solving triangular-constrained integration..."<<std::endl;
    intData.set_triangular_symmetry(N);
    directional::setup_integration(rawField,intData, meshCut, combedField);
    directional::integrate(combedField,  intData, meshCut, NFunctionTri, NCornerFunc);
    std::cout<<"Done!"<<std::endl;

    viewer.init();
    viewer.set_mesh(meshWhole,0);
    viewer.set_field(rawField);
    viewer.set_seams(combedField.matching);
    viewer.set_isolines(meshCut, NFunctionSign);
    viewer.set_callback(callbackFunc);
    viewer.launch();

}


