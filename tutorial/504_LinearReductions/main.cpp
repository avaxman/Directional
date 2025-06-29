#include <iostream>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/PCFaceTangentBundle.h>
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
directional::PCFaceTangentBundle ftb;
directional::CartesianField rawField, combedField;
Eigen::MatrixXd NFunctionSign, NFunctionTri, NCornerFunc;
directional::DirectionalViewer viewer;

Eigen::MatrixXd P1Sign, P2Sign, P1Tri, P2Tri;
Eigen::VectorXi funcNumSign, funcNumTri;

typedef enum {SIGN_SYMMETRY, TRI_SYMMETRY} ViewingModes;
ViewingModes viewingMode=SIGN_SYMMETRY;


void callbackFunc(){
    
    const char* items[] = {"Sign Symmetry", "Triangular Symmetry"};
    static const char* current_item = NULL;

    ImGui::PushItemWidth(300);
    static float combo_width = 0.0f;
    if (combo_width == 0.0f) {
        ImGuiStyle& style = ImGui::GetStyle();
        for (auto& item : items)
            combo_width = std::max(combo_width, ImGui::CalcTextSize(item).x);
        combo_width += style.FramePadding.x * 5.0 + ImGui::GetFontSize() + style.ItemInnerSpacing.x;
    }
    
    ImGui::PushItemWidth(combo_width);
    if (ImGui::BeginCombo("Viewing Mode", current_item)) // The second parameter is the label previewed before opening the combo.
    {
        for (int n = 0; n < IM_ARRAYSIZE(items); n++)
        {
            bool is_selected = (current_item == items[n]); // You can store your selection however you want, outside or inside your objects
            if (ImGui::Selectable(items[n], is_selected)) {
                switch (n) {
                    case 0:
                        viewingMode = SIGN_SYMMETRY;
                        viewer.set_isolines(meshCut, NFunctionSign);
                        break;
                    case 1:
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
    viewer.set_surface_mesh(meshWhole,0);
    viewer.set_cartesian_field(rawField);
    //viewer.set_seams(combedField.matching);
    viewer.set_isolines(meshCut, NFunctionSign);
    viewer.set_callback(callbackFunc);
    viewer.launch();

}


