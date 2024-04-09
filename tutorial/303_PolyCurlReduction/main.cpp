#include <iostream>
#include <fstream>
#include <directional/readOFF.h>
#include <directional/readOBJ.h>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/read_raw_field.h>
#include <directional/curl_matching.h>
#include <directional/combing.h>
#include <directional/polycurl_reduction.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>

using namespace std;

directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField rawFieldOrig, rawFieldCF;
directional::CartesianField combedFieldOrig, combedFieldCF;
Eigen::VectorXd curlOrig, curlCF; // norm of curl per edge
directional::DirectionalViewer viewer;

double curlMax, curlMaxOrig;
int N;
int iter=0;

// The set of parameters for calculating the curl-free fields
directional::polycurl_reduction_parameters params;

// Solver data (needed for precomputation)
directional::PolyCurlReductionSolverData pcrdata;


typedef enum {ORIGINAL_FIELD, ORIGINAL_CURL, OPTIMIZED_FIELD, OPTIMIZED_CURL} ViewingModes;
ViewingModes viewingMode=ORIGINAL_FIELD;



void update_visualization()
{
  using namespace std;
  using namespace Eigen;

  Eigen::VectorXd currCurl = (viewingMode==ORIGINAL_CURL ? curlOrig: curlCF);
  viewer.set_edge_data(currCurl, 0.0,curlMaxOrig);
  //viewer.toggle_edge_data((viewingMode==ORIGINAL_CURL) || (viewingMode==OPTIMIZED_CURL))
}

void callbackFunc() {
    ImGui::PushItemWidth(100);

    if (ImGui::Button("Reduce Curl")){
        for (int bi = 0; bi<5; ++bi)
        {
            directional::polycurl_reduction_solve(pcrdata, params, rawFieldCF, iter ==0);
            iter++;
            params.wSmooth *= params.redFactor_wsmooth;
        }

        Eigen::VectorXi prinIndices;
        directional::curl_matching(rawFieldCF, curlCF);
        directional::combing(rawFieldCF, combedFieldCF);
        directional::curl_matching(combedFieldCF,curlCF);
        viewer.set_field(combedFieldCF,"", 0, 1);
        viewer.set_seams(combedFieldCF.matching, 0,1);
    }
    curlMax= curlCF.maxCoeff();
    ImGui::SameLine();
    ImGui::Text("Maxmimum absolute curl: %ld", curlMax);

    const char* items[] = { "Original field", "Original curl", "Optimized field", "Optimized curl"};
    static const char* current_item = NULL;

    if (ImGui::BeginCombo("##combo", current_item)) // The second parameter is the label previewed before opening the combo.
    {
        for (int n = 0; n < IM_ARRAYSIZE(items); n++)
        {
            bool is_selected = (current_item == items[n]); // You can store your selection however you want, outside or inside your objects
            if (ImGui::Selectable(items[n], is_selected)){
                switch (n){
                    case 0:
                        viewingMode = ORIGINAL_FIELD;
                        break;
                    case 1:
                        viewingMode = ORIGINAL_CURL;

                        break;
                    case 2:
                        viewingMode = OPTIMIZED_FIELD;
                        break;
                    case 3:
                        viewingMode = OPTIMIZED_CURL;
                        break;
                }
                update_visualization();
            }
            current_item = items[n];
            if (is_selected)
                ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
        }
        ImGui::EndCombo();
    }

    if (ImGui::Button("Save optimized Field"))
        directional::write_raw_field(TUTORIAL_DATA_PATH "/polycurl.rawfield", combedFieldCF);

    ImGui::PopItemWidth();
}


int main(int argc, char *argv[])
{
  

  // Load a mesh
  directional::readOFF(TUTORIAL_DATA_PATH "/cheburashka.off", mesh);
  ftb.init(mesh);
  rawFieldOrig.init(ftb, directional::fieldTypeEnum::RAW_FIELD, N);
  rawFieldCF.init(ftb, directional::fieldTypeEnum::RAW_FIELD, N);
  directional::read_raw_field(TUTORIAL_DATA_PATH "/cheburashka.rawfield", ftb, N, rawFieldOrig);
  
  //combing the field in a way that minimizes curl
  directional::curl_matching(rawFieldOrig,curlOrig);
  curlMaxOrig= curlOrig.maxCoeff();
  curlMax = curlMaxOrig;

  directional::combing(rawFieldOrig, combedFieldOrig);

  //trivial constraints
  Eigen::VectorXi b; b.resize(1); b<<0;
  Eigen::MatrixXd bc; bc.resize(1,6); bc<<rawFieldOrig.extField.row(0).head(6);
  Eigen::VectorXi blevel; blevel.resize(1); b<<1;
  directional::polycurl_reduction_precompute(mesh, b, bc, blevel, rawFieldOrig , pcrdata);
  
  rawFieldCF = rawFieldOrig;
  combedFieldCF = combedFieldOrig;
  curlCF = curlOrig;
  
  //triangle mesh setup
  viewer.init();
  viewer.set_mesh(mesh);
  viewer.set_field(combedFieldOrig, "", 0, 0);
  viewer.set_seams((combedFieldOrig.matching, 0,0);
  viewer.set_field(combedFieldCF, "", 0, 1);
  viewer.set_seams((combedFieldCF.matching, 0,1);
  update_visualization();
  
  viewer.set_callback(callbackFunc);
  viewer.launch();
  
  return 0;
}
