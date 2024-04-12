#include <iostream>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
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
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField rawField, combedField;
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


void callbackFunc(){
    ImGui::PushItemWidth(100);

    const char* items[] = { "Original field", "Rotationally seamless", "Full seamless"};
    static const char* current_item = NULL;

    if (ImGui::BeginCombo("##combo", current_item)) // The second parameter is the label previewed before opening the combo.
    {
        for (int n = 0; n < IM_ARRAYSIZE(items); n++)
        {
            bool is_selected = (current_item == items[n]); // You can store your selection however you want, outside or inside your objects
            if (ImGui::Selectable(items[n], is_selected)){
                switch (n){
                    case 0:
                        viewingMode = FIELD;
                        break;
                    case 1:
                        viewingMode = ROT_INTEGRATION;
                        viewer.set_uv(cutUVRot,1);
                        break;
                    case 2:
                        viewingMode = FULL_INTEGRATION;
                        viewer.set_uv(cutUVFull,1);
                        break;
                }
            }
            current_item = items[n];
            if (is_selected)
                ImGui::SetItemDefaultFocus();   // You may set the initial focus when opening the combo (scrolling + for keyboard navigation support)
        }
        ImGui::EndCombo();
    }

    if (ImGui::Button("Save Textured OBJ File")){
        Eigen::MatrixXd emptyMat;
        directional::writeOBJ(TUTORIAL_DATA_PATH "/horsers-param-rot-seamless.obj", meshCut, cutUVRot, meshCut.F);
        directional::writeOBJ(TUTORIAL_DATA_PATH "/horsers-param-full-seamless.obj", meshCut, cutUVFull, meshCut.F);
    }

    ImGui::PopItemWidth();
}


int main()
{
  std::cout <<
  "  1  Loaded field" << std::endl <<
  "  2  Show textured rotationally-seamless parameterization mesh" << std::endl <<
  "  3  Show textured fully-seamless parameterization mesh" << std::endl <<
  "  W  Save parameterized OBJ file "<< std::endl;
  
  setup_line_texture();
  directional::readOFF(TUTORIAL_DATA_PATH "/horsers.off", meshWhole);
  ftb.init(meshWhole);
  directional::read_raw_field(TUTORIAL_DATA_PATH "/horsers-cf.rawfield", ftb, N, rawField);
  
  //combing and cutting
  Eigen::VectorXd curlNorm;
  directional::curl_matching(rawField, curlNorm);
  std::cout<<"curlNorm max: "<<curlNorm.maxCoeff()<<std::endl;

  //testing cutting

  directional::IntegrationData intData(N);
  cut_mesh_with_singularities(meshWhole, rawField.singLocalCycles, intData.face2cut);
  Eigen::VectorXi seams = Eigen::VectorXi::Zero(meshWhole.EV.rows());
  for (int i=0;i<meshWhole.F.rows();i++)
      for (int j=0;j<3;j++)
          seams(meshWhole.HE(meshWhole.FH(i,j)))=intData.face2cut(i,j);


  std::cout<<"Setting up Integration"<<std::endl;
  directional::setup_integration(rawField, intData, meshCut, combedField);
  
  intData.verbose=true;
  intData.integralSeamless=false;
  
  std::cout<<"Solving for permutationally-seamless integration"<<std::endl;
  //directional::integrate(combedField, intData, meshCut, cutUVRot ,cornerWholeUV);
  //Extracting the UV from [U,V,-U, -V];
  //cutUVRot=cutUVRot.block(0,0,cutUVRot.rows(),2);
  std::cout<<"Done!"<<std::endl;
  
  intData.integralSeamless = true;  //do not do translational seamless.
  std::cout<<"Solving for integrally-seamless integration"<<std::endl;
  //directional::integrate(combedField,  intData, meshCut, cutUVFull,cornerWholeUV);
  //cutUVFull=cutUVFull.block(0,0,cutUVFull.rows(),2);
  std::cout<<"Done!"<<std::endl;
  
  //viewer cut (texture) and whole (field) meshes
  viewer.init();
  viewer.set_mesh(meshWhole,0);
  //viewer.set_mesh(meshCut,1);
  viewer.set_field(combedField,"", 0,0);
  viewer.set_seams(seams, 0);
  //viewer.set_texture(texture_R,texture_G,texture_B,1);
  //update_viewer();
  viewer.set_callback(callbackFunc);
  viewer.launch();
}


