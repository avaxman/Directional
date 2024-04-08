#include <iostream>
#include <Eigen/Core>
#include <directional/readOFF.h>
#include <directional/write_raw_field.h>
#include <directional/polyvector_to_raw.h>
#include <directional/raw_to_polyvector.h>
#include <directional/polyvector_field.h>
#include <directional/principal_matching.h>
#include <directional/write_raw_field.h>
#include <directional/directional_viewer.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>

Eigen::VectorXi constFaces;
directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField pvFieldHard, pvFieldSoft, rawFieldHard, rawFieldSoft,constraintsField;
Eigen::MatrixXd constVectors;

double smoothWeight, roSyWeight, alignWeight;

directional::DirectionalViewer viewer;

int N = 4;
typedef enum {CONSTRAINTS, HARD_PRESCRIPTION, SOFT_PRESCRIPTION} ViewingModes;
ViewingModes viewingMode=CONSTRAINTS;


void recompute_field()
{
  directional::polyvector_field(ftb, constFaces, constVectors, smoothWeight, roSyWeight, Eigen::VectorXd::Constant(constFaces.size(),-1), N, pvFieldHard);
  directional::polyvector_field(ftb, constFaces, constVectors, smoothWeight, roSyWeight, alignWeight*Eigen::VectorXd::Constant(constFaces.size(),1.0), N, pvFieldSoft);

  directional::polyvector_to_raw(pvFieldHard, rawFieldHard, N%2==0);
  directional::principal_matching(rawFieldHard);

  directional::polyvector_to_raw(pvFieldSoft, rawFieldSoft, N%2==0);
  directional::principal_matching(rawFieldSoft);

}

void update_visualization()
{
    viewer.set_field(constraintsField,"",0, 0);
    viewer.highlight_faces(constFaces,0);
    if (viewingMode==CONSTRAINTS) {
        viewer.toggle_field(true, 0);
        viewer.toggle_field(false, 1);
    }
    if (viewingMode==HARD_PRESCRIPTION){
        viewer.set_field(rawFieldHard,"", 0, 1);
        viewer.toggle_field(true, 1);
    }

    if (viewingMode==SOFT_PRESCRIPTION){
        viewer.set_field(rawFieldSoft,"", 0, 1);
        viewer.toggle_field(true, 1);
    }

}


void callbackFunc() {
    ImGui::PushItemWidth(100);

    bool recomputeField = false;
    recomputeField = recomputeField || ImGui::InputDouble("RoSy Weight", &roSyWeight);
    recomputeField = recomputeField || ImGui::InputDouble("Alignment Weight", &alignWeight);

    if (recomputeField) {
        recompute_field();
        update_visualization();
    }

    const char* items[] = { "Only Constraints", "Hard prescription", "Soft prescription"};
    static const char* current_item = NULL;

    if (ImGui::BeginCombo("##combo", current_item)) // The second parameter is the label previewed before opening the combo.
    {
        for (int n = 0; n < IM_ARRAYSIZE(items); n++)
        {
            bool is_selected = (current_item == items[n]); // You can store your selection however you want, outside or inside your objects
            if (ImGui::Selectable(items[n], is_selected)){
                switch (n){
                    case 0:
                        viewingMode = CONSTRAINTS;
                        break;
                    case 1:
                        viewingMode = HARD_PRESCRIPTION;
                        break;
                    case 2:
                        viewingMode = SOFT_PRESCRIPTION;
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

    if (ImGui::Button("Save Raw Field"))
        directional::write_raw_field(TUTORIAL_DATA_PATH "/polyvector.rawfield", (viewingMode==HARD_PRESCRIPTION ? rawFieldHard : rawFieldSoft));

    ImGui::PopItemWidth();
}


int main()
{

  // Load mesh
  directional::readOFF(TUTORIAL_DATA_PATH "/fandisk.off", mesh);
  ftb.init(mesh);
  pvFieldHard.init(ftb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);
  pvFieldSoft.init(ftb, directional::fieldTypeEnum::POLYVECTOR_FIELD,N);

  //discovering and constraining sharp edges
  std::vector<int> constFaceslist;
  std::vector<Eigen::Vector3d> constVectorslist;
  for (int i=0;i<mesh.EF.rows();i++){
    if (mesh.faceNormals.row(mesh.EF(i,0)).dot(mesh.faceNormals.row(mesh.EF(i,1)))<0.5){
      constFaceslist.push_back(mesh.EF(i,0));
      constFaceslist.push_back(mesh.EF(i,1));
      constVectorslist.push_back((mesh.V.row(mesh.EV(i,0))-mesh.V.row(mesh.EV(i,1))).normalized());
      constVectorslist.push_back((mesh.V.row(mesh.EV(i,0))-mesh.V.row(mesh.EV(i,1))).normalized());
    }
  }
  
  constFaces.resize(constFaceslist.size());
  constVectors.resize(constVectorslist.size(),3);
  for (int i=0;i<constFaces.size();i++){
    constFaces(i)=constFaceslist[i];
    constVectors.row(i)=constVectorslist[i];
  }
                                
  //generating the viewing fields
  Eigen::MatrixXd rawFieldConstraints=Eigen::MatrixXd::Zero(mesh.F.rows(),N*3);
  Eigen::VectorXi posInFace=Eigen::VectorXi::Zero(mesh.F.rows());
  for (int i=0;i<constFaces.size();i++){
    rawFieldConstraints.block(constFaces(i),3*posInFace(constFaces(i)), 1,3)=constVectors.row(i);
    posInFace(constFaces(i))++;
  }
  
  //Just to show the other direction if N is even, since we are by default constraining it
  if (N%2==0)
    rawFieldConstraints.middleCols(rawFieldConstraints.cols()/2, rawFieldConstraints.cols()/2)=-rawFieldConstraints.middleCols(0, rawFieldConstraints.cols()/2);

  constraintsField.init(ftb, directional::fieldTypeEnum::RAW_FIELD, N);
  constraintsField.set_extrinsic_field(rawFieldConstraints);
  
  smoothWeight = 1.0;
  roSyWeight = 1.0;
  alignWeight = 1.0;
  
  //triangle mesh setup
  viewer.init();
  viewer.set_mesh(mesh,0);
  recompute_field();
  update_visualization();
  viewer.set_callback(callbackFunc);

  viewer.launch();
}
