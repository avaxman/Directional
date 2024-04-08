#include <iostream>
#include <Eigen/Core>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/combing.h>
#include <directional/directional_viewer.h>
#include <directional/readOBJ.h>

int N;
directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField origField, combedField;
directional::DirectionalViewer viewer;

bool showCombed=false;


void callbackFunc() {
    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,

    if (showCombed) {
        if (ImGui::Button("Show Original Field")) {
            showCombed=false;
            viewer.set_field(origField);
            //viewer.toggle_seams(false);
        }
    } else {
        if (ImGui::Button("Show Combed Field")) {
            showCombed=true;
            viewer.set_field(combedField);
            //viewer.toggle_seams(true);
        }

    }

    ImGui::PopItemWidth();
}



int main()
{

    directional::readOBJ(TUTORIAL_DATA_PATH "/lilium.obj", mesh);
    ftb.init(mesh);
    directional::read_raw_field(TUTORIAL_DATA_PATH "/lilium.rawfield", ftb, N, origField);
    std::cout<<"Showing raw field"<<std::endl;

    //computing the principal matching and then combing the field to trivialize the matching anywhere but a set of seams
    directional::principal_matching(origField);
    directional::combing(origField, combedField);
    Eigen::VectorXd edgeData=origField.effort;
    viewer.init();
    viewer.set_mesh(mesh);
    viewer.set_field(origField);
    showCombed = false;
    viewer.set_seams(combedField.matching);
    //viewer.toggle_seams(false);

    viewer.set_edge_data(edgeData, edgeData.minCoeff(),  edgeData.maxCoeff());

    viewer.set_callback(callbackFunc);
    viewer.launch();
}


