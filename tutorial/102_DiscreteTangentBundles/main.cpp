#include <iostream>
#include <fstream>
#include <directional/readOBJ.h>
#include <directional/CartesianField.h>
#include <directional/TriMesh.h>
#include <directional/IntrinsicVertexTangentBundle.h>
#include <directional/PCFaceTangentBundle.h>
#include <directional/read_raw_field.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/directional_viewer.h>
#include <directional/principal_matching.h>
#include "tutorial_shared_path.h"
#include <polyscope/polyscope.h>

using namespace std;

int N = 2;
directional::TriMesh mesh;
directional::IntrinsicVertexTangentBundle vtb;
directional::PCFaceTangentBundle ftb;
directional::CartesianField rawFaceField, powerFaceField;
directional::CartesianField rawVertexField, powerVertexField;
directional::DirectionalViewer viewer;


typedef enum {FACE_FIELD, VERTEX_FIELD} ViewingModes;
ViewingModes viewingMode=FACE_FIELD;


void callbackFunc()
{
    ImGui::PushItemWidth(100); // Make ui elements 100 pixels wide,
    // instead of full width. Must have
    // matching PopItemWidth() below.

    if (viewingMode==FACE_FIELD) {
        if (ImGui::Button("Toggle Vertex Field")) {
            viewingMode=VERTEX_FIELD;
        }
    }

    if (viewingMode==VERTEX_FIELD) {
        if (ImGui::Button("Toggle Face Field")) {
            viewingMode=FACE_FIELD;
        }
    }

    viewer.toggle_singularities(viewingMode==FACE_FIELD, 0);
    viewer.toggle_singularities(viewingMode==VERTEX_FIELD, 1);
    viewer.toggle_cartesian_field(viewingMode==FACE_FIELD,0);
    viewer.toggle_cartesian_field(viewingMode==VERTEX_FIELD,1);

    ImGui::PopItemWidth();
}



int main(int argc, char *argv[])
{
    directional::readOBJ(TUTORIAL_DATA_PATH "/elephant.obj", mesh);
    viewer.init();
    viewer.set_callback(&callbackFunc);

    ftb.init(mesh);
    vtb.init(mesh);

    Eigen::VectorXi constFaces, constVertices;
    Eigen::MatrixXd constVectors;
    constFaces.resize(1);
    constFaces<<0;
    constVectors.resize(1,3);
    constVectors<<mesh.V.row(mesh.F(0,2))-mesh.V.row(mesh.F(0,1));
    constVertices.resize(1);
    constVertices<<mesh.F(0,1);

    directional::power_field(vtb, constVertices, constVectors, Eigen::VectorXd::Constant(constVertices.size(),-1.0), N, powerVertexField);
    directional::power_field(ftb, constFaces, constVectors, Eigen::VectorXd::Constant(constFaces.size(),-1.0), N, powerFaceField);

    //computing power fields
    directional::power_to_raw(powerFaceField, N, rawFaceField,true);
    directional::power_to_raw(powerVertexField, N, rawVertexField,true);
  
    directional::principal_matching(rawFaceField);
    directional::principal_matching(rawVertexField);

    viewer.set_surface_mesh(mesh,0);
    viewer.set_cartesian_field(rawFaceField, "Face-Based Field",  0);
    viewer.set_cartesian_field(rawVertexField, "Vertex-Based Field", 1);
  
    // Update view
    viewer.launch();

    return 0;

}
