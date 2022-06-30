#include <iostream>
#include <fstream>
#include <unordered_set>
#include <igl/readOBJ.h>
#include <directional/readOBJ.h>
#include <directional/TriMesh.h>
#include <directional/FaceField.h>
#include <directional/VertexField.h>
#include <directional/read_raw_field.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/directional_viewer.h>
#include <directional/cut_mesh_with_singularities.h>
#include "tutorial_shared_path.h"



using namespace std;

int N = 2;
directional::TriMesh mesh;
directional::FaceField rawFaceField, powerFaceField;
directional::VertexField rawVertexField, powerVertexField;
directional::DirectionalViewer viewer;


typedef enum {FACE_FIELD, VERTEX_FIELD} ViewingModes;
ViewingModes viewingMode=FACE_FIELD;


void update_viewer()
{
  viewer.toggle_mesh(viewingMode==FACE_FIELD, 0);
  viewer.toggle_mesh(viewingMode==VERTEX_FIELD, 1);
  viewer.toggle_singularities(viewingMode==FACE_FIELD, 0);
  viewer.toggle_singularities(viewingMode==VERTEX_FIELD, 1);
  viewer.toggle_field(viewingMode==FACE_FIELD,0);
  viewer.toggle_field(viewingMode==VERTEX_FIELD,1);
}

// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': viewingMode = FACE_FIELD; break;
    case '2': viewingMode = VERTEX_FIELD; break;
  }
  update_viewer();
  return true;
}



int main(int argc, char *argv[])
{
  auto keyAction = [](const std::string& key, const std::string& description)
  {
    std::cout << "  " << key << "      " << description << std::endl;
  };
  //keyAction("A", "Optimize 60 batches for curl reduction.");
  keyAction("1", "Show face-based field.");
  keyAction("2", "Show vertex-based field.");

  directional::readOBJ(TUTORIAL_SHARED_PATH "/elephant.obj", mesh);
  powerFaceField.init_field(mesh, POWER_FIELD, N);
  powerVertexField.init_field(mesh, POWER_FIELD, N);
  //cout<<"powerVertexField.stiffnessWeights:"<<powerVertexField.stiffnessWeights<<endl;
  //cout<<"powerVertexField.connection:"<<powerVertexField.connection<<endl;
  
  Eigen::VectorXi constFaces, constVertices;
  Eigen::MatrixXd constVectors;
  constFaces.resize(1);
  constFaces<<0;
  constVectors.resize(1,3);
  constVectors<<mesh.V.row(mesh.F(0,2))-mesh.V.row(mesh.F(0,1));
  constVertices.resize(1);
  constVertices<<mesh.F(0,1);
  
  directional::power_field(powerFaceField, constFaces, constVectors, Eigen::VectorXd::Constant(constFaces.size(),-1.0), N);
  directional::power_field(powerVertexField, constVertices, constVectors, Eigen::VectorXd::Constant(constVertices.size(),-1.0), N);

  //computing power fields
  directional::power_to_raw(powerFaceField, N, rawFaceField,true);
  directional::power_to_raw(powerVertexField, N, rawVertexField,true);
  
  directional::principal_matching(rawFaceField);
  directional::principal_matching(rawVertexField);
 
  viewer.set_mesh(mesh,0);
  viewer.set_mesh(mesh,1);
  
  viewer.set_field(rawFaceField,Eigen::MatrixXd(),0, 1.5, 0, 0.3);
  viewer.set_field(rawVertexField,Eigen::MatrixXd(),1, 1.5, 0, 0.4);
  
  // Update view
  update_viewer();
  viewer.callback_key_down = &key_down;
  viewer.launch();
  
  return 0;
}
