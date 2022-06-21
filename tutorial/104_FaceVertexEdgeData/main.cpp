#include <math.h>
#include <directional/readOFF.h>
#include <directional/directional_viewer.h>
#include <directional/read_raw_field.h>
#include <directional/read_singularities.h>
#include <directional/principal_matching.h>
#include <directional/TriMesh.h>
#include <directional/FaceField.h>

int N;
directional::TriMesh mesh;
directional::FaceField field;
Eigen::VectorXd vertexData, faceData, edgeData;

directional::DirectionalViewer viewer;

bool key_down(igl::opengl::glfw::Viewer& iglViewer, int key, int modifiers)
{
  igl::opengl::glfw::Viewer* iglViewerPointer=&iglViewer;
  directional::DirectionalViewer* directionalViewer = static_cast<directional::DirectionalViewer*>(iglViewerPointer);
  switch (key)
  {
    case '1': directionalViewer->set_face_data(faceData, faceData.minCoeff(),  faceData.maxCoeff()); break;
    case '2': directionalViewer->set_vertex_data(vertexData, vertexData.minCoeff(),  vertexData.maxCoeff()); break;
    case '3': directionalViewer->set_edge_data(edgeData, edgeData.minCoeff(),  edgeData.maxCoeff()); break;
    default: return false;
  }
  directionalViewer->toggle_edge_data(key=='3');
  directionalViewer->toggle_mesh(key!='3');
  return true;
}

int main()
{
  std::cout <<"1  Show Face-based values" << std::endl;
  std::cout <<"2  Show Verex-based values" << std::endl;
  std::cout <<"3  Show Edge-based values" << std::endl;
  
  directional::readOFF(TUTORIAL_SHARED_PATH "/eight.off", mesh);
  field.set_mesh(mesh);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/eight.rawfield", N, field);
 
  //Face data - the x component of the face normals
  faceData=mesh.faceNormals.col(0);
  
  //vertex data: sin(z coordinate)
  vertexData=sin(10.0*mesh.V.col(2).array());
  
  //Edge data - the (squared) effort of the field (under principal matching)
  directional::principal_matching(field);
  edgeData=field.effort.array()*field.effort.array();
  
  viewer.set_mesh(mesh);
  viewer.set_field(field, Eigen::RowVector3d::Constant(1.0));
  viewer.toggle_mesh_edges(false);
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

