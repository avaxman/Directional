#include <math.h>
#include <igl/edge_topology.h>
#include <igl/per_face_normals.h>
#include <directional/directional_viewer.h>
#include <directional/read_raw_field.h>
#include <directional/read_singularities.h>
#include <directional/principal_matching.h>

int N;
Eigen::MatrixXi F, EV, FE, EF;
Eigen::MatrixXd V, rawField, normals;
Eigen::VectorXi singVertices, singIndices, matching;
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
    case '3': directionalViewer->set_edge_data(edgeData, edgeData.minCoeff(),  edgeData.maxCoeff(), EV, FE, EF); break;
    default: return false;
  }
  return true;
}

int main()
{
  std::cout <<"1  Show Face-based values" << std::endl;
  std::cout <<"2  Show Verex-based values" << std::endl;
  std::cout <<"3  Show Edge-based values" << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/eight.off", V, F);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/eight.rawfield", N, rawField);
  igl::edge_topology(V, F, EV, FE, EF);

  //Face data - the x component of normals
  igl::per_face_normals(V, F, normals);
  faceData=normals.col(0);
  
  //vertex data: sin(z coordinate)
  vertexData=sin(10.0*V.col(2).array());
  
  //Edge data - the (squared) effort of the field (under principal matching)
  directional::principal_matching(V, F, EV, EF, FE, rawField, matching, edgeData, singVertices, singIndices);
  edgeData.array()*=edgeData.array();
  
  viewer.set_mesh(V,F);
  viewer.set_field(rawField, Eigen::RowVector3d::Constant(1.0));
  viewer.set_singularities(singVertices, singIndices);
  viewer.toggle_mesh_edges(false);
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

