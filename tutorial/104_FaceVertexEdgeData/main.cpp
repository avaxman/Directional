#include <math.h>
#include <directional/readOFF.h>
#include <directional/directional_viewer.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>

int N;
directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField field;
Eigen::VectorXd vertexData, faceData, edgeData;

directional::DirectionalViewer viewer;


int main()
{

  directional::readOFF(TUTORIAL_DATA_PATH "/eight.off", mesh);
  ftb.init(mesh);
  directional::read_raw_field(TUTORIAL_DATA_PATH "/eight.rawfield", ftb, N, field);
 
  //Face data - the x component of the face normals
  faceData=mesh.faceNormals.col(0);
  
  //vertex data: sin(z coordinate)
  vertexData=sin(10.0*mesh.V.col(2).array());
  
  //Edge data - the (squared) effort of the field (under principal matching)
  directional::principal_matching(field);
  edgeData=field.effort.cwiseAbs2();
  viewer.init();
  viewer.set_mesh(mesh);
  viewer.set_field(field);
  viewer.toggle_mesh_edges(false);
  viewer.set_face_data(faceData, faceData.minCoeff(),  faceData.maxCoeff());
  viewer.set_vertex_data(vertexData, vertexData.minCoeff(),  vertexData.maxCoeff());
  viewer.set_edge_data(edgeData, edgeData.minCoeff(),  edgeData.maxCoeff());
  viewer.launch();
}

