#include <igl/edge_topology.h>
#include <directional/readOFF.h>
#include <directional/streamlines.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/directional_viewer.h>
#include <directional/TriMesh.h>
#include <directional/FaceField.h>

// Mesh
Eigen::MatrixXcd powerField;
directional::TriMesh mesh;
directional::FaceField field;

int N=3;         // degree of the vector field
int anim_t = 0;
int anim_t_dir = 1;

bool pre_draw(igl::opengl::glfw::Viewer &iglViewer)
{
  igl::opengl::glfw::Viewer* iglViewerPointer=&iglViewer;
  directional::DirectionalViewer* directional_viewer = static_cast<directional::DirectionalViewer*>(iglViewerPointer);
  using namespace Eigen;
  using namespace std;
  
  if (!iglViewer.core().is_animating)
    return false;
  
  directional_viewer->advance_streamlines();
  
  anim_t += anim_t_dir;
  
  return false;
}

bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int modifier)
{
  if (key == ' ')
  {
    viewer.core().is_animating = !viewer.core().is_animating;
    return true;
  }
  return false;
}

int main(int argc, char *argv[])
{
  using namespace Eigen;
  using namespace std;
  
  directional::DirectionalViewer viewer;
  
  directional::readOFF(TUTORIAL_SHARED_PATH "/lion.off", mesh);
  // Create a Vector Field
  Eigen::VectorXi constFaces(1); constFaces(0) = 0;
  Eigen::MatrixXd constVectors(1, 3); constVectors.row(0) <<(mesh.V.row(mesh.F(0, 1)) - mesh.V.row(mesh.F(0, 0))).normalized();
  Eigen::VectorXd alignWeights(1); alignWeights(0) = -1.0;
  directional::power_field(mesh, constFaces, constVectors, alignWeights ,N, powerField);
  
  // Convert it to raw field
  directional::power_to_raw(powerField,N,field, true);
  
  //triangle mesh
  viewer.set_mesh(mesh);
  viewer.set_field(field);
  Eigen::MatrixXd fieldColors=directional::DirectionalViewer::indexed_glyph_colors(field.extField);
  viewer.set_field_colors(fieldColors);
  viewer.toggle_field(false);
  viewer.init_streamlines();
  viewer.advance_streamlines();  //to get the initial step

  // Viewer Settings
  viewer.callback_pre_draw = &pre_draw;
  viewer.callback_key_down = &key_down;
  viewer.core().is_animating = false;
  viewer.core().animation_max_fps = 30.;
  
  cout << "Press [space] to toggle animation" << endl;
  viewer.launch();
}
