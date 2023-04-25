#include <directional/readOBJ.h>
#include <directional/streamlines.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/directional_viewer.h>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>

directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField field, powerField;

int N = 3;
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
  
  directional_viewer->advance_streamlines(0.5);
  
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
  
  directional::readOBJ(TUTORIAL_SHARED_PATH "/spherers.obj", mesh);
  ftb.init(mesh);
  // Create a power field
  Eigen::VectorXi constFaces(1); constFaces(0) = 0;
  Eigen::MatrixXd constVectors(1, 3); constVectors.row(0) <<(mesh.V.row(mesh.F(0, 1)) - mesh.V.row(mesh.F(0, 0))).normalized();
  Eigen::VectorXd alignWeights(1); alignWeights(0) = -1.0;
  directional::power_field(ftb, constFaces, constVectors, alignWeights ,N, powerField);
  directional::power_to_raw(powerField,N,field, true);

  viewer.set_mesh(mesh);
  viewer.set_field(field);
  Eigen::MatrixXd fieldColors=directional::DirectionalViewer::indexed_glyph_colors(field.extField);
  viewer.set_field_colors(fieldColors);
  viewer.toggle_field(false);
  viewer.init_streamlines(0, Eigen::VectorXi(), 3);
  viewer.advance_streamlines(0.5);  //to get the initial step

  // Viewer Settings
  viewer.callback_pre_draw = &pre_draw;
  viewer.callback_key_down = &key_down;
  viewer.core().is_animating = false;
  viewer.core().animation_max_fps = 30.;
  
  cout << "Press [space] to toggle animation" << endl;
  viewer.launch();
}
