#include <igl/edge_topology.h>
#include <igl/read_triangle_mesh.h>
#include <igl/readOFF.h>
#include <directional/streamlines.h>
#include <directional/power_field.h>
#include <directional/power_to_raw.h>
#include <directional/directional_viewer.h>


// Mesh
Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXcd powerField;
Eigen::MatrixXd rawField;
Eigen::MatrixXd P1,P2;

/*directional::StreamlineData sl_data;
directional::StreamlineState sl_state;*/

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
  
  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/lion.off", V, F);
  // Create a Vector Field
  directional::power_field(V, F, Eigen::VectorXi(),  Eigen::MatrixXd() , Eigen::VectorXd() ,3, powerField);
  
  // Convert it to raw field
  directional::power_to_raw(V,F,powerField,3,rawField, true);
  
  //triangle mesh
  viewer.set_mesh(V,F);
  viewer.set_field(rawField);
  Eigen::MatrixXd fieldColors=directional::DirectionalViewer::indexed_glyph_colors(rawField);
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
