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

directional::StreamlineData sl_data;
directional::StreamlineState sl_state;

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
  
  directional::streamlines_next(V, F, sl_data, sl_state);
  Eigen::RowVector3d color = Eigen::RowVector3d::Zero();
  double value = ((anim_t) % 100) / 100.;
  
  if (value > 0.5)
    value = 1 - value;
  value = value / 0.5;
  igl::parula(value, color[0], color[1], color[2]);
  
  P1.conservativeResize(P1.rows()+sl_state.start_point.rows(),3);
  P2.conservativeResize(P2.rows()+sl_state.end_point.rows(),3);
  P1.block(P1.rows()-sl_state.start_point.rows(),0,sl_state.start_point.rows(),3)=sl_state.start_point;
  P2.block(P2.rows()-sl_state.end_point.rows(),0,sl_state.end_point.rows(),3)=sl_state.end_point;
  
  directional_viewer->set_streamlines(P1, P2, color.replicate(P2.rows(),1));
  
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
  directional::power_field(V, F, Eigen::VectorXi(),  Eigen::MatrixXd() , 3, powerField);
  
  // Convert it to raw field
  directional::power_to_raw(V,F,powerField,3,rawField, true);
  
  directional::streamlines_init(V, F, rawField, sl_data, sl_state);
  
  //triangle mesh
  viewer.set_mesh(V,F);
  viewer.data().show_lines=false;
  
  // Viewer Settings
  viewer.callback_pre_draw = &pre_draw;
  viewer.callback_key_down = &key_down;
  viewer.core().is_animating = false;
  viewer.core().animation_max_fps = 30.;
  
  // Draw initial seeds on sample points
  directional::StreamlineState sl_state0;
  sl_state0 = sl_state;
  directional::streamlines_next(V, F, sl_data, sl_state0);
  Eigen::MatrixXd v = sl_state0.end_point - sl_state0.start_point;
  v.rowwise().normalize();
  
  //streamline mesh
  P1=sl_state0.start_point;
  P2=sl_state0.start_point + 0.0005 * v;
  viewer.set_streamlines(P1, P2, Eigen::MatrixXd::Constant(sl_state0.start_point.rows(),3,1.0));
  
  
  cout << "Press [space] to toggle animation" << endl;
  viewer.launch();
}
