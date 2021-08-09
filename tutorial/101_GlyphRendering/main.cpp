#include <directional/directional_viewer.h>
#include <directional/read_raw_field.h>
#include <directional/read_singularities.h>

int N;
Eigen::MatrixXi F, EV, FE, EF;
Eigen::MatrixXd V, rawField;
Eigen::VectorXi singVertices, singIndices;
directional::DirectionalViewer viewer;
bool showField=true, showSingularities=true;

bool key_down(igl::opengl::glfw::Viewer& iglViewer, int key, int modifiers)
{
  igl::opengl::glfw::Viewer* iglViewerPointer=&iglViewer;
  directional::DirectionalViewer* directional_viewer = static_cast<directional::DirectionalViewer*>(iglViewerPointer);
  switch (key)
  {
    case '1': showField=!showField; directional_viewer->toggle_field(showField); break;
    case '2': showSingularities=!showSingularities; directional_viewer->toggle_singularities(showSingularities); break;;
    default: return false;
  }
  return true;
}

int main()
{
  std::cout <<"1    Show/hide Field" << std::endl;
  std::cout <<"2    Show/hide Singularities" << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/bumpy.off", V, F);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/bumpy.rawfield", N, rawField);
  directional::read_singularities(TUTORIAL_SHARED_PATH "/bumpy.sings", N, singVertices, singIndices);
  
  directional::DirectionalViewer viewer;
  
  viewer.set_mesh(V,F);
  viewer.set_field(rawField);
  viewer.set_singularities(singVertices, singIndices);
  viewer.toggle_mesh_edges(false);
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

