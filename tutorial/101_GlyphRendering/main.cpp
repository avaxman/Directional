#include <directional/directional_viewer.h>
#include <directional/read_raw_field.h>
#include <directional/read_singularities.h>

int N;
Eigen::MatrixXi F, EV, FE, EF;
Eigen::MatrixXd V, rawField;
Eigen::VectorXi singVertices, singIndices;
directional::DirectionalViewer viewer;

bool key_down(igl::opengl::glfw::Viewer& iglViewer, int key, int modifiers)
{
  igl::opengl::glfw::Viewer* iglViewerPointer=&iglViewer;
  directional::DirectionalViewer* directional_viewer = static_cast<directional::DirectionalViewer*>(iglViewerPointer);
  switch (key)
  {
    case '1': directional_viewer->toggle_field(); return true;
    case '2': directional_viewer->toggle_singularities(); return true;
    default: return false;
  }
  
}

int main()
{
  std::cout <<"<space bar>  Show/hide Singularities" << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/bumpy.off", V, F);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/bumpy.rawfield", N, rawField);
  directional::read_singularities(TUTORIAL_SHARED_PATH "/bumpy.sings", N, singVertices, singIndices);
  
  directional::DirectionalViewer viewer;
  
  viewer.set_mesh(V,F);
  viewer.set_field(rawField);
  viewer.set_singularities(N, singVertices, singIndices);

  viewer.callback_key_down = &key_down;
  viewer.launch();
}

