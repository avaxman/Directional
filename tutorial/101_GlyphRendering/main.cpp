#include <directional/directional_viewer.h>
#include <directional/readOFF.h>
#include <directional/read_raw_field.h>
#include <directional/read_singularities.h>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>

int N;
directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField field;
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
  
  directional::readOFF(TUTORIAL_SHARED_PATH "/bumpy.off",mesh);
  ftb.init(mesh);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/bumpy.rawfield", ftb, N, field);
  directional::read_singularities(TUTORIAL_SHARED_PATH "/bumpy.sings", field);
  directional::DirectionalViewer viewer;
  
  viewer.set_mesh(mesh);
  viewer.set_field(field);
  viewer.toggle_mesh_edges(false);
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

