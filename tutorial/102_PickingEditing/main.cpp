#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/unproject_onto_mesh.h>
#include <directional/directional_viewer.h>
#include <directional/read_raw_field.h>
#include <directional/TriMesh.h>
#include <directional/FaceField.h>

int currF, currVec, N;
directional::TriMesh mesh;
directional::FaceField field;
directional::DirectionalViewer viewer;

//User input variables
bool zeroPressed = false;


bool key_up(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '0': zeroPressed=false; break;
  }
  return true;
}

// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& iglViewer, int key, int modifiers)
{
  igl::opengl::glfw::Viewer* iglViewerPointer=&iglViewer;
  directional::DirectionalViewer* directionalViewer = static_cast<directional::DirectionalViewer*>(iglViewerPointer);
  switch (key)
  {
      // Select vector
    case '0': zeroPressed=true; break;
    case '1':
      currVec = (currVec+1)%N;
      directionalViewer->set_selected_vector(currF, currVec);
      break;
  }
  return true;
}

//Select vertices using the mouse
bool mouse_down(igl::opengl::glfw::Viewer& iglViewer, int button, int modifiers)
{
  igl::opengl::glfw::Viewer* iglViewerPointer=&iglViewer;
  directional::DirectionalViewer* directionalViewer = static_cast<directional::DirectionalViewer*>(iglViewerPointer);
  if (!zeroPressed)
    return false;
  int fid;
  Eigen::Vector3d baryInFace;
  
  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core().viewport(3) - viewer.current_mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core().view,
                               viewer.core().proj, viewer.core().viewport, mesh.V, mesh.F, fid, baryInFace))
  {
    
    //choosing face
    if ((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Left){
      currF=fid;
      Eigen::VectorXi selectedFaces(1); selectedFaces(0)=currF;
      directionalViewer->set_selected_faces(selectedFaces);
      directionalViewer->set_selected_vector(currF, currVec);
      return true;
    }
    //moving vector
    if (((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Right)&&(fid==currF)){
      // Calculate direction from the center of the face to the mouse
      Eigen::RowVector3d newVec =(mesh.V.row(mesh.F(fid, 0)) * baryInFace(0) +
                                  mesh.V.row(mesh.F(fid, 1)) * baryInFace(1) +
                                  mesh.V.row(mesh.F(fid, 2)) * baryInFace(2) - mesh.barycenters.row(fid)).normalized();
      
      field.extField.block(currF, currVec*3, 1,3)=newVec;
      directionalViewer->set_field(field);
      directionalViewer->set_selected_vector(currF, currVec);
      return true;
      
    }
  }
  return false;
};

int main()
{
  Eigen::MatrixXd V,rawField;
  Eigen::MatrixXi F;
  igl::readOBJ(TUTORIAL_SHARED_PATH "/torus.obj", V, F);
  mesh.set_mesh(V,F);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/torus.rawfield", N, rawField);
  field.set_field(rawField,mesh);
  
  std::cout <<
  "  1                Choose vector in current face." << std::endl <<
  "  0+Left button    Choose face" << std::endl <<
  "  0+Right button   Edit vector in current face" << std::endl;
  
  viewer.set_mesh(mesh);
  viewer.set_field(field);
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
