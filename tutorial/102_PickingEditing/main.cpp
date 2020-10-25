#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <directional/directional_viewer.h>
#include <directional/read_raw_field.h>

int currF, currVec, N;
Eigen::MatrixXi F;
Eigen::MatrixXd V, barycenters;
Eigen::MatrixXd rawField;
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
      //directionalViewer->set_selected_vector(rawField, currF, currVec);
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
                               viewer.core().proj, viewer.core().viewport, V, F, fid, baryInFace))
  {
    
    //choosing face
    if ((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Left){
      currF=fid;
      Eigen::VectorXi selectedFaces(1); selectedFaces(0)=currF;
      viewer.set_selected_faces(selectedFaces);
      //directionalViewer->set_selected_vector(rawField, currF, currVec);
      return true;
    }
    //choosing face
    if (((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Right)&&(fid==currF)){
      // Calculate direction from the center of the face to the mouse
      Eigen::RowVector3d newVec =(V.row(F(fid, 0)) * baryInFace(0) +
                                  V.row(F(fid, 1)) * baryInFace(1) +
                                  V.row(F(fid, 2)) * baryInFace(2) - barycenters.row(fid)).normalized();
      
      rawField.block(currF, currVec*3, 1,3)=newVec;
      //directionalViewer->set_selected_vector(rawField, currF, currVec);
      return true;
      
    }
  }
  return false;
};

int main()
{
  igl::readOBJ(TUTORIAL_SHARED_PATH "/torus.obj", V, F);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/torus.rawfield", N, rawField);
  
  std::cout <<
  "  1                Choose vector in current face." << std::endl <<
  "  0+Left button    Choose face" << std::endl <<
  "  0+Right button   Edit vector in current face" << std::endl;
  
  igl::barycenter(V, F, barycenters);
  
  viewer.set_mesh(V,F);
  viewer.set_field(rawField);
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
