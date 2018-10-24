#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <directional/visualization_schemes.h>
#include <directional/glyph_lines_raw.h>
#include <directional/read_raw_field.h>

int currF, currVec, N;
Eigen::MatrixXi FMesh, FField;
Eigen::MatrixXd VMesh, VField,  barycenters;
Eigen::MatrixXd CMesh, CField;
Eigen::MatrixXd rawField;
igl::opengl::glfw::Viewer viewer;

//User input variables
bool zeroPressed = false;

void update_triangle_mesh()
{
  
  CMesh=directional::default_mesh_color().replicate(FMesh.rows(),1);
  CMesh.row(currF)=directional::selected_face_color();
  viewer.data_list[0].set_colors(CMesh);
}

void update_raw_field_mesh()
{
  Eigen::MatrixXd glyphColors=directional::default_glyph_color().replicate(FMesh.rows(),N);
  glyphColors.row(currF)=directional::selected_face_glyph_color().replicate(1,N);
  glyphColors.block(currF,3*currVec,1,3)=directional::selected_vector_glyph_color();
  
  directional::glyph_lines_raw(VMesh, FMesh, rawField, glyphColors, VField, FField, CField);
  
  viewer.data_list[1].set_mesh(VField, FField);
  viewer.data_list[1].set_colors(CField);
  viewer.data_list[1].show_lines=false;
}


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
bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '0': zeroPressed=true; break;
    case '1':
    case '2':
    case '3':
    case '4':
      currVec = key - '1';
      update_raw_field_mesh();
      break;
  }
  return true;
}

//Select vertices using the mouse
bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifiers)
{
  if (!zeroPressed)
    return false;
  int fid;
  Eigen::Vector3d bc;
  
  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view,
                               viewer.core.proj, viewer.core.viewport, VMesh, FMesh, fid, bc))
  {
    
    //choosing face
    if ((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Left){
      currF=fid;
      update_triangle_mesh();
      update_raw_field_mesh();
      return true;
    }
    //choosing face
    if (((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Right)&&(fid==currF)){
      // Calculate direction from the center of the face to the mouse
      Eigen::RowVector3d newVec =(VMesh.row(FMesh(fid, 0)) * bc(0) +
                                  VMesh.row(FMesh(fid, 1)) * bc(1) +
                                  VMesh.row(FMesh(fid, 2)) * bc(2) - barycenters.row(fid)).normalized();
      
      rawField.block(currF, currVec*3, 1,3)=newVec;
      update_raw_field_mesh();
      return true;
      
    }
  }
  return false;
};

int main()
{
  igl::readOBJ(TUTORIAL_SHARED_PATH "/torus.obj", VMesh, FMesh);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/torus.rawfield", N, rawField);
  
  std::cout <<
  "  1-"<< N <<"  Choose vector in current face." << std::endl <<
  "  0+Left button    Choose face" << std::endl <<
  "  0+Right button   Edit vector in current face" << std::endl;
  
  igl::barycenter(VMesh, FMesh, barycenters);
  
  //triangle mesh setup
  viewer.data().set_face_based(true);
  viewer.data().set_mesh(VMesh, FMesh);
  update_triangle_mesh();
  //apending and updating raw field mesh
  viewer.append_mesh();
  update_raw_field_mesh();
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
