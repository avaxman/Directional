#include <iostream>
#include <Eigen/Core>
#include <igl/unproject_onto_mesh.h>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/readOBJ.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/directional_viewer.h>

int currF=0, N;
directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField field;
directional::DirectionalViewer viewer;

bool zeroPressed=false;

void update_triangle_mesh()
{
  Eigen::VectorXi selectedFaces(1); selectedFaces(0)=currF;
  viewer.set_selected_faces(selectedFaces);
}

void update_raw_field_mesh()
{
  //configuring just part of the faces to have combed coloring
  Eigen::Vector3i otherFaces;
  Eigen::Vector3i zeroInFace;
  for (int i=0;i<3;i++){
    otherFaces(i)=(mesh.EF(mesh.FE(currF,i),0)==currF ? mesh.EF(mesh.FE(currF,i),1) : mesh.EF(mesh.FE(currF,i),0));
    zeroInFace(i)=(mesh.EF(mesh.FE(currF,i),0)==currF ? field.matching(mesh.FE(currF,i)) : -field.matching(mesh.FE(currF,i)));
  }
  
  Eigen::MatrixXd glyphColors=directional::DirectionalViewer::default_glyph_color().replicate(mesh.F.rows(),N);
  glyphColors.row(currF)=directional::DirectionalViewer::indexed_glyph_colors(field.extField.row(currF), false);
  for (int i=0;i<N;i++)
    for (int j=0;j<3;j++)
      glyphColors.block(otherFaces(j),3*((i+zeroInFace(j)+N)%N),1,3)<<glyphColors.row(currF).segment(3*i,3);
  
  viewer.set_field(field, glyphColors);
  
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
bool key_down(igl::opengl::glfw::Viewer& iglViewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '0': zeroPressed=true; break;
  }
  return true;
}

//Select vertices using the mouse
bool mouse_down(igl::opengl::glfw::Viewer& iglViewer, int button, int modifiers)
{
  //TODO: encapsulated picking
  igl::opengl::glfw::Viewer* iglViewerPointer=&iglViewer;
   directional::DirectionalViewer* directional_viewer = static_cast<directional::DirectionalViewer*>(iglViewerPointer);
  if (!zeroPressed)
    return false;
  int fid;
  Eigen::Vector3d baryInFace;
  
  // Cast a ray in the view direction starting from the mouse position
  double x = directional_viewer->current_mouse_x;
  double y = directional_viewer->core().viewport(3) - directional_viewer->current_mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), directional_viewer->core().view,
                               directional_viewer->core().proj, directional_viewer->core().viewport, mesh.V, mesh.F, fid, baryInFace))
  {
    
    //choosing face
    if ((igl::opengl::glfw::Viewer::MouseButton)button==igl::opengl::glfw::Viewer::MouseButton::Left){
      currF=fid;
      update_triangle_mesh();
      update_raw_field_mesh();
      return true;
    }
  }
  return false;
};

int main()
{
  std::cout <<
  "  0+Left button    Choose face" << std::endl <<
  directional::readOBJ(TUTORIAL_SHARED_PATH "/lilium.obj", mesh);
  ftb.init(mesh);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/lilium.rawfield", ftb, N, field);
  
  directional::principal_matching(field);

  //triangle mesh setup
  viewer.set_mesh(mesh);
  viewer.set_field(field);
  update_triangle_mesh();
  update_raw_field_mesh();

  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}

