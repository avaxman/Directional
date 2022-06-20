#include <iostream>
#include <Eigen/Core>
#include <igl/read_triangle_mesh.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/directional_viewer.h>


int currF=0, N;
Eigen::MatrixXi F, FField, FSings;
Eigen::MatrixXd V, VField, VSings;
Eigen::MatrixXd CMesh, CField, CSings;
Eigen::MatrixXd rawField, barycenters;
Eigen::VectorXd effort;
Eigen::VectorXi matching;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi singIndices, singVertices;

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
    otherFaces(i)=(EF(FE(currF,i),0)==currF ? EF(FE(currF,i),1) : EF(FE(currF,i),0));
    zeroInFace(i)=(EF(FE(currF,i),0)==currF ? matching(FE(currF,i)) : -matching(FE(currF,i)));
  }
  
  Eigen::MatrixXd glyphColors=directional::DirectionalViewer::default_glyph_color().replicate(F.rows(),N);
  glyphColors.row(currF)=directional::DirectionalViewer::indexed_glyph_colors(rawField.row(currF), false);
  for (int i=0;i<N;i++)
    for (int j=0;j<3;j++)
      glyphColors.block(otherFaces(j),3*((i+zeroInFace(j)+N)%N),1,3)<<glyphColors.row(currF).segment(3*i,3);
  
  viewer.set_field(rawField, glyphColors);
  
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
                               directional_viewer->core().proj, directional_viewer->core().viewport, V, F, fid, baryInFace))
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
  igl::readOBJ(TUTORIAL_SHARED_PATH "/lilium.obj", V, F);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/lilium.rawfield", N, rawField);
  igl::edge_topology(V, F, EV, FE, EF);
  igl::barycenter(V, F, barycenters);
  
  directional::principal_matching(V, F,EV, EF, FE, rawField, matching, effort, singVertices, singIndices);

  //triangle mesh setup
  viewer.set_mesh(V, F);
  update_triangle_mesh();
  update_raw_field_mesh();
  
  //singularity mesh
  viewer.set_singularities(singVertices, singIndices);
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}

