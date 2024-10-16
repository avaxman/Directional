#include <iostream>
#include <Eigen/Core>
#include <igl/unproject_onto_mesh.h>
#include <directional/TriMesh.h>
#include <directional/IntrinsicFaceTangentBundle.h>
#include <directional/CartesianField.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/combing.h>
#include <directional/directional_viewer.h>
#include <directional/readOBJ.h>

int currF=0, N;
directional::TriMesh mesh;
directional::IntrinsicFaceTangentBundle ftb;
directional::CartesianField rawField, combedField;
directional::DirectionalViewer viewer;

bool showCombed=false;
bool showSingularities=true;

void update_raw_field_mesh()
{
  viewer.set_field((showCombed ? combedField : rawField),directional::DirectionalViewer::indexed_glyph_colors(rawField.extField, false));
  viewer.toggle_seams(showCombed);
}


// Handle keyboard input
bool key_down(igl::opengl::glfw::Viewer& iglViewer, int key, int modifiers)
{
  switch (key)
  {
      // Select vector
    case '1': showCombed = !showCombed; update_raw_field_mesh(); break;
  }
  if (showCombed)
    std::cout<<"Showing combed field"<<std::endl;
  else
    std::cout<<"Showing raw field"<<std::endl;
  return true;
}


int main()
{
  std::cout <<
  "  1        Toggle raw field/Combed field" << std::endl <<
  directional::readOBJ(TUTORIAL_SHARED_PATH "/lilium.obj", mesh);
  ftb.init(mesh);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/lilium.rawfield", ftb, N, rawField);
  std::cout<<"Showing raw field"<<std::endl;

  //computing the principal matching and then combing the field to trivialize the matching anywhere but a set of seams
  directional::principal_matching(rawField);
  directional::combing(rawField, combedField);

  viewer.set_mesh(mesh);
  viewer.toggle_mesh_edges(false);
  update_raw_field_mesh();
  viewer.set_seams(combedField.matching);
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


