#include <igl/opengl/glfw/Viewer.h>
#include <directional/visualization_schemes.h>
#include <directional/read_raw_field.h>
#include <directional/read_singularities.h>
#include <directional/singularity_spheres.h>
#include <directional/glyph_lines_raw.h>

int N;
Eigen::MatrixXi FMesh, FField, FSings;
Eigen::MatrixXi EV, FE, EF;
Eigen::MatrixXd VMesh, VField, VSings;
Eigen::MatrixXd CMesh, CField, CSings;
Eigen::MatrixXd rawField;
Eigen::VectorXi singVertices, singIndices;
igl::opengl::glfw::Viewer viewer;

Eigen::MatrixXd positiveIndexColors, negativeIndexColors;
Eigen::RowVector3d rawGlyphColor;
bool drawSingularities=false;


void create_meshes()
{
  directional::glyph_lines_raw(VMesh, FMesh, rawField, directional::default_glyph_color(), VField, FField, CField);
  directional::singularity_spheres(VMesh, FMesh, N, singVertices, singIndices, VSings, FSings, CSings);
  
  //triangle mesh
  viewer.data().set_mesh(VMesh, FMesh);
  viewer.data().set_colors(directional::default_mesh_color());
  viewer.data().set_face_based(true);
  viewer.data().show_lines=true;
  
  //Raw field mesh
  viewer.append_mesh();
  viewer.data().set_mesh(VField, FField);
  viewer.data().set_colors(CField);
  viewer.data().set_face_based(true);
  viewer.data().show_lines=false;
  
  //Singularities mesh
  viewer.append_mesh();
  viewer.data().set_mesh(VSings, FSings);
  viewer.data().set_colors(CSings);
  viewer.data().set_face_based(true);
  viewer.data().show_lines=false;
  viewer.data().show_faces=false;
  
  viewer.selected_data_index=0;  //for all generic libigl commands.
  
}


bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  
  switch (key)
  {
    case GLFW_KEY_SPACE: drawSingularities=!drawSingularities; viewer.data_list[2].show_faces=!viewer.data_list[2].show_faces; /*update_mesh();*/ break;
    default: return false;
  }
  return true;
}



int main()
{
  std::cout <<"<space bar>  Show/hide Singularities" << std::endl;
  
  igl::readOFF(TUTORIAL_SHARED_PATH "/bumpy.off", VMesh, FMesh);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/bumpy.rawfield", N, rawField);
  directional::read_singularities(TUTORIAL_SHARED_PATH "/bumpy.sings", N, singVertices, singIndices);
  
  // Set colors for Singularities
  
  
  create_meshes();
  viewer.callback_key_down = &key_down;
  viewer.launch();
}

