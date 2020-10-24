#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/edge_topology.h>
#include <directional/seam_lines.h>
#include <directional/glyph_lines_raw.h>
#include <directional/read_raw_field.h>
#include <directional/principal_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/singularity_spheres.h>
#include <directional/combing.h>
#include <directional/line_cylinders.h>
#include <directional/directional_viewer.h>


int currF=0, N;
Eigen::MatrixXi F;
Eigen::MatrixXd V;
Eigen::MatrixXd rawField, combedField, barycenters;
directional::DirectionalViewer viewer;
Eigen::VectorXi matching, combedMatching;
Eigen::VectorXd effort, combedEffort;
Eigen::MatrixXi EV, FE, EF;
Eigen::VectorXi singIndices, singVertices;
bool showCombed=false;
bool showSingularities=true;

void update_raw_field_mesh()
{
  Eigen::MatrixXd currField = (showCombed ? combedField : rawField);
  viewer.set_field(currField,directional::DirectionalViewer::indexed_glyph_colors(currField, false));
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
  return true;
}


int main()
{
  std::cout <<
  "  1        Toggle raw field/Combed field" << std::endl <<
  igl::readOBJ(TUTORIAL_SHARED_PATH "/lilium.obj", V, F);
  directional::read_raw_field(TUTORIAL_SHARED_PATH "/lilium.rawfield", N, rawField);
  igl::edge_topology(V, F, EV, FE, EF);
  igl::barycenter(V, F, barycenters);
  
  //computing
  directional::principal_matching(V, F,EV, EF, FE, rawField, matching, effort,singVertices, singIndices);
  directional::combing(V,F, EV, EF, FE, rawField, matching, combedField);
  directional::principal_matching(V, F,EV, EF, FE, combedField, combedMatching, combedEffort,singVertices, singIndices);
  
  //Mesh setup
  viewer.set_mesh(V, F);
  viewer.toggle_mesh_edges(false);
  update_raw_field_mesh();
  viewer.set_singularities(N, singVertices, singIndices);
  viewer.set_seams(EV, combedMatching);  //TODO: allow to define seams in several ways
  
  viewer.callback_key_down = &key_down;
  viewer.launch();
}


