#include <iostream>
#include <Eigen/Core>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/read_triangle_mesh.h>
#include <igl/per_face_normals.h>
#include <igl/triangle_triangle_adjacency.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/boundary_loop.h>
#include <igl/edge_topology.h>
#include <directional/visualization_schemes.h>
#include <directional/power_field.h>
#include <directional/power_to_representative.h>
#include <directional/power_to_raw.h>
#include <directional/representative_to_raw.h>
#include <directional/principal_matching.h>
#include <directional/effort_to_indices.h>
#include <directional/glyph_lines_raw.h>
#include <directional/singularity_spheres.h>
#include <directional/write_raw_field.h>


Eigen::VectorXi b, matching, singVertices, singIndices;
Eigen::VectorXd effort;
Eigen::MatrixXi FMesh, FField, FSings;
Eigen::MatrixXi EV, EF, FE;
Eigen::MatrixXd VMesh, VField, VSings;
Eigen::MatrixXd CMesh, CField, CSings;
Eigen::MatrixXd rawField,representative, bc, barycenters;
Eigen::MatrixXcd powerField;
igl::opengl::glfw::Viewer viewer;

int N = 4;
bool normalized = false;
bool zeroPressed = false;

void update_triangle_mesh()
{
  
  Eigen::MatrixXd CMesh=directional::default_mesh_color().replicate(FMesh.rows(),1);
  for (int i = 0; i < b.rows(); i++)
    CMesh.row(b(i)) = directional::selected_face_color();
  
  viewer.data_list[0].set_colors(CMesh);
}

void recompute_field()
{
  Eigen::VectorXi bcSoft;
  Eigen::MatrixXd wSoft;
  Eigen::MatrixXd bSoft;
  directional::power_field(VMesh, FMesh, b, bc, bcSoft, wSoft, bSoft, N, powerField);
}

void update_raw_field_mesh()
{
  directional::power_to_representative(VMesh, FMesh, powerField, N, representative);
  if (normalized)
    representative.rowwise().normalize();
  
  directional::representative_to_raw(VMesh,FMesh,representative, N, rawField);
  directional::principal_matching(VMesh, FMesh, EV, EF, FE, rawField, matching, effort);
  directional::effort_to_indices(VMesh,FMesh,EV, EF, effort,matching, N, singVertices, singIndices);
  
  directional::singularity_spheres(VMesh, FMesh, N, singVertices, singIndices, VSings, FSings, CSings);
  
  Eigen::MatrixXd glyphColors=directional::default_glyph_color().replicate(FMesh.rows(),N);
  if (b.rows()!=0){
    glyphColors.row(b(b.rows()-1))=directional::selected_face_glyph_color().replicate(1,N);
    //glyphColors.block(b(b.rows()-1),0,1,3)=directional::selected_vector_glyph_color();
  }
  
  directional::glyph_lines_raw(VMesh, FMesh, rawField, glyphColors,  VField, FField, CField);
  
  if (viewer.data_list.size()<2){
    
    //apending and updating raw field mesh
    viewer.append_mesh();
    viewer.data_list[1].show_faces = true;
    viewer.data_list[1].show_lines = false;
    
    viewer.append_mesh();
    viewer.data_list[2].show_faces = true;
    viewer.data_list[2].show_lines = false;
    
    viewer.selected_data_index = 0;
  }
  
  viewer.data_list[1].clear();
  viewer.data_list[1].set_mesh(VField, FField);
  viewer.data_list[1].set_colors(CField);
  
  viewer.data_list[2].clear();
  viewer.data_list[2].set_mesh(VSings, FSings);
  viewer.data_list[2].set_colors(CSings);
  
  
}

bool key_up(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
    case '0': zeroPressed=false; break;
  }
  return true;
}


bool key_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  switch (key)
  {
      // Toggle field drawing for easier rotation
      
    case '0': zeroPressed=true; break;
      
      // Reset the constraints
    case 'R':
      b.resize(0);
      bc.resize(0, 3);
      recompute_field();
      update_triangle_mesh();
      update_raw_field_mesh();
      break;
      
      // Toggle normalization
    case 'N':
      normalized = !normalized;
      update_triangle_mesh();
      update_raw_field_mesh();
      break;
      
      
    case 'W':
      if (directional::write_raw_field(TUTORIAL_SHARED_PATH "/rocker-arm2500.rawfield", rawField))
        std::cout << "Saved raw field" << std::endl;
      else
        std::cout << "Unable to save raw field. Error: " << errno << std::endl;
      
  }
  
  return true;
}

//Select vertices using the mouse
bool mouse_down(igl::opengl::glfw::Viewer& viewer, int key, int modifiers)
{
  if (!zeroPressed)
    return false;
  
  int fid;
  Eigen::Vector3d baryInFace;
  
  // Cast a ray in the view direction starting from the mouse position
  double x = viewer.current_mouse_x;
  double y = viewer.core.viewport(3) - viewer.current_mouse_y;
  if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y), viewer.core.view,
                               viewer.core.proj, viewer.core.viewport, VMesh, FMesh, fid, baryInFace))
  {
   
    int i;
    for (i = 0; i < b.rows(); i++)
      if (b(i) == fid)
        break;
    if (i == b.rows())
    {
      b.conservativeResize(b.rows() + 1);
      b(i) = fid;
      bc.conservativeResize(bc.rows() + 1, 3);
    }
    
    // Compute direction from the center of the face to the mouse
    bc.row(i) =(VMesh.row(FMesh(fid, 0)) * baryInFace(0) +
                VMesh.row(FMesh(fid, 1)) * baryInFace(1) +
                VMesh.row(FMesh(fid, 2)) * baryInFace(2) - barycenters.row(fid)).normalized();
    recompute_field();
    update_triangle_mesh();
    update_raw_field_mesh();
    return true;
  }
  return false;
};

int main()
{
  
  std::cout <<
  "  R        Reset the constraints" << std::endl <<
  "  N        Toggle field normalization" << std::endl <<
  "  W        Save raw field" << std::endl <<
  "  0+L-bttn Place constraint pointing from the center of face to the cursor" << std::endl;
  
  // Load mesh
  igl::readOBJ(TUTORIAL_SHARED_PATH "/rocker-arm2500.obj", VMesh, FMesh);
  igl::edge_topology(VMesh, FMesh, EV,FE,EF);
  igl::barycenter(VMesh, FMesh, barycenters);
  
  b.resize(0);
  bc.resize(0, 3);
  
  //triangle mesh setup
  viewer.data_list[0].set_mesh(VMesh, FMesh);
  viewer.data_list[0].set_colors(directional::default_mesh_color());
  
  viewer.selected_data_index = 0;
  recompute_field();
  update_raw_field_mesh();
  
  viewer.callback_key_down = &key_down;
  viewer.callback_key_up = &key_up;
  viewer.callback_mouse_down = &mouse_down;
  viewer.launch();
}
